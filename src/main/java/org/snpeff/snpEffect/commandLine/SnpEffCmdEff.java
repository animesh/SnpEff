package org.snpeff.snpEffect.commandLine;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.snpeff.SnpEff;
import org.snpeff.annotate.AnnotateVcf;
import org.snpeff.fileIterator.VcfFileIterator;
import org.snpeff.filter.VariantEffectFilter;
import org.snpeff.interval.tree.IntervalForest;
import org.snpeff.outputFormatter.VcfOutputFormatter;
import org.snpeff.snpEffect.EffectType;
import org.snpeff.snpEffect.SnpEffectPredictor;
import org.snpeff.snpEffect.VcfAnnotator;
import org.snpeff.stats.CountByType;
import org.snpeff.stats.VariantEffectStats;
import org.snpeff.stats.VcfStats;
import org.snpeff.util.Gpr;
import org.snpeff.vcf.EffFormatVersion;
import org.snpeff.vcf.PedigreeEnrty;
import org.snpeff.vcf.VcfEntry;

/**
 * Command line program: Predict variant effects
 *
 * @author Pablo Cingolani
 */
public class SnpEffCmdEff extends SnpEff implements VcfAnnotator {

	public static final String DEFAULT_SUMMARY_CSV_FILE = "snpEff_summary.csv";
	public static final String DEFAULT_SUMMARY_GENES_FILE = "snpEff_genes.txt";
	public static final String DEFAULT_SUMMARY_HTML_FILE = "snpEff_summary.html";

	AnnotateVcf annotateVcf;
	boolean anyCancerSample;
	boolean cancer = false; // Perform cancer comparisons
	String cancerSamples = null;
	boolean chromoPlots = true; // Create mutations by chromosome plots?
	String chrStr = "";
	String commandLineStr, commandLineStrReport;
	boolean createSummaryCsv = false;
	boolean createSummaryHtml = true;
	CountByType errByType, warnByType;
	String fastaProt = null;
	IntervalForest filterIntervals; // Filter only variants that match these intervals
	EffFormatVersion formatVersion = EffFormatVersion.DEFAULT_FORMAT_VERSION;
	boolean gatk = false; // Use GATK compatibility mode
	String inputFile = ""; // Input file
	ArrayList<String> inputFiles;
	boolean lossOfFunction = true; // Create loss of function LOF tag?
	List<PedigreeEnrty> pedigree;
	SnpEffectPredictor snpEffectPredictor;
	String summaryFileCsv; // HTML Summary file name
	String summaryFileHtml; // CSV Summary file name
	String summaryGenesFile; // Gene table file
	boolean useGeneId = false; // Use gene ID instead of gene name (VCF output)
	boolean useLocalTemplate = false; // Use template from 'local' file instead of 'jar' (this is only used for development and debugging)
	boolean useOicr = false; // Use OICR tag
	boolean useSequenceOntology = true; // Use Sequence Ontology terms
	VariantEffectFilter variantEffectResutFilter; // Filter prediction results
	List<VcfEntry> vcfEntriesDebug = null; // Use for debugging or testing (in some test-cases)
	VcfOutputFormatter vcfOutputFormatter = null;
	VcfStats vcfStats;

	public SnpEffCmdEff() {
		super();
		chrStr = ""; // Default: Don't show 'chr' before chromosome
		inputFile = ""; // variant input file
		variantEffectResutFilter = new VariantEffectFilter(); // Filter prediction results
		summaryFileHtml = DEFAULT_SUMMARY_HTML_FILE;
		summaryFileCsv = DEFAULT_SUMMARY_CSV_FILE;
		summaryGenesFile = DEFAULT_SUMMARY_GENES_FILE;
		annotateVcf = new AnnotateVcf();
	}

	@Override
	public boolean addHeaders(VcfFileIterator vcfFile) {
		return annotateVcf.addHeaders(vcfFile);
	}

	@Override
	public boolean annotate(VcfEntry vcfEntry) {
		return annotateVcf.annotate(vcfEntry);
	}

	@Override
	public boolean annotateFinish(VcfFileIterator vcfFile) {
		return annotateVcf.annotateFinish(vcfFile);
	}

	@Override
	public boolean annotateInit(VcfFileIterator vcfFile) {
		// Set all parsed parameters into annotateVcf
		copyValuesToAnnotateVcf();
		return annotateVcf.annotateInit(vcfFile);
	}

	void copyValuesToAnnotateVcf() {
		commandLineStr = commandLineStr(false);
		commandLineStrReport = commandLineStr(createSummaryCsv ? false : true);

		// Set all parsed parameters into annotateVcf
		annotateVcf.set(this);
	}

	public VariantEffectStats getChangeEffectResutStats() {
		return annotateVcf.getChangeEffectResutStats();
	}

	public int getTotalErrs() {
		return annotateVcf.getTotalErrs();
	}

	/**
	 * Parse command line arguments
	 */
	@Override
	public void parseArgs(String[] args) {
		if (args == null) return;

		boolean isFileList = false;
		this.args = args;

		for (int i = 0; i < args.length; i++) {
			String arg = args[i];

			// Is it a command line option?
			// Note: Generic options (such as config, verbose, debug, quiet, etc.) are parsed by SnpEff class
			//---
			if (isOpt(arg)) {
				if (arg.equalsIgnoreCase("-fileList")) isFileList = true;
				else {
					arg = arg.toLowerCase();

					switch (arg) {
					//---
					// Output options
					//---
					case "-chr":
						chrStr = args[++i];
						break;

					case "-csvstats":
						createSummaryCsv = true; // Create a CSV formatted summary file.
						if ((i + 1) < args.length) {
							summaryFileCsv = args[++i];
							String base = Gpr.baseName(summaryFileCsv, ".csv");
							String dir = Gpr.dirName(summaryFileCsv);
							summaryGenesFile = (dir != null ? dir + "/" : "") + base + ".genes.txt";
						} else usage("Missing parameter: CSV stats file name ");
						break;

					case "-fastaprot":
						if ((i + 1) < args.length) fastaProt = args[++i]; // Output protein sequences in fasta files
						else usage("Missing -cancerSamples argument");
						break;

					case "-gatk": // GATK compatibility mode
						gatk = true;
						useSequenceOntology = false;
						hgvs = false;
						nextProt = false;
						motif = false;
						// GATK doesn't support SPLICE_REGION at the moment.
						// Set parameters to zero so that splcie regions are not created.
						spliceRegionExonSize = spliceRegionIntronMin = spliceRegionIntronMax = 0;
						break;

					case "-nochromoplots":
						chromoPlots = false;
						break;

					case "-nostats":
						createSummaryHtml = createSummaryCsv = false;
						break;

					case "-o": // Output format
						usage("Command line option '-o' no longer supported: Only VCF format can be used for input and output.");
						break;

					case "-s":
					case "-stats":
					case "-htmlstats":
						createSummaryHtml = true;
						if ((i + 1) < args.length) {
							summaryFileHtml = args[++i];
							String base = Gpr.baseName(summaryFileHtml, ".html");
							String dir = Gpr.dirName(summaryFileHtml);
							summaryGenesFile = (dir != null ? dir + "/" : "") + base + ".genes.txt";
						} else usage("Missing parameter: HTML stats file name ");
						break;

					case "-uselocaltemplate": // Undocumented option (only used for development & debugging)
						useLocalTemplate = true;
						break;

					//---
					// Annotation options
					//---
					case "-cancer":
						cancer = true; // Perform cancer comparisons
						break;

					case "-cancersamples":
						if ((i + 1) < args.length) cancerSamples = args[++i]; // Read cancer samples from TXT files
						else usage("Missing -cancerSamples argument");
						break;

					case "-classic":
						useSequenceOntology = false;
						formatVersion = EffFormatVersion.FORMAT_EFF_4;
						hgvs = hgvsForce;
						break;

					case "-formateff":
						formatVersion = EffFormatVersion.FORMAT_EFF_4;
						break;

					case "-geneid":
						useGeneId = true; // Use gene ID instead of gene name
						break;

					case "-lof":
						lossOfFunction = true; // Add LOF tag
						break;

					case "-nohgvs":
						hgvs = false; // Do not use HGVS notation
						hgvsShift = false;
						break;

					case "-noshifthgvs":
					case "-no_shift_hgvs":
						hgvsShift = false;
						break;

					case "-nolof":
						lossOfFunction = false; // Do not add LOF tag
						break;

					case "-oicr":
						useOicr = true; // Use OICR tag
						break;

					case "-sequenceontology":
						useSequenceOntology = true; // Use SO temrs
						break;

					//---
					// Input options
					//---
					case "-fi":
					case "-filterinterval":
						if ((i + 1) < args.length) filterIntervalFiles.add(args[++i]);
						else usage("Option '-fi' without config filter_interval_file argument");
						break;

					case "-i":
						// Input format
						usage("Command line option '-i' no longer supported: Only VCF format can be used for input and output.");
						break;

					//---
					// Filters
					//---
					case "-no-downstream":
						variantEffectResutFilter.add(EffectType.DOWNSTREAM);
						break;

					case "-no-upstream":
						variantEffectResutFilter.add(EffectType.UPSTREAM);
						break;

					case "-no-intergenic":
						variantEffectResutFilter.add(EffectType.INTERGENIC);
						break;

					case "-no-intron":
						variantEffectResutFilter.add(EffectType.INTRON);
						break;

					case "-no-utr":
						variantEffectResutFilter.add(EffectType.UTR_3_PRIME);
						variantEffectResutFilter.add(EffectType.UTR_3_DELETED);
						variantEffectResutFilter.add(EffectType.UTR_5_PRIME);
						variantEffectResutFilter.add(EffectType.UTR_5_DELETED);
						break;

					case "-no":
						String filterOut = "";
						if ((i + 1) < args.length) filterOut = args[++i];

						String filterOutArray[] = filterOut.split(",");
						for (String filterStr : filterOutArray) {
							if (filterStr.equalsIgnoreCase("utr")) {
								variantEffectResutFilter.add(EffectType.UTR_3_PRIME);
								variantEffectResutFilter.add(EffectType.UTR_3_DELETED);
								variantEffectResutFilter.add(EffectType.UTR_5_PRIME);
								variantEffectResutFilter.add(EffectType.UTR_5_DELETED);
							} else if (filterStr.equalsIgnoreCase("None")) ; // OK, nothing to do
							else variantEffectResutFilter.add(EffectType.valueOf(filterStr.toUpperCase()));
						}
						break;

					default:
						usage("Unknown option '" + arg + "'");
					}
				}
			} else if (genomeVer.isEmpty()) genomeVer = arg;
			else if (inputFile.isEmpty()) inputFile = arg;
			else usage("Unknown parameter '" + arg + "'");
		}

		//---
		// Sanity checks
		//---

		// Check: Do we have all required parameters?
		if (genomeVer == null || genomeVer.isEmpty()) usage("Missing genomer_version parameter");

		// Check input file
		if (inputFile.isEmpty()) inputFile = "-"; // Use STDIN as default
		else if (!Gpr.canRead(inputFile)) usage("Cannot read input file '" + inputFile + "'");

		// Read input files from file list?
		if (isFileList) {
			inputFiles = new ArrayList<>();
			for (String file : Gpr.readFileFatalError(inputFile).split("\n"))
				inputFiles.add(file);
		}
	}

	@Override
	public HashMap<String, String> reportValues() {
		HashMap<String, String> report = super.reportValues();
		annotateVcf.addReportValues(report);
		return report;
	}

	/**
	 * Run according to command line options
	 */
	@Override
	public boolean run() {
		return run(false) != null;
	}

	/**
	 * Run according to command line options
	 */
	public List<VcfEntry> run(boolean createList) {
		loadConfig(); // Read config file
		loadDb(); // Load database
		copyValuesToAnnotateVcf();
		return annotateVcf.run(createList);
	}

	public void setFormatVersion(EffFormatVersion formatVersion) {
		annotateVcf.setFormatVersion(formatVersion);
	}

	/**
	 * Show 'usage;' message and exit with an error code '-1'
	 * @param message
	 */
	@Override
	public void usage(String message) {
		if (message != null) {
			System.err.println("Error        :\t" + message);
			System.err.println("Command line :\t" + commandLineStr(false) + "\n");
		}

		System.err.println("snpEff version " + VERSION);
		System.err.println("Usage: snpEff [eff] [options] genome_version [input_file]");
		System.err.println("\n");
		System.err.println("\tvariants_file                   : Default is STDIN");
		System.err.println("\n");
		System.err.println("\nOptions:");
		System.err.println("\t-chr <string>                   : Prepend 'string' to chromosome name (e.g. 'chr1' instead of '1'). Only on TXT output.");
		System.err.println("\t-classic                        : Use old style annotations instead of Sequence Ontology and Hgvs.");
		System.err.println("\t-csvStats <file>                : Create CSV summary file.");
		System.err.println("\t-download                       : Download reference genome if not available. Default: " + download);
		System.err.println("\t-fileList                       : Input actually contains a list of files to process.");
		System.err.println("\t-s , -stats, -htmlStats         : Create HTML summary file.  Default is '" + DEFAULT_SUMMARY_HTML_FILE + "'");
		System.err.println("\t-noStats                        : Do not create stats (summary) file");
		System.err.println("\nResults filter options:");
		System.err.println("\t-fi , -filterInterval  <file>   : Only analyze changes that intersect with the intervals specified in this file (you may use this option many times)");
		System.err.println("\t-no-downstream                  : Do not show DOWNSTREAM changes");
		System.err.println("\t-no-intergenic                  : Do not show INTERGENIC changes");
		System.err.println("\t-no-intron                      : Do not show INTRON changes");
		System.err.println("\t-no-upstream                    : Do not show UPSTREAM changes");
		System.err.println("\t-no-utr                         : Do not show 5_PRIME_UTR or 3_PRIME_UTR changes");
		System.err.println("\t-no <effectType>                : Do not show 'EffectType'. This option can be used several times.");
		System.err.println("\nAnnotations options:");
		System.err.println("\t-cancer                         : Perform 'cancer' comparisons (Somatic vs Germline). Default: " + cancer);
		System.err.println("\t-cancerSamples <file>           : Two column TXT file defining 'oringinal \\t derived' samples.");
		System.err.println("\t-formatEff                      : Use 'EFF' field compatible with older versions (instead of 'ANN').");
		System.err.println("\t-geneId                         : Use gene ID instead of gene name (VCF output). Default: " + useGeneId);
		System.err.println("\t-hgvs                           : Use HGVS annotations for amino acid sub-field. Default: " + hgvs);
		System.err.println("\t-hgvsOld                        : Use old HGVS notation. Default: " + hgvsOld);
		System.err.println("\t-hgvs1LetterAa                  : Use one letter Amino acid codes in HGVS notation. Default: " + hgvsOneLetterAa);
		System.err.println("\t-hgvsTrId                       : Use transcript ID in HGVS notation. Default: " + hgvsTrId);
		System.err.println("\t-lof                            : Add loss of function (LOF) and Nonsense mediated decay (NMD) tags.");
		System.err.println("\t-noHgvs                         : Do not add HGVS annotations.");
		System.err.println("\t-noLof                          : Do not add LOF and NMD annotations.");
		System.err.println("\t-noShiftHgvs                    : Do not shift variants according to HGVS notation (most 3prime end).");
		System.err.println("\t-oicr                           : Add OICR tag in VCF file. Default: " + useOicr);
		System.err.println("\t-sequenceOntology               : Use Sequence Ontology terms. Default: " + useSequenceOntology);

		usageGenericAndDb();

		System.exit(-1);
	}

}
