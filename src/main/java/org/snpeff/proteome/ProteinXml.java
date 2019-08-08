package org.snpeff.proteome;

import org.snpeff.interval.*;
import org.snpeff.interval.codonChange.CodonChange;
import org.snpeff.snpEffect.EffectType;
import org.snpeff.snpEffect.VariantEffect;
import org.snpeff.util.Timer;
import org.snpeff.util.Tuple;
import org.snpeff.vcf.EffFormatVersion;
import org.snpeff.vcf.VcfEffect;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class ProteinXml {

    private static int xmlDepth = 0;

    public static void writeProteinXml(String xmlProtLocation, String organism, HashMap<Transcript, List<VariantEffect>> transcriptVariants, List<String> sampleNames) throws IOException, XMLStreamException {
        XMLStreamWriter writer = XMLOutputFactory.newInstance().createXMLStreamWriter(new FileWriter(xmlProtLocation));
        writer.writeStartDocument();
        writeStartElement(writer, "mzLibProteinDb");

        for (Transcript tr : transcriptVariants.keySet())
        {
            if (!tr.isProteinCoding() || tr.isErrorStartCodon() ||
                    tr.proteinTrimmed().isEmpty() || tr.proteinTrimmed().endsWith("?"))
                continue;

            writeStartElement(writer, "entry");
            writeStartElement(writer, "accession");
            writer.writeCharacters(tr.getId());
            writer.writeEndElement(); xmlDepth--;

            if (tr.getId() != null) {
                writeStartElement(writer, "name");
                writer.writeCharacters(tr.getId());
                writer.writeEndElement(); xmlDepth--;
            }

            Gene gene = (Gene)tr.findParent(Gene.class);
            String name = gene.getGeneName();
            String id = gene.getId();
            String primary = name != null ? name : id;
            writeStartElement(writer, "gene");
            writeStartElement(writer, "name");
            writer.writeAttribute("type", "primary");
            writer.writeCharacters(primary);
            writer.writeEndElement(); xmlDepth--;
            writeStartElement(writer, "name");
            writer.writeAttribute("type", "accession");
            writer.writeCharacters(id);
            writer.writeEndElement(); xmlDepth--;
            writer.writeEndElement(); xmlDepth--;

            if (organism != null) {
                writeStartElement(writer, "organism");
                writeStartElement(writer, "name");
                writer.writeAttribute("type", "scientific");
                writer.writeCharacters(organism);
                writer.writeEndElement(); xmlDepth--;
                writer.writeEndElement(); xmlDepth--;
            }

            // TODO: I can't get at VcfGenotype from Variant, so will need to do that in Spritz. (Tried including VcfEntry in Variant initialization, and that back reference didn't work...)
            // TODO: I can also do the depth stuff there.
            // That would also make this output more compatible with UniProt
            // TODO: Implement a better Transcript.protein() method for this output: trim to stop codon and possibly go beyond coding domain if there's a frameshift that extends
            for (VariantEffect var : transcriptVariants.get(tr))
            {
//                if (var.getEffectType() == EffectType.STOP_LOST) Timer.showStdErr("stop_lost variant in protein xml writer: "+ (new VcfEffect(var, EffFormatVersion.DEFAULT_FORMAT_VERSION, false, false)));
//                if (var.getEffectType() == EffectType.STOP_GAINED) Timer.showStdErr("stop_gained variant in protein xml writer: "+ (new VcfEffect(var, EffFormatVersion.DEFAULT_FORMAT_VERSION, false, false)));
                if (var.getFunctionalClass().compareTo(VariantEffect.FunctionalClass.MISSENSE) < 0 && var.getAaRef() == var.getAaAlt()) { // only annotate indels and nonsynonymous variations
                    continue;
                }
                if (var.getVariant() == null){
                    Timer.showStdErr("getVariant() returned null for: " + var.toString());
                    continue;
                }
                else{
                    // cases to avoid
                    if (var.getVariant().line ==  null){
                        Timer.showStdErr("line for getVariant() returned null for: " + (new VcfEffect(var, EffFormatVersion.DEFAULT_FORMAT_VERSION, false, false)).toString());
                        continue;
                    }
                    if (var.getCodonNum() < 0){
                        Timer.showStdErr("line for getCodonNum() was negative for: " + (new VcfEffect(var, EffFormatVersion.DEFAULT_FORMAT_VERSION, false, false)).toString());
                        continue;
                    }

                    // multiple stop codons or weird indel case, but keep stop gains
                    String protAlt = tr.apply(var.getVariant()).proteinTrimmed();
                    boolean longerThanProtein = var.getCodonNum() >= protAlt.length();
                    boolean stopGain = var.getCodonNum() == protAlt.length() && var.getEffectType() == EffectType.STOP_GAINED;
                    if (longerThanProtein && !stopGain)
                    {
                        Timer.showStdErr("line where getCodonNum() " + Integer.toString(var.getCodonNum()) + " was longer than the protein length " + Integer.toString(tr.apply(var.getVariant()).proteinTrimmed().length()) +  ": " + (new VcfEffect(var, EffFormatVersion.DEFAULT_FORMAT_VERSION, false, false)).toString());
                        continue;
                    }

                    // indel cases to avoid
                    boolean isIndel = var.getAaAlt().endsWith("?") || var.getAaAlt() == "";
                    if (isIndel && var.getCodonNum() >= tr.apply(var.getVariant()).proteinTrimmed().length())
                    {
                        Timer.showStdErr("line where getCodonNum() was longer than the variant protein length: " + (new VcfEffect(var, EffFormatVersion.DEFAULT_FORMAT_VERSION, false, false)).toString());
                        continue;
                    }
                    else if (isIndel && tr.apply(var.getVariant()).proteinTrimmed().substring(var.getCodonNum()).endsWith("?"))
                    {
                        Timer.showStdErr("line where indel variant transcript does not have a stop codon: " + (new VcfEffect(var, EffFormatVersion.DEFAULT_FORMAT_VERSION, false, false)).toString());
                        continue;
                    }
                }

                // write the xml
                writeStartElement(writer, "feature");
                writer.writeAttribute("type", "sequence variant");
                String[] descArr = var.getVariant().line.split("\t");
                descArr[7] = "ANN=" + (new VcfEffect(var, EffFormatVersion.DEFAULT_FORMAT_VERSION, false, false)).toString();// var.toString().replace("\t", "|");
                writer.writeAttribute("description", String.join("\\t", descArr));
                writeStartElement(writer, "original");
                String refProt = tr.proteinTrimmed();
                String altProt = tr.apply(var.getVariant()).proteinTrimmed();
                int refSeqStart = var.getCodonNum() < refProt.length() ? var.getCodonNum() : refProt.length() - 1; // for stop loss
                int altSeqStart = var.getCodonNum() < altProt.length() ? var.getCodonNum() : altProt.length() - 1; // for stop gain
                String refSeq = var.getAaRef();
                String varSeq = var.getAaAlt();
                int pos = var.getCodonNum();
                if (var.getAaAlt().endsWith("?") || var.getAaAlt() == "" || var.getCodonNum() >= refProt.length()){
                    refSeq = refProt.substring(refSeqStart);
                    varSeq = altProt.substring(refSeqStart);
                    pos = refSeqStart;
                }
                else if (var.getCodonNum() >= altProt.length()){
                    refSeq = refProt.substring(altSeqStart);
                    varSeq = altProt.substring(altSeqStart);
                    pos = altSeqStart;
                }
                writer.writeCharacters(refSeq);
                writer.writeEndElement(); xmlDepth--; // reference
                writeStartElement(writer, "variation");
                writer.writeCharacters(varSeq);
                writer.writeEndElement(); xmlDepth--; // alternate

                writeStartElement(writer, "location");
                if (var.getAaNetChange() != null)
                {
                    writeStartElement(writer, "position");
                    writer.writeAttribute("position", Integer.toString(pos + 1));
                    writer.writeEndElement(); xmlDepth--;
                }
                else
                {
                    writeStartElement(writer, "begin");
                    writer.writeAttribute("position", Integer.toString(pos + 1));
                    writer.writeEndElement(); xmlDepth--;
                    writeStartElement(writer, "end");
                    writer.writeAttribute("position", Integer.toString(pos + var.getAaNetChange().length()));
                    writer.writeEndElement(); xmlDepth--;
                }
                writer.writeEndElement(); xmlDepth--; // location
                xmlDepth--; writePretty(writer); writer.writeEndElement();  // feature
            }

            float codonsSoFar = 0;
            List<Cds> cdsList = tr.getCds();
            for (int i = 0; i < cdsList.size() - 1; i++) {
                Cds cds = cdsList.get(i);
                float codonCount = (float)(cds.getEnd() - cds.getStart() + 1) / (float)CodonChange.CODON_SIZE;
                float codonNum = codonCount + codonsSoFar;
                int startCodonNum = (int)Math.ceil(codonNum); // which codon is the start of this CDS in?
                int endCodonNum = startCodonNum == codonNum ? startCodonNum + 1 : startCodonNum; // completed the last codon?
                codonsSoFar = codonNum;
                if (endCodonNum >= tr.proteinTrimmed().length())
                    continue;

                writeStartElement(writer, "feature");
                writer.writeAttribute("type", "splice site");
                writer.writeAttribute("description", cds.getChromosomeName() + "\\t" + cds.getStrand() + "\\t" + cds.getStart()+ "\\t" + cds.getEnd() + "\\t" + cdsList.get(i+1).getStart()+ "\\t" + cdsList.get(i+1).getEnd());
                writeStartElement(writer, "location");
                if (endCodonNum - startCodonNum == 0)
                {
                    writeStartElement(writer, "position");
                    writer.writeAttribute("position", Integer.toString(startCodonNum));
                    writer.writeEndElement(); xmlDepth--;
                }
                else
                {
                    writeStartElement(writer, "begin");
                    writer.writeAttribute("position", Integer.toString(startCodonNum));
                    writer.writeEndElement(); xmlDepth--;
                    writeStartElement(writer, "end");
                    writer.writeAttribute("position", Integer.toString(endCodonNum));
                    writer.writeEndElement(); xmlDepth--;
                }
                writer.writeEndElement(); xmlDepth--; // location
                xmlDepth--; writePretty(writer); writer.writeEndElement();  // feature
            }

            writeStartElement(writer, "sequence");
            writer.writeAttribute("length", Integer.toString(tr.proteinTrimmed().length()));
            writer.writeCharacters(tr.proteinTrimmed());
            writer.writeEndElement(); xmlDepth--; // sequence
            xmlDepth--; writePretty(writer); writer.writeEndElement(); // entry
        }
        xmlDepth--; writePretty(writer); writer.writeEndElement(); // mzLibProteinDb
        writer.writeEndDocument();
        writer.flush();
    }

    private static void writePretty(XMLStreamWriter writer) throws XMLStreamException {
        writer.writeCharacters("\n");
        for (int x = 0; x < xmlDepth; x++){
            writer.writeCharacters("  ");
        }
    }

    private static void writeStartElement(XMLStreamWriter writer, String localName) throws XMLStreamException {
        writePretty(writer);
        writer.writeStartElement(localName);
        xmlDepth++;
    }
}
