package org.snpeff.proteome;

import org.snpeff.interval.Gene;
import org.snpeff.interval.Transcript;
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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class ProteinXml {

    private static int xmlDepth = 0;

    public static void writeVariantProteinXml(String xmlProtLocation, String organism, ArrayList<Tuple<Transcript, VariantEffect>> altTranscriptVariants) throws IOException, XMLStreamException {
        XMLStreamWriter writer = XMLOutputFactory.newInstance().createXMLStreamWriter(new FileWriter(xmlProtLocation));
        writer.writeStartDocument();
        writeStartElement(writer, "mzLibProteinDb");

        for (Tuple<Transcript, VariantEffect> trve : altTranscriptVariants)
        {
            Transcript tr = trve.first;
            VariantEffect var = trve.second;
            if (!tr.isProteinCoding() || tr.protein().isEmpty()  || tr.isErrorStartCodon() || !tr.protein().contains("*")
                    || (var.getFunctionalClass().compareTo(VariantEffect.FunctionalClass.MISSENSE) < 0 && var.getAaRef() == var.getAaAlt())) // only annotate nonsynonymous variations
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

            writeStartElement(writer, "feature");
            writer.writeAttribute("type", "sequence variant");
            writer.writeAttribute("description", var.getVariant().line);
            writeStartElement(writer, "original");
            writer.writeCharacters(var.getAaRef());
            writer.writeEndElement(); xmlDepth--; // reference
            writeStartElement(writer, "variation");
            writer.writeCharacters(var.getAaAlt());
            writer.writeEndElement(); xmlDepth--; // alternate

            writeStartElement(writer, "location");
            if (var.getAaNetChange() != null)
            {
                writeStartElement(writer, "position");
                writer.writeAttribute("position", Integer.toString(var.getCodonNum()));
                writer.writeEndElement(); xmlDepth--;
            }
            else
            {
                writeStartElement(writer, "begin");
                writer.writeAttribute("position", Integer.toString(var.getCodonNum()));
                writer.writeEndElement(); xmlDepth--;
                writeStartElement(writer, "end");
                writer.writeAttribute("position", Integer.toString(var.getCodonNum() + var.getAaNetChange().length() - 1));
                writer.writeEndElement(); xmlDepth--;
            }
            writer.writeEndElement(); xmlDepth--; // location
            xmlDepth--; writePretty(writer); writer.writeEndElement();  // feature

            writeStartElement(writer, "sequence");
            writer.writeAttribute("length", Integer.toString(tr.protein().length()));
            writer.writeCharacters(tr.proteinTrimmed());
            writer.writeEndElement(); xmlDepth--; // sequence
            xmlDepth--; writePretty(writer); writer.writeEndElement(); // entry
        }
        xmlDepth--; writePretty(writer); writer.writeEndElement(); // mzLibProteinDb
        writer.writeEndDocument();
        writer.flush();
    }

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
                if (var.getFunctionalClass().compareTo(VariantEffect.FunctionalClass.MISSENSE) < 0 && var.getAaRef() == var.getAaAlt()) { // only annotate indels and nonsynonymous variations
                    continue;
                }
                if (var.getVariant() == null){
                    Timer.showStdErr("getVariant() returned null for: " + var.toString());
                    continue;
                }
                else if (var.getVariant().line ==  null){
                    Timer.showStdErr("line for getVariant() returned null for: " + (new VcfEffect(var, EffFormatVersion.DEFAULT_FORMAT_VERSION, false, false)).toString());
                    continue;
                }

                writeStartElement(writer, "feature");
                writer.writeAttribute("type", "sequence variant");
                String[] descArr = var.getVariant().line.split("\t");
                descArr[7] = "ANN=" + (new VcfEffect(var, EffFormatVersion.DEFAULT_FORMAT_VERSION, false, false)).toString();// var.toString().replace("\t", "|");
                writer.writeAttribute("description", String.join("\\t", descArr));
                writeStartElement(writer, "original");
                String refSeq = var.getAaAlt().endsWith("?") || var.getAaAlt() == "" ? tr.proteinTrimmed().substring(var.getCodonNum()) : var.getAaRef();
                writer.writeCharacters(refSeq);
                writer.writeEndElement(); xmlDepth--; // reference
                writeStartElement(writer, "variation");
                String varSeq = var.getAaAlt().endsWith("?") || var.getAaAlt() == "" ? tr.apply(var.getVariant()).proteinTrimmed().substring(var.getCodonNum()) : var.getAaAlt();
                writer.writeCharacters(varSeq);
                writer.writeEndElement(); xmlDepth--; // alternate

                writeStartElement(writer, "location");
                if (var.getAaNetChange() != null)
                {
                    writeStartElement(writer, "position");
                    writer.writeAttribute("position", Integer.toString(var.getCodonNum() + 1));
                    writer.writeEndElement(); xmlDepth--;
                }
                else
                {
                    writeStartElement(writer, "begin");
                    writer.writeAttribute("position", Integer.toString(var.getCodonNum() + 1));
                    writer.writeEndElement(); xmlDepth--;
                    writeStartElement(writer, "end");
                    writer.writeAttribute("position", Integer.toString(var.getCodonNum() + var.getAaNetChange().length()));
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
