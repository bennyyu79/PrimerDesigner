package nhs.genetics.cardiff;

import com.google.gson.Gson;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;

import java.io.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;

//TODO Make splitting framgents configurable

public class Main {

    private static final double version = 0.2;
    private static final Logger log = Logger.getLogger(Main.class.getName());

    public static void main(String[] args) {

        if (args.length != 1) {
            System.err.println("Usage: <RegionOfInterest.bed>");
            System.exit(1);
        }

        log.log(Level.INFO, "Primer designer v" + version);
        if (Configuration.isDebug()) log.log(Level.INFO, "Debugging mode");

        int n;
        StringBuilder bedOutput = new StringBuilder();
        Output output = new Output();
        Gson gson = new Gson();

        //read regions of interest BED
        try (AbstractFeatureReader reader = AbstractFeatureReader.getFeatureReader(args[0], new BEDCodec(BEDCodec.StartOffset.ZERO), false)){
            Iterable<BEDFeature> iter = reader.iterator();

            //loop over regions of interest
            for (BEDFeature feature : iter) {

                log.log(Level.INFO, "Designing primer pair to cover supplied region of interest " + feature.getName() + " " + feature.getContig() + ":" + feature.getStart() + "-" + feature.getEnd());

                GenomicLocation roi = new GenomicLocation(feature.getContig(), feature.getStart(), feature.getEnd(), feature.getName());

                //find overlapping regions of interest with exons
                HashSet<GenomicLocation> exonicRegions = new HashSet<>();
                for (String line : BedtoolsWrapper.getOverlappingFeatures(Configuration.getBedtoolsFilePath(), Configuration.getExonsBed(), roi)){

                    String[] fields = line.split("\t");
                    exonicRegions.add(new GenomicLocation(fields[3], Integer.parseInt(fields[4]), Integer.parseInt(fields[5])));

                    log.log(Level.INFO, "Target overlaps with exon " + fields[3] + ":" + fields[4] + "-" + fields[5]);
                }

                //loop over exonic overlaps and merge
                HashSet<GenomicLocation> finalRegionOfInterest = new HashSet<>();
                if (exonicRegions.size() > 0){

                    //merge exonic overlaps
                    for (GenomicLocation mergedTarget : BedtoolsWrapper.mergeOverlappingFeatures(Configuration.getBedtoolsFilePath(), exonicRegions)) {
                        finalRegionOfInterest.add(mergedTarget);

                        log.log(Level.INFO, "Exonic target(s) were merged into " + mergedTarget.getChromosome() + ":" + mergedTarget.getStartPosition() + "-" + mergedTarget.getEndPosition());
                    }

                } else {
                    log.log(Level.INFO, "Target does not overlap with any supplied exons");
                    finalRegionOfInterest.add(roi); //could not find overlapping exons
                }

                //loop over final ROIs and split into amplifible amplicons
                ArrayList<GenomicLocation> splitFinalRegionsOfInterest = new ArrayList<>();
                for (GenomicLocation finalROI : finalRegionOfInterest){

                    //split target
                    if ((finalROI.getEndPosition() - finalROI.getStartPosition()) + 1 > Configuration.getMaxTargetLength()){

                        int numberOfWindows = (((finalROI.getEndPosition() - finalROI.getStartPosition()) + 1) / Configuration.getMaxTargetLength()) + 1;

                        log.log(Level.WARNING, "Target " + finalROI.getChromosome() + ":" + finalROI.getStartPosition() + "-" + finalROI.getEndPosition() + " exceeds max target length. Splitting into " + numberOfWindows + " fragments");

                        for (String line : BedtoolsWrapper.splitRegionIntoWindows(Configuration.getBedtoolsFilePath(), finalROI, numberOfWindows)){
                            String[] fields = line.split("\t");
                            splitFinalRegionsOfInterest.add(new GenomicLocation(fields[0], Integer.parseInt(fields[1]), Integer.parseInt(fields[2])));
                        }

                    } else {
                        splitFinalRegionsOfInterest.add(finalROI);
                    }

                }

                //exonic and split ROIs
                for (GenomicLocation finalROI : splitFinalRegionsOfInterest) {

                    //pair number
                    n = 0;

                    //convert to 1-based
                    finalROI.convertTo1Based();

                    log.log(Level.INFO, "Designing amplicon for target " + finalROI.getChromosome() + ":" + finalROI.getStartPosition() + "-" + finalROI.getEndPosition());

                    //get sequence
                    ReferenceSequence sequence = new ReferenceSequence(finalROI, Configuration.getReferenceGenomeFasta(), new File(Configuration.getReferenceGenomeFasta() + ".fai"), Configuration.getPadding());
                    sequence.populateReferenceSequence();

                    if (sequence.isRefAllNSites()) {
                        log.log(Level.WARNING, "Could not design primer for target containing all N-sites: " + finalROI.getChromosome() + ":" + finalROI.getStartPosition() + "-" + finalROI.getEndPosition());
                        break;
                    }

                    //design primers
                    Primer3 primer3 = new Primer3(
                            sequence,
                            finalROI,
                            Configuration.getPadding(),
                            Configuration.getMaxPrimerDistance(),
                            Configuration.getPrimer3FilePath(),
                            Configuration.getPrimerMisprimingLibrary(),
                            Configuration.getPrimer3Settings(),
                            Configuration.getPrimerThermodynamicPararmetersPath()
                    );
                    primer3.setExcludedRegions(Configuration.getExcludedVariants(), Configuration.getMaxIndelLength());
                    primer3.callPrimer3();

                    if (Configuration.isDebug()){
                        try (PrintWriter p = new PrintWriter(finalROI.getChromosome() + "_" + finalROI.getStartPosition() + "_" + finalROI.getEndPosition() + "_primer3out.txt")) {
                            for (String line : primer3.getPrimer3Output()) {
                                p.println(line);
                            }
                            p.close();
                        } catch (IOException e) {
                            log.log(Level.SEVERE, e.getMessage());
                        }
                    }

                    if (!Configuration.isDebug()){

                        primer3.splitPrimer3Output();
                        primer3.checkPrimerAlignments();

                        roi.convertTo1Based();

                        for (PrimerPair primerPair : primer3.getFilteredPrimerPairs()){
                            ++n;

                            //print primers to JSON
                            output.setChromosome(primerPair.getAmplifiableRegion().getChromosome());
                            output.setStartPosition(primerPair.getAmplifiableRegion().getStartPosition());
                            output.setEndPosition(primerPair.getAmplifiableRegion().getEndPosition());
                            output.setLeftSequence(primerPair.getLeftSequence());
                            output.setRightSequence(primerPair.getRightSequence());
                            output.setLeftTm(primerPair.getLeftTm());
                            output.setRightTm(primerPair.getRightTm());
                            output.setRightTm(primerPair.getPairPenalty());

                            //print primers to bed
                            bedOutput.append(roi.getChromosome() + "\t");
                            bedOutput.append((primerPair.getAmplifiableRegion().getStartPosition() - primerPair.getLeftSequence().length()) - 1 + "\t");
                            bedOutput.append((primerPair.getAmplifiableRegion().getEndPosition() + primerPair.getRightSequence().length()) + "\t");
                            bedOutput.append("\t");
                            bedOutput.append(Math.round(primerPair.getPairPenalty()) + "\t");
                            if (primerPair.getAmplifiableRegion().getStrand() == 1) bedOutput.append("+\t"); else bedOutput.append("-\t");
                            bedOutput.append((primerPair.getAmplifiableRegion().getStartPosition() - 1) + "\t");
                            bedOutput.append(primerPair.getAmplifiableRegion().getEndPosition());
                            bedOutput.append("\n");

                        }

                    }

                }

                if (!Configuration.isDebug()){
                    ReferenceSequence paddedSequence = new ReferenceSequence(roi, Configuration.getReferenceGenomeFasta(), new File(Configuration.getReferenceGenomeFasta() + ".fai"), Configuration.getPadding());
                    paddedSequence.populateReferenceSequence();

                    MutationSurveyorReference mutationSurveyorReference = new MutationSurveyorReference(roi, paddedSequence, Configuration.getPadding());
                    mutationSurveyorReference.writeMutationSurveyorSeqFile();
                }

            }// done reading supplied ROIs

            reader.close();
        } catch (IOException e){
            log.log(Level.SEVERE, e.toString());
        }

        if (!Configuration.isDebug()){

            //write output to stdout
            System.out.print(gson.toJson(output));

            //print primer pairs to BED
            try (PrintWriter p = new PrintWriter("Primers.bed")){
                p.print(bedOutput.toString());
                p.close();
            } catch (IOException e){
                log.log(Level.SEVERE, e.toString());
            }

        }

    }

}
