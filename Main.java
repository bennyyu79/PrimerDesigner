package nhs.genetics.cardiff;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;

import java.io.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Main {

    //test

    private static final double version = 0.1;
    private static final Logger log = Logger.getLogger(Main.class.getName());

    public static void main(String[] args) {

        if (args.length != 4) {
            System.err.println("Usage: <ReferenceSequence> <ExcludedVariants> <Exons> <RegionOfInterest>");
            System.exit(1);
        }

        log.log(Level.INFO, "Primer designer v" + version);
        if (Configuration.isDebug()) log.log(Level.INFO, "Debugging mode");

        boolean hasAmpliconInDB;
        StringBuilder output = new StringBuilder();
        StringBuilder bedOutput = new StringBuilder();

        //write headers
        output.append("SuppliedTargetChr\tSuppliedTargetStart\tSuppliedTargetEnd\tSuppliedName\tDesignStart\tDesignEnd\tLeftPrimer\tRightPrimer\tLeftTm\tRightTm\tDesignPenalty\tSize\tAmpliconStart\tAmpliconEnd\n");

        //read regions of interest BED
        try (AbstractFeatureReader reader = AbstractFeatureReader.getFeatureReader(args[3], new BEDCodec(BEDCodec.StartOffset.ZERO), false)){
            Iterable<BEDFeature> iter = reader.iterator();

            //loop over regions of interest
            for (BEDFeature feature : iter) {

                log.log(Level.INFO, "Designing primer pair to cover supplied region of interest " + feature.getName() + " " + feature.getContig() + ":" + feature.getStart() + "-" + feature.getEnd());

                GenomicLocation roi = new GenomicLocation(feature.getContig(), feature.getStart(), feature.getEnd(), feature.getName());
                hasAmpliconInDB = false;

                //lookup roi in primer database; print amplicon(s) if exist
                for (GenomicLocation inHouseAmplicon : PrimerDatabase.isROICoveredByAmplicon(Configuration.getBedtoolsFilePath(), Configuration.getPrimerDatabaseFile(), roi)){

                    log.log(Level.INFO, "Target is already covered by amplicon(s) in primer database");

                    //write to temp
                    output.append(roi.getChromosome());
                    output.append("\t");
                    output.append(roi.getStartPosition());
                    output.append("\t");
                    output.append(roi.getEndPosition());
                    output.append("\t\t\t\t\t\t");
                    output.append(inHouseAmplicon.getEndPosition() - inHouseAmplicon.getStartPosition());
                    output.append("\t");
                    output.append(inHouseAmplicon.getChromosome());
                    output.append("\t");
                    output.append(inHouseAmplicon.getStartPosition());
                    output.append("\t");
                    output.append(inHouseAmplicon.getEndPosition());
                    output.append("\t");
                    output.append(inHouseAmplicon.getName());
                    output.append("\n");

                    hasAmpliconInDB = true;
                }

                //skip design if amplicon already exists
                if (hasAmpliconInDB) continue;

                //find overlapping regions of interest with exons
                HashSet<GenomicLocation> exonicRegions = new HashSet<>();
                for (String line : BedtoolsWrapper.getOverlappingFeatures(Configuration.getBedtoolsFilePath(), new File(args[2]), roi)){

                    String[] fields = line.split("\t");
                    exonicRegions.add(new GenomicLocation(fields[3], Integer.parseInt(fields[4]) - Configuration.getSpliceSitePadding(), Integer.parseInt(fields[5]) + Configuration.getSpliceSitePadding())); //splice site padding added

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

                        int numberOfWindows = ((finalROI.getEndPosition() - finalROI.getStartPosition()) + 1 / Configuration.getMaxTargetLength() + 1);
                        System.out.println(numberOfWindows);System.exit(0);


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

                    //convert to 1-based
                    finalROI.convertTo1Based();

                    log.log(Level.INFO, "Designing amplicon for target " + finalROI.getChromosome() + ":" + finalROI.getStartPosition() + "-" + finalROI.getEndPosition());

                    //get sequence
                    ReferenceSequence sequence = new ReferenceSequence(finalROI, new File(args[0]), new File(args[0] + ".fai"), Configuration.getPadding());
                    sequence.populateReferenceSequence();

                    if (sequence.isRefAllNSites()) {
                        log.log(Level.WARNING, "Could not design primer for target containing all N-sites: " + finalROI.getChromosome() + ":" + finalROI.getStartPosition() + "-" + finalROI.getEndPosition());
                        break;
                    }

                    //design primers
                    Primer3 primer3 = new Primer3(
                            sequence.getReferenceSequence(),
                            finalROI,
                            Configuration.getPadding(),
                            Configuration.getMaxPrimerDistance(),
                            Configuration.getMaxIndelExclusionLength(),
                            Configuration.getPrimer3RootFilePath(),
                            new File(args[1])
                    );
                    primer3.setExcludedRegions();
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

                        //no design available
                        if (primer3.getFilteredPrimerPairs().size() == 0){
                            output.append(roi.getChromosome() + "\t");
                            output.append(roi.getStartPosition() + "\t");
                            output.append(roi.getEndPosition() + "\t");
                            output.append("\n");
                        }

                        for (PrimerPair primerPair : primer3.getFilteredPrimerPairs()){

                            //print primers to table

                            //print target info
                            output.append(roi.getChromosome() + "\t");
                            output.append(roi.getStartPosition() + "\t");
                            output.append(roi.getEndPosition() + "\t");
                            output.append(roi.getName() + "\t");

                            //print design info
                            output.append(finalROI.getStartPosition() + "\t");
                            output.append(finalROI.getEndPosition() + "\t");

                            //print primers
                            output.append(primerPair.getLeftSequence() + "\t");
                            output.append(primerPair.getRightSequence() + "\t");
                            output.append(primerPair.getLeftTm() + "\t");
                            output.append(primerPair.getRightTm() + "\t");
                            output.append(primerPair.getPairPenalty() + "\t");
                            output.append(primerPair.getProductSize() + "\t");

                            //print amplicon coordinates
                            output.append(primerPair.getAmplifiableRegion().getStartPosition() + "\t");
                            output.append(primerPair.getAmplifiableRegion().getEndPosition() + "\t");

                            output.append("\n");

                            //print primers to bed
                            bedOutput.append(roi.getChromosome() + "\t");
                            bedOutput.append((primerPair.getAmplifiableRegion().getStartPosition() - primerPair.getLeftSequence().length() - 1) + "\t");
                            bedOutput.append((primerPair.getAmplifiableRegion().getEndPosition() + primerPair.getRightSequence().length()) + "\t");
                            bedOutput.append(roi.getName() + "\t");
                            bedOutput.append(Math.round(primerPair.getPairPenalty()) + "\t");
                            if (primerPair.getAmplifiableRegion().getStrand() == 1) bedOutput.append("+"); else bedOutput.append("-");
                            bedOutput.append(primerPair.getAmplifiableRegion().getStartPosition() - 1 + "\t");
                            bedOutput.append(primerPair.getAmplifiableRegion().getEndPosition() + "\t");
                            bedOutput.append("\n");

                        }

                    }

                }

                if (!Configuration.isDebug()){
                    ReferenceSequence paddedSequence = new ReferenceSequence(roi, new File (args[0]), new File (args[0] + ".fai"), Configuration.getPadding());
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

            //print primer pairs
            try (PrintWriter p = new PrintWriter("Primers.txt")){
                p.print(output.toString());
                p.close();
            } catch (IOException e){
                log.log(Level.SEVERE, e.toString());
            }

            //print primer pairs
            try (PrintWriter p = new PrintWriter("Primers.bed")){
                p.print(bedOutput.toString());
                p.close();
            } catch (IOException e){
                log.log(Level.SEVERE, e.toString());
            }

        }

    }

}
