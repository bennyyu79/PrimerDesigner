package nhs.genetics.cardiff;

import com.google.gson.Gson;

import java.io.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Main {

    private static final double version = 0.6;
    private static final Logger log = Logger.getLogger(Main.class.getName());

    public static void main(String[] args) {

        if (args.length != 5) {
            System.err.println("Usage: <Chrom> <Start> <Stop> <ConfigFilePath> <OutputType>");
            System.err.println("Coordinates should be 1-based");
            System.err.println("OutputType is JSON or BED");
            System.exit(1);
        }

        log.log(Level.INFO, "Primer designer v" + version);

        Configuration configuration = new Configuration(new File(args[3]));
        try {
            configuration.parseConfigurationFile();
        } catch (IOException e){
            log.log(Level.SEVERE, "Could not read config file: " + e.getMessage());
            System.exit(-1);
        }

        if (configuration.isDebug()) {
            log.log(Level.INFO, "Debugging mode");
        }

        StringBuilder bedOutput = new StringBuilder();
        Output output = new Output();
        Gson gson = new Gson();
        GenomicLocation suppliedROI = new GenomicLocation(args[0], Integer.parseInt(args[1]), Integer.parseInt(args[2]));
        suppliedROI.convertTo0Based();
        ArrayList<GenomicLocation> overlappingExonicRegionsOfInterest = new ArrayList<>();
        HashSet<GenomicLocation> mergedOverlappingExonicRegionsOfInterest = new HashSet<>();
        ArrayList<GenomicLocation> splitFinalRegionsOfInterest = new ArrayList<>();

        log.log(Level.INFO, "Designing primer pair to cover supplied region of interest " + suppliedROI.getContig() + ":" + suppliedROI.getStartPosition() + "-" + suppliedROI.getEndPosition());

        //find overlapping exons with ROI
        for (String feature : BedtoolsWrapper.getOverlappingFeatures(configuration.getBedtoolsFilePath(), configuration.getExonsBed(), suppliedROI)){

            String[] fields = feature.split("\t");
            overlappingExonicRegionsOfInterest.add(new GenomicLocation(fields[3], Integer.parseInt(fields[4]), Integer.parseInt(fields[5])));

            log.log(Level.INFO, "Target overlaps with exon " + fields[3] + ":" + fields[4] + "-" + fields[5]);
        }

        //loop over exonic overlaps and merge
        if (overlappingExonicRegionsOfInterest.size() > 0){

            //merge exonic overlaps
            for (GenomicLocation mergedOverlappingExonicROI : BedtoolsWrapper.mergeOverlappingFeatures(configuration.getBedtoolsFilePath(), overlappingExonicRegionsOfInterest)) {

                mergedOverlappingExonicRegionsOfInterest.add(mergedOverlappingExonicROI);

                log.log(Level.INFO, "Exonic target(s) were merged into " + mergedOverlappingExonicROI.getContig() + ":" + mergedOverlappingExonicROI.getStartPosition() + "-" + mergedOverlappingExonicROI.getEndPosition());
            }

        } else {
            log.log(Level.INFO, "Target does not overlap with any supplied exons");
            mergedOverlappingExonicRegionsOfInterest.add(suppliedROI); //could not find overlapping exons
        }

        //loop over final ROIs and split into amplifible amplicons
        //todo
        for (GenomicLocation finalROI : mergedOverlappingExonicRegionsOfInterest){

            /*split target
            if ((finalROI.getEndPosition() - finalROI.getStartPosition()) + 1 > Configuration.getMaxTargetLength()){

                int numberOfWindows = (((finalROI.getEndPosition() - finalROI.getStartPosition()) + 1) / Configuration.getMaxTargetLength()) + 1;

                log.log(Level.WARNING, "Target " + finalROI.getChromosome() + ":" + finalROI.getStartPosition() + "-" + finalROI.getEndPosition() + " exceeds max target length. Splitting into " + numberOfWindows + " fragments");

                for (String line : BedtoolsWrapper.splitRegionIntoWindows(Configuration.getBedtoolsFilePath(), finalROI, numberOfWindows)){
                    String[] fields = line.split("\t");
                    splitFinalRegionsOfInterest.add(new GenomicLocation(fields[0], Integer.parseInt(fields[1]), Integer.parseInt(fields[2])));
                }

                splitFinalRegionsOfInterest.add(suppliedROI);

            } else {
                splitFinalRegionsOfInterest.add(finalROI);
            }*/

            splitFinalRegionsOfInterest.add(finalROI);
        }

        //exonic and split ROIs
        for (GenomicLocation finalROI : splitFinalRegionsOfInterest) {

            //convert to 1-based
            finalROI.convertTo1Based();

            log.log(Level.INFO, "Designing amplicon for target " + finalROI.getContig() + ":" + finalROI.getStartPosition() + "-" + finalROI.getEndPosition());

            //get sequence
            ReferenceSequence sequence = new ReferenceSequence(finalROI, configuration.getReferenceGenomeFasta(), new File(configuration.getReferenceGenomeFasta() + ".fai"), configuration.getPadding());
            sequence.populateReferenceSequence();

            if (sequence.isRefAllNSites()) {
                log.log(Level.WARNING, "Could not design primer for target containing all N-sites: " + finalROI.getContig() + ":" + finalROI.getStartPosition() + "-" + finalROI.getEndPosition());
                break;
            }

            //design primers
            Primer3 primer3 = new Primer3(
                    sequence,
                    finalROI,
                    configuration
            );
            primer3.setExcludedRegions(configuration.getExcludedVariants(), configuration.getMaxIndelLength());
            primer3.callPrimer3();

            if (configuration.isDebug()){
                try (PrintWriter p = new PrintWriter(finalROI.getContig() + "_" + finalROI.getStartPosition() + "_" + finalROI.getEndPosition() + "_primer3out.txt")) {
                    for (String line : primer3.getPrimer3Output()) {
                        p.println(line);
                    }
                    p.close();
                } catch (IOException e) {
                    log.log(Level.SEVERE, e.getMessage());
                }
            }

            if (!configuration.isDebug()){

                primer3.splitPrimer3Output();
                primer3.checkPrimerAlignments();

                suppliedROI.convertTo1Based();

                for (PrimerPair primerPair : primer3.getFilteredPrimerPairs()){

                    //print primers to JSON
                    output.setChromosome(primerPair.getAmplifiableRegion().getContig());
                    output.setStartPosition(primerPair.getAmplifiableRegion().getStartPosition());
                    output.setEndPosition(primerPair.getAmplifiableRegion().getEndPosition());
                    output.setLeftSequence(primerPair.getLeftSequence());
                    output.setRightSequence(primerPair.getRightSequence());
                    output.setLeftTm(primerPair.getLeftTm());
                    output.setRightTm(primerPair.getRightTm());

                    //print primers to bed
                    bedOutput.append(suppliedROI.getContig());
                    bedOutput.append("\t");
                    bedOutput.append((primerPair.getAmplifiableRegion().getStartPosition() - primerPair.getLeftSequence().length()) - 1);
                    bedOutput.append("\t");
                    bedOutput.append((primerPair.getAmplifiableRegion().getEndPosition() + primerPair.getRightSequence().length()));
                    bedOutput.append("\t");
                    bedOutput.append(Math.round(primerPair.getPairPenalty()));
                    bedOutput.append("\t");
                    if (primerPair.getAmplifiableRegion().getStrand() == 1) bedOutput.append("+\t"); else bedOutput.append("-\t");
                    bedOutput.append((primerPair.getAmplifiableRegion().getStartPosition() - 1));
                    bedOutput.append("\t");
                    bedOutput.append(primerPair.getAmplifiableRegion().getEndPosition());
                    bedOutput.append("\n");

                }

            }

        }

        //write primers to stout
        if (!configuration.isDebug()){

            if (args[4].toUpperCase().equals("JSON")){
                System.out.print(gson.toJson(output));
            } else if (args[4].toUpperCase().equals("BED")){
                System.out.println(bedOutput.toString());
            }

        }

    }

}
