package nhs.genetics.cardiff;

import java.io.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * A class for running bedtools functions http://bedtools.readthedocs.io/en/latest
 *
 * @author  Matt Lyon
 * @version 1.0
 * @since   2015-06-04
 */
public class BedtoolsWrapper {
    private static final Logger log = Logger.getLogger(BedtoolsWrapper.class.getName());

    public static ArrayList<String> getOverlappingFeatures(File bedtoolsFilePath, File bedFilePath, GenomicLocation lookup){

        StringBuilder targetBedInput = new StringBuilder();
        ArrayList<String> bedtoolsOutput = new ArrayList<>();

        //convert target to bed record
        targetBedInput.append(lookup.getContig());
        targetBedInput.append("\t");
        targetBedInput.append(lookup.getStartPosition());
        targetBedInput.append("\t");
        targetBedInput.append(lookup.getEndPosition());
        targetBedInput.append("\n");

        //intersect target with exon bed
        try{
            ProcessBuilder exeBuilder = new ProcessBuilder(
                    bedtoolsFilePath.toString(),
                    "intersect",
                    "-a",
                    "-",
                    "-b",
                    bedFilePath.toString(),
                    "-wo"
            );

            Process process = exeBuilder.start();

            OutputStream stdin = process.getOutputStream();
            InputStream stdout = process.getInputStream();

            BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(stdin));

            writer.write(targetBedInput.toString());
            writer.flush();
            writer.close();

            Scanner scanner = new Scanner(stdout);
            while (scanner.hasNextLine()) {
                bedtoolsOutput.add(scanner.nextLine());
            }

            if (process.waitFor() != 0){
                throw new RuntimeException("Problem invoking bedtools intersect, exit code: " + process.exitValue());
            }

        } catch (IOException e){
            log.log(Level.SEVERE, e.toString());
        } catch (InterruptedException e){
            log.log(Level.SEVERE, e.toString());
        }

        return  bedtoolsOutput;
    }
    public static HashSet<GenomicLocation> mergeOverlappingFeatures(File bedtoolsFilePath, ArrayList<GenomicLocation> features){

        HashSet<GenomicLocation> mergedFeatures = new HashSet<>();
        ArrayList<String> bedtoolsOutput = new ArrayList<>();
        StringBuilder bedFeaturesToMerge = new StringBuilder();

        for (GenomicLocation loc : features){
            bedFeaturesToMerge.append(loc.getContig());
            bedFeaturesToMerge.append("\t");
            bedFeaturesToMerge.append(loc.getStartPosition());
            bedFeaturesToMerge.append("\t");
            bedFeaturesToMerge.append(loc.getEndPosition());
            bedFeaturesToMerge.append("\n");
        }

        //merge overlapping bed features
        try{
            ProcessBuilder exeBuilder = new ProcessBuilder(
                    bedtoolsFilePath.toString(),
                    "merge"
            );

            Process process = exeBuilder.start();

            OutputStream stdin = process.getOutputStream();
            InputStream stdout = process.getInputStream();

            BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(stdin));

            writer.write(bedFeaturesToMerge.toString());
            writer.flush();
            writer.close();

            Scanner scanner = new Scanner(stdout);
            while (scanner.hasNextLine()) {
                bedtoolsOutput.add(scanner.nextLine());
            }

            if (process.waitFor() != 0){
                throw new RuntimeException("Problem invoking bedtools merge, exit code: " + process.exitValue());
            }

        } catch (IOException e){
            log.log(Level.SEVERE, e.toString());
        } catch (InterruptedException e){
            log.log(Level.SEVERE, e.toString());
        }

        //convert output to GenomicLocations; make unique
        for (String i : bedtoolsOutput){
            String[] fields = i.split("\t");
            mergedFeatures.add(new GenomicLocation(fields[0], Integer.parseInt(fields[1]), Integer.parseInt(fields[2])));
        }

        return mergedFeatures;
    }

    public static ArrayList<String> splitRegionIntoWindows(File bedtoolsFilePath, GenomicLocation regionOfInterest, int numberOfWindows){

        String numberOfWindowsString = "" + numberOfWindows;
        StringBuilder targetBedInput = new StringBuilder();
        ArrayList<String> bedtoolsOutput = new ArrayList<>();

        //convert target to bed record
        targetBedInput.append(regionOfInterest.getContig());
        targetBedInput.append("\t");
        targetBedInput.append(regionOfInterest.getStartPosition());
        targetBedInput.append("\t");
        targetBedInput.append(regionOfInterest.getEndPosition());
        targetBedInput.append("\n");

        //intersect target with exon bed
        try{
            ProcessBuilder exeBuilder = new ProcessBuilder(
                    bedtoolsFilePath.toString(),
                    "makewindows",
                    "-b",
                    "-",
                    "-n",
                    numberOfWindowsString
            );

            Process process = exeBuilder.start();

            OutputStream stdin = process.getOutputStream();
            InputStream stdout = process.getInputStream();

            BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(stdin));

            writer.write(targetBedInput.toString());
            writer.flush();
            writer.close();

            Scanner scanner = new Scanner(stdout);
            while (scanner.hasNextLine()) {
                bedtoolsOutput.add(scanner.nextLine());
            }

            if (process.waitFor() != 0){
                throw new RuntimeException("Problem invoking bedtools makewindows, exit code: " + process.exitValue());
            }

        } catch (IOException e){
            log.log(Level.SEVERE, e.toString());
        } catch (InterruptedException e){
            log.log(Level.SEVERE, e.toString());
        }

        return  bedtoolsOutput;
    }

}
