package nhs.genetics.cardiff;

import java.io.*;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by ml on 27/04/15.
 */

public class Blast {

    private static final Logger log = Logger.getLogger(Blast.class.getName());

    public static ArrayList<GenomicLocation> callShortQueryBlast(String query, File blastnFilePath, File blastnRefPath){

        log.log(Level.FINE, "Calling short query blast ...");
        log.log(Level.FINE, "Sequence: " + query);

        ArrayList<String> output = new ArrayList<>();
        ArrayList<GenomicLocation> alignments = new ArrayList<>();

        try{

            ProcessBuilder builder = new ProcessBuilder(
                    blastnFilePath.toString(),
                    "-db", blastnRefPath.toString(),
                    "-task" , "blastn-short",
                    "-outfmt", "6"
            );

            Process process = builder.start();

            OutputStream stdin = process.getOutputStream();
            InputStream stdout = process.getInputStream();

            BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(stdin));

            writer.write(query);
            writer.flush();
            writer.close();

            Scanner scanner = new Scanner(stdout);
            while (scanner.hasNextLine()) {
                output.add(scanner.nextLine());
            }

            if (process.waitFor() != 0){
                throw new RuntimeException("Problem invoking Blast, exit code: " + process.exitValue());
            }

        } catch (IOException e){
            log.log(Level.SEVERE, e.toString());
        } catch (RuntimeException e){
            log.log(Level.SEVERE, e.toString());
        } catch (InterruptedException e){
            log.log(Level.SEVERE, e.toString());
        }

        //extract genome coordinates from blast output
        for (String s : output){
            String[] fields = s.split("\t");
            alignments.add(new GenomicLocation(fields[1], Integer.parseInt(fields[8]), Integer.parseInt(fields[9])));
        }

        return alignments;

    }

}
