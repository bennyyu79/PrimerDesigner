package nhs.genetics.cardiff;

import java.io.*;
import java.nio.charset.MalformedInputException;
import java.util.MissingFormatArgumentException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Program configuration class
 *
 * @author  Matt Lyon
 * @version 1.0
 * @since   2015-05-26
 */
public class Configuration {

    private final Logger log = Logger.getLogger(Configuration.class.getName());
    private File configurationPath;

    //parameters
    private int maxTargetLength = 450; //maximum sequence length to attempt a primer design before splitting
    private int maxPrimerDistance = 5000; //maxmimum distance between two blastn alignements to consider a viable amplicon
    private int padding = 350; //extra reference sequence surrounding target
    private int maxIndelLength = 10; //maximum length of an indel to be excluded
    private int maxExactMatches = 1;
    private double minSimilarity = 0.95;
    private boolean debug = false;

    private File exonsBed, blastnRefPath, referenceGenomeFasta, primerDatabaseFile, excludedVariants, primerMisprimingLibrary, primer3Settings, primer3FilePath, blastnFilePath, bedtoolsFilePath, primerThermodynamicPararmetersPath;

    public Configuration(File configurationPath) {
        this.configurationPath = configurationPath;
    }

    public void parseConfigurationFile() throws IOException {

        log.log(Level.INFO, "Parsing config file");

        String line;

        //read config file
        try (BufferedReader reader = new BufferedReader(new FileReader(configurationPath))){

            while ((line = reader.readLine()) != null) {

                if (!line.equals("")) {
                    String[] fields = line.split("=");

                    if (fields[0].equals("exonsBed")){
                        exonsBed = new File(fields[1]);
                    } else if (fields[0].equals("blastnRefPath")){
                        blastnRefPath = new File(fields[1]);
                    } else if (fields[0].equals("referenceGenomeFasta")){
                        referenceGenomeFasta = new File(fields[1]);
                    } else if (fields[0].equals("primerDatabaseFile")){
                        primerDatabaseFile = new File(fields[1]);
                    } else if (fields[0].equals("excludedVariants")){
                        excludedVariants = new File(fields[1]);
                    } else if (fields[0].equals("primerMisprimingLibrary")){
                        primerMisprimingLibrary = new File(fields[1]);
                    } else if (fields[0].equals("primer3Settings")){
                        primer3Settings = new File(fields[1]);
                    } else if (fields[0].equals("primer3FilePath")){
                        primer3FilePath = new File(fields[1]);
                    } else if (fields[0].equals("blastnFilePath")){
                        blastnFilePath = new File(fields[1]);
                    } else if (fields[0].equals("bedtoolsFilePath")){
                        bedtoolsFilePath = new File(fields[1]);
                    } else if (fields[0].equals("primerThermodynamicPararmetersPath")){
                        primerThermodynamicPararmetersPath = new File(fields[1]);
                    }

                }

            }

            reader.close();
        } catch (NullPointerException e){
            throw new IllegalArgumentException("Config file malformed");
        }

    }

    public int getMaxTargetLength() {
        return maxTargetLength;
    }
    public int getMaxPrimerDistance() {
        return maxPrimerDistance;
    }
    public int getPadding() {
        return padding;
    }
    public int getMaxIndelLength() {
        return maxIndelLength;
    }
    public int getMaxExactMatches() {
        return maxExactMatches;
    }
    public double getMinSimilarity() {
        return minSimilarity;
    }
    public boolean isDebug() {
        return debug;
    }
    public File getExonsBed() {
        return exonsBed;
    }
    public File getBlastnRefPath() {
        return blastnRefPath;
    }
    public File getReferenceGenomeFasta() {
        return referenceGenomeFasta;
    }
    public File getPrimerDatabaseFile() {
        return primerDatabaseFile;
    }
    public File getExcludedVariants() {
        return excludedVariants;
    }
    public File getPrimerMisprimingLibrary() {
        return primerMisprimingLibrary;
    }
    public File getPrimer3Settings() {
        return primer3Settings;
    }
    public File getPrimer3FilePath() {
        return primer3FilePath;
    }
    public File getBlastnFilePath() {
        return blastnFilePath;
    }
    public File getBedtoolsFilePath() {
        return bedtoolsFilePath;
    }
    public File getPrimerThermodynamicPararmetersPath() {
        return primerThermodynamicPararmetersPath;
    }

}