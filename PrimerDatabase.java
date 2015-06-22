package nhs.genetics.cardiff;

import java.io.File;
import java.util.HashSet;
import java.util.logging.Logger;

/**
 * Created by ml on 04/06/15.
 */
public class PrimerDatabase { //TODO: link to actual primer database: placeholder

    private static final Logger log = Logger.getLogger(PrimerDatabase.class.getName());

    public static HashSet<GenomicLocation> isROICoveredByAmplicon(File bedtoolsFilePath, File primerDatabaseFile, GenomicLocation target){

        HashSet<GenomicLocation> amplicons = new HashSet<>();

        for (String line : BedtoolsWrapper.getOverlappingFeatures(bedtoolsFilePath, primerDatabaseFile, target)) {
            String[] fields = line.split("\t");

            //check if target is covered 100%
            if (Integer.parseInt(fields[13]) == target.getEndPosition() - target.getStartPosition()) {
                amplicons.add(new GenomicLocation(fields[3], Integer.parseInt(fields[4]), Integer.parseInt(fields[5]), fields[6]));
            }

        }

        return amplicons;
    }

}
