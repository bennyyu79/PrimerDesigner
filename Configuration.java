package nhs.genetics.cardiff;

import java.io.File;

/**
 * Created by ml on 26/05/15.
 */
public class Configuration {

    //parameters
    private static final int maxTargetLength = 400; //maximum sequence length to attempt a primer design before splitting
    private static final int maxPrimerDistance = 5000; //maxmimum distance between two blastn alignements to consider a viable amplicon
    private static final int padding = 350; //extra reference sequence surrounding target
    private static final int maxIndelExclusionLength = 10; //maximum length of an indel to be excluded
    private static final int spliceSitePadding = 20; //extra reference sequence surrounding the exon
    private static final File primer3RootFilePath = new File("/Users/ml/Downloads/Primer3-2.3.6");
    private static final File blastnFilePath = new File("/usr/local/ncbi/blast/bin/blastn");
    private static final File blastnRefPath = new File("/Users/ml/Downloads/blastn_db/b37_blast");
    private static final File bedtoolsFilePath = new File("/share/apps/bedtools-distros/bedtools2/bin/bedtools");
    private static final File primerDatabaseFile = new File("/Users/ml/Downloads/Primer3-2.3.6/amplicons.bed");
    private static final boolean debug = false;


    public static int getMaxPrimerDistance() {
        return maxPrimerDistance;
    }
    public static int getPadding() {
        return padding;
    }
    public static File getPrimer3RootFilePath() {
        return primer3RootFilePath;
    }
    public static File getBlastnFilePath() {
        return blastnFilePath;
    }
    public static File getBlastnRefPath() {
        return blastnRefPath;
    }
    public static File getBedtoolsFilePath() {
        return bedtoolsFilePath;
    }
    public static int getMaxTargetLength() {
        return maxTargetLength;
    }
    public static File getPrimerDatabaseFile() {
        return primerDatabaseFile;
    }
    public static int getMaxIndelExclusionLength(){
        return  maxIndelExclusionLength;
    }
    public static boolean isDebug() {
        return debug;
    }
    public static int getSpliceSitePadding() {
        return spliceSitePadding;
    }
}
