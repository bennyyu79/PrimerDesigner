package nhs.genetics.cardiff;

import java.io.File;
import java.util.logging.Logger;

/**
 * Created by ml on 26/05/15.
 */
public class Configuration {

    private static final Logger log = Logger.getLogger(Configuration.class.getName());

    //parameters
    private static final int maxTargetLength = 450; //maximum sequence length to attempt a primer design before splitting
    private static final int maxPrimerDistance = 5000; //maxmimum distance between two blastn alignements to consider a viable amplicon
    private static final int padding = 350; //extra reference sequence surrounding target
    private static final int maxIndelLength = 10; //maximum length of an indel to be excluded
    private static final int spliceSitePadding = 20; //extra reference sequence surrounding the exon

    //exe paths
    private static final File primer3FilePath = new File("/share/apps/primer3-distros/release-2.3.6/primer3_core");
    private static final File blastnFilePath = new File("/usr/local/ncbi/blast/bin/blastn");
    private static final File bedtoolsFilePath = new File("/share/apps/bedtools-distros/bedtools2/bin/bedtools");

    //references
    private static final File exonsBed = new File("/Users/ml/Documents/Projects/PrimerDesigner/Homo_sapiens.GRCh37.75_cds.bed");
    private static final File blastnRefPath = new File("/data/db/human/gatk/2.8/b37/human_g1k_v37.fasta");
    private static final File referenceGenomeFasta = new File("/data/db/human/gatk/2.8/b37/human_g1k_v37.fasta");
    private static final File primerDatabaseFile = new File("amplicons.bed");
    private static final File excludedVariants = new File("/data/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf");
    private static final File primerMisprimingLibrary = new File("/Users/ml/Documents/Projects/PrimerDesigner/humrep_and_simple.txt");
    private static final File primer3Settings = new File("/Users/ml/Documents/Projects/PrimerDesigner/primer3_settings.txt");
    private static final String primerThermodynamicPararmetersPath = "/share/apps/primer3-distros/release-2.3.6/primer3_config/";

    private static final boolean debug = false;

    public static int getMaxPrimerDistance() {
        return maxPrimerDistance;
    }
    public static int getPadding() {
        return padding;
    }
    public static File getPrimer3FilePath() {
        return primer3FilePath;
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
    public static int getMaxIndelLength(){
        return  maxIndelLength;
    }
    public static boolean isDebug() {
        return debug;
    }
    public static int getSpliceSitePadding() {
        return spliceSitePadding;
    }
    public static File getReferenceGenomeFasta() {
        return referenceGenomeFasta;
    }
    public static File getExcludedVariants() {
        return excludedVariants;
    }
    public static File getExonsBed() {
        return exonsBed;
    }
    public static File getPrimerMisprimingLibrary() {
        return primerMisprimingLibrary;
    }
    public static File getPrimer3Settings() {
        return primer3Settings;
    }
    public static String getPrimerThermodynamicPararmetersPath() {
        return primerThermodynamicPararmetersPath;
    }
}
