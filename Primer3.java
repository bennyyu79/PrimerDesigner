package nhs.genetics.cardiff;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;

/**
 * Primer3 wrapper. Call to design primers around a target
 *
 * @author  Matt Lyon
 * @version 1.0
 * @since   2015-04-17
 */
public class Primer3 {

    private static final Logger log = Logger.getLogger(Primer3.class.getName());

    private ReferenceSequence referenceSequence;
    private ArrayList<PrimerPair> candidatePrimerPairs = new ArrayList<>();
    private ArrayList<PrimerPair> filteredPrimerPairs = new ArrayList<>();
    private ArrayList<String> primer3Output = new ArrayList<>();
    private StringBuilder excludedRegions = new StringBuilder();
    private HashMap<String, ArrayList<GenomicLocation>> primerAlignments = new HashMap<>();
    private GenomicLocation targetLocation;
    private Configuration configuration;

    //TODO: ligate M13 adapters
    //TODO: Re-calculate primer hairpin with M13 adapter

    public Primer3(ReferenceSequence referenceSequence,
                   GenomicLocation targetLocation,
                   Configuration configuration) {

        this.referenceSequence = referenceSequence;
        this.targetLocation = targetLocation;
        this.configuration = configuration;
    }

    public void callPrimer3(){
        log.log(Level.INFO, "Calling Primer3 ...");

        StringBuilder primer3input = new StringBuilder();

        primer3input.append("SEQUENCE_TEMPLATE=");
        primer3input.append(referenceSequence.getReferenceSequence());
        primer3input.append("\n");
        primer3input.append("SEQUENCE_TARGET=");
        primer3input.append((configuration.getPadding() + 1) + "," + (targetLocation.getEndPosition() - targetLocation.getStartPosition() + 1));
        primer3input.append("\n");
        primer3input.append("SEQUENCE_EXCLUDED_REGION=");
        primer3input.append(excludedRegions.toString());
        primer3input.append("\n");
        primer3input.append("PRIMER_MISPRIMING_LIBRARY=");
        primer3input.append(configuration.getPrimerMisprimingLibrary());
        primer3input.append("\n");
        primer3input.append("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=");
        primer3input.append(configuration.getPrimerThermodynamicPararmetersPath());
        primer3input.append("/\n");
        primer3input.append("PRIMER_EXPLAIN_FLAG=1\n");
        primer3input.append("=");

        log.log(Level.FINE, "Passing Primer3 args");
        log.log(Level.FINE, primer3input.toString());
        try{
            ProcessBuilder exeBuilder;

            if (configuration.isDebug()){
                exeBuilder = new ProcessBuilder(
                        configuration.getPrimer3FilePath().toString(),
                        "-p3_settings_file=" + configuration.getPrimer3Settings().getAbsolutePath(),
                        "-echo_settings_file",
                        "-format_output"
                );
            } else {
                exeBuilder = new ProcessBuilder(
                        configuration.getPrimer3FilePath().toString(),
                        "-p3_settings_file=" + configuration.getPrimer3Settings().getAbsolutePath(),
                        "-echo_settings_file"
                );
            }

            Process process = exeBuilder.start();

            OutputStream stdin = process.getOutputStream();
            InputStream stdout = process.getInputStream();

            BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(stdin));

            writer.write(primer3input.toString());
            writer.flush();
            writer.close();

            Scanner scanner = new Scanner(stdout);
            while (scanner.hasNextLine()) {
                primer3Output.add(scanner.nextLine());
            }

            if (process.waitFor() != 0){
                throw new RuntimeException("Problem invoking Primer3, exit code: " + process.exitValue());
            }

        } catch (IOException e){
            log.log(Level.SEVERE, e.toString());
        } catch (InterruptedException e){
            log.log(Level.SEVERE, e.toString());
        }

        if (configuration.isDebug()) {
            try(PrintWriter p = new PrintWriter(targetLocation.getContig() + "_" + targetLocation.getStartPosition() + "_" + targetLocation.getEndPosition() + "_primer3in.txt")){
                p.write(primer3input.toString());
            }catch (IOException e){
                log.log(Level.SEVERE, e.getMessage());
            }
        }

    }
    public void splitPrimer3Output(){

        int primerNo = 0;
        ArrayList<String> primerPairOutput = new ArrayList<>();

        //Split primer3 output into individual designs and parse
        for (String line : primer3Output){

            //check record is primer
            if (Pattern.matches("^PRIMER_PAIR_[0-9]+[_|=].*", line) || Pattern.matches("^PRIMER_LEFT_[0-9]+[_|=].*", line) || Pattern.matches("^PRIMER_RIGHT_[0-9]+[_|=].*", line)) {

                //split primer fields
                String[] fields = line.split("_|=");

                //check if this is a new record
                if (Integer.parseInt(fields[2]) == primerNo) {

                    //current primer
                    primerPairOutput.add(line);

                } else {

                    //new primer

                    //dispose of last primer & skip null primers
                    try {
                        PrimerPair primerPair = new PrimerPair(primerPairOutput);
                        primerPair.populatePrimerMetrics();

                        candidatePrimerPairs.add(primerPair);
                    } catch (MissingFormatArgumentException e){
                        log.log(Level.FINE, e.getMessage());
                    }

                    primerPairOutput.clear();

                    //load new primer
                    primerPairOutput.add(line);

                    //reset primerNo value for new primer
                    primerNo = Integer.parseInt(fields[2]);
                }

            } else if (line.equals("=")){

                //dispose of last primer & skip null primers
                try {
                    PrimerPair primerPair = new PrimerPair(primerPairOutput);
                    primerPair.populatePrimerMetrics();

                    candidatePrimerPairs.add(primerPair);
                } catch (MissingFormatArgumentException e){
                    log.log(Level.FINE, e.getMessage());
                }

            } else {
                //primer summary
                log.log(Level.FINE, line);
            }

        }

    }

    public void checkPrimerAlignments(){

        boolean hasCorrectAlignment;
        int alignments;

        log.log(Level.INFO, "Testing " + candidatePrimerPairs.size() + " candidate primer pair(s).");

        //loop over candidate primer pairs
        for (int j = 0; j < candidatePrimerPairs.size(); ++j){

            log.log(Level.INFO, "Checking primer specificity for candidate pair: " + (j + 1));

            hasCorrectAlignment = false;
            alignments = 0;

            try {

                //blast primers if not already done
                if (!primerAlignments.containsKey(candidatePrimerPairs.get(j).getLeftSequence())){
                    primerAlignments.put(candidatePrimerPairs.get(j).getLeftSequence(), Blast.callShortQueryBlast(candidatePrimerPairs.get(j).getLeftSequence(), configuration.getBlastnFilePath(), configuration.getBlastnRefPath(), configuration.getMaxExactMatches(), configuration.getMinSimilarity()));
                }

                if (!primerAlignments.containsKey(candidatePrimerPairs.get(j).getRightSequence())){
                    primerAlignments.put(candidatePrimerPairs.get(j).getRightSequence(), Blast.callShortQueryBlast(candidatePrimerPairs.get(j).getRightSequence(), configuration.getBlastnFilePath(), configuration.getBlastnRefPath(), configuration.getMaxExactMatches(), configuration.getMinSimilarity()));
                }

            } catch (MaxAlignmentExceededException e){
                log.log(Level.INFO, "Skipping pair: " + e.getMessage());
                continue;
            }

            //loop over all primer alignments for this pair
            for (GenomicLocation leftAlignment : primerAlignments.get(candidatePrimerPairs.get(j).getLeftSequence())) {
                for (GenomicLocation rightAlignment : primerAlignments.get(candidatePrimerPairs.get(j).getRightSequence())) {

                    //skip alignments on different contigs
                    if (!leftAlignment.getContig().equals(rightAlignment.getContig())) {
                        continue;
                    }

                    //check primers are orientated correctly for amplification
                    if (
                            leftAlignment.getStartPosition() < leftAlignment.getEndPosition() && //check orientation
                            rightAlignment.getStartPosition() > rightAlignment.getEndPosition() &&//check orientation
                            rightAlignment.getStartPosition() - leftAlignment.getStartPosition() > 0 && //check primers point towards each other
                            rightAlignment.getStartPosition() - leftAlignment.getStartPosition() < configuration.getMaxPrimerDistance() //check amplicon is less than maxSize;

                    ) {

                        alignments++; //+ strand

                        //check primer alignment start and length
                        String[] leftPrimerOffsetAndLength = candidatePrimerPairs.get(j).getLeftPosition().split(",");
                        String[] rightPrimerOffsetAndLength = candidatePrimerPairs.get(j).getRightPosition().split(",");

                        if (
                                Integer.parseInt(leftPrimerOffsetAndLength[0]) + (targetLocation.getStartPosition() - configuration.getPadding()) == leftAlignment.getStartPosition() + 1 &&
                                Integer.parseInt(leftPrimerOffsetAndLength[1]) == ((leftAlignment.getEndPosition() - leftAlignment.getStartPosition()) + 1) &&
                                Integer.parseInt(rightPrimerOffsetAndLength[0]) + (targetLocation.getStartPosition() - configuration.getPadding()) == rightAlignment.getStartPosition() + 1 &&
                                Integer.parseInt(rightPrimerOffsetAndLength[1]) == ((rightAlignment.getStartPosition() - rightAlignment.getEndPosition()) + 1) &&
                                targetLocation.getContig().equals(leftAlignment.getContig())) {

                            hasCorrectAlignment = true;

                            GenomicLocation amplifibleRegion = new GenomicLocation(leftAlignment.getContig(), leftAlignment.getStartPosition() + Integer.parseInt(leftPrimerOffsetAndLength[1]), rightAlignment.getStartPosition() - Integer.parseInt(rightPrimerOffsetAndLength[1]));
                            amplifibleRegion.setStrand(1);

                            candidatePrimerPairs.get(j).setAmplifiableRegion(amplifibleRegion);
                        }

                    } else if (

                            rightAlignment.getEndPosition() < rightAlignment.getStartPosition() && //check orientation
                            leftAlignment.getEndPosition() > leftAlignment.getStartPosition() &&//check orientation
                            leftAlignment.getEndPosition() - rightAlignment.getEndPosition() > 0 && //check primers point towards each other
                            leftAlignment.getEndPosition() - rightAlignment.getEndPosition() < configuration.getMaxPrimerDistance()) {//check amplicon is less than maxSize;

                        alignments++; //- strand

                        //check primer alignment start and length
                        String[] leftPrimerOffsetAndLength = candidatePrimerPairs.get(j).getLeftPosition().split(",");
                        String[] rightPrimerOffsetAndLength = candidatePrimerPairs.get(j).getRightPosition().split(",");

                        if (
                                Integer.parseInt(leftPrimerOffsetAndLength[0]) + (targetLocation.getStartPosition() - configuration.getPadding()) == leftAlignment.getStartPosition() + 1 &&
                                Integer.parseInt(leftPrimerOffsetAndLength[1]) == ((leftAlignment.getEndPosition() - leftAlignment.getStartPosition()) + 1) &&
                                Integer.parseInt(rightPrimerOffsetAndLength[0]) + (targetLocation.getStartPosition() - configuration.getPadding()) == rightAlignment.getStartPosition() + 1 &&
                                Integer.parseInt(rightPrimerOffsetAndLength[1]) == ((rightAlignment.getStartPosition() - rightAlignment.getEndPosition()) + 1) &&
                                targetLocation.getContig().equals(leftAlignment.getContig())) {

                            hasCorrectAlignment = true;

                            GenomicLocation amplifibleRegion = new GenomicLocation(leftAlignment.getContig(), leftAlignment.getStartPosition() + Integer.parseInt(leftPrimerOffsetAndLength[1]), rightAlignment.getStartPosition() - Integer.parseInt(rightPrimerOffsetAndLength[1]));
                            amplifibleRegion.setStrand(-1);

                            candidatePrimerPairs.get(j).setAmplifiableRegion(amplifibleRegion);
                        }

                    }

                }
            } //done looping over primer alignments

            //bank specific primer pairs
            if (!hasCorrectAlignment){
                log.log(Level.WARNING, "Could not find correct alignment for: " + (j + 1));
            } else {
                if (alignments == 1){
                    filteredPrimerPairs.add(candidatePrimerPairs.get(j));
                    break; //only deliver one good primer pair
                } else {
                    log.log(Level.INFO, "Could not find specific alignment for: " + (j + 1));
                }
            }

        }

    }

    public void setExcludedRegions(File vcfFilePath, int maxIndelLength){

        //get nearby dbSNP entries
        VCFFileReader vcfFile = new VCFFileReader(vcfFilePath, new File(vcfFilePath + ".idx"));
        Iterator<VariantContext> it = vcfFile.query(targetLocation.getContig(), targetLocation.getStartPosition() - configuration.getPadding(), targetLocation.getEndPosition() + configuration.getPadding());
        ArrayList<Long> sortedExcludedPositions = new ArrayList<>();
        HashSet<Long> excludedPositions = new HashSet<>();

        while(it.hasNext()) {

            VariantContext poly = it.next();

            if (poly.isSNP()){

                //loop over known variant bases
                for (long n = poly.getStart() + 1; n < poly.getEnd() + 2; ++n){ //1-based

                    //exclude regions (convert from chrom to seq pos)
                    excludedPositions.add(n - (targetLocation.getStartPosition() - configuration.getPadding()));

                }

            } else if (poly.isIndel()){

                //loop over known variant bases
                for (long n = poly.getStart() + 2; n < poly.getEnd() + 2; ++n){ //1-based

                    //skip indels greater than 10bp
                    if (poly.getEnd() - poly.getStart() > maxIndelLength){
                        continue;
                    }

                    //exclude regions (convert from chrom to seq pos)
                    excludedPositions.add(n - (targetLocation.getStartPosition() - configuration.getPadding()));

                }

            }

        }

        //convert hashset to array
        for (long n : excludedPositions){
            sortedExcludedPositions.add(n);
        }

        //sort excluded bases
        Collections.sort(sortedExcludedPositions);

        for (long n : sortedExcludedPositions){
            excludedRegions.append(n);
            excludedRegions.append(",");
            excludedRegions.append(1); //length
            excludedRegions.append(" ");
        }

        vcfFile.close();
    }

    public ArrayList<PrimerPair> getFilteredPrimerPairs() {
        return filteredPrimerPairs;
    }
    public ArrayList<String> getPrimer3Output() {
        return primer3Output;
    }
}
