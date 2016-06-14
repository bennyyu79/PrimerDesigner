package nhs.genetics.cardiff;

import java.util.ArrayList;
import java.util.MissingFormatArgumentException;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
/**
 * Class for holding PCR primer pairs
 *
 * @author  Matt Lyon
 * @version 1.0
 * @since   2015-04-17
 */
public class PrimerPair {

    private static final Logger log = Logger.getLogger(PrimerPair.class.getName());

    private ArrayList<String> primer3Output = new ArrayList<>();

    private GenomicLocation amplifiableRegion;

    private double pairPenalty, leftPenalty, rightPenalty, leftTm, rightTm, leftGC, rightGC, leftSelfAnyTh, rightSelfAnyTh, leftSelfEnd,
            rightSelfEnd, leftHairpin, rightHairpin, leftEndStability, rightEndStability, complAny, complEnd;
    private String leftSequence, rightSequence, leftPosition, rightPosition;
    private int productSize;

    public PrimerPair (ArrayList<String> primer3Output){
        this.primer3Output = primer3Output;
    }

    public void populatePrimerMetrics(){

        log.log(Level.FINE, "Parsing Primer3 output ...");

        //loop over primer3 output
        for (String line : primer3Output){

            String[] fields = line.split("=");

            //extract fields
            if (Pattern.matches("^PRIMER_PAIR_[0-9]+_PENALTY=.*", line)){
                this.pairPenalty = Double.parseDouble(fields[1]);
            } else if (Pattern.matches("^PRIMER_LEFT_[0-9]+_PENALTY=.*", line)){
                this.leftPenalty = Double.parseDouble(fields[1]);
            } else if (Pattern.matches("^PRIMER_RIGHT_[0-9]+_PENALTY=.*", line)){
                this.rightPenalty = Double.parseDouble(fields[1]);
            } else if (Pattern.matches("^PRIMER_LEFT_[0-9]+_SEQUENCE=.*", line)){
                this.leftSequence = fields[1];
            } else if (Pattern.matches("^PRIMER_RIGHT_[0-9]+_SEQUENCE=.*", line)){
                this.rightSequence = fields[1];
            } else if (Pattern.matches("^PRIMER_LEFT_[0-9]+=.*", line)){
                this.leftPosition = fields[1];
            } else if (Pattern.matches("^PRIMER_RIGHT_[0-9]+=.*", line)){
                this.rightPosition = fields[1];
            } else if (Pattern.matches("^PRIMER_LEFT_[0-9]+_TM=.*", line)){
                this.leftTm = Double.parseDouble(fields[1]);
            } else if (Pattern.matches("^PRIMER_RIGHT_[0-9]+_TM=.*", line)){
                this.rightTm = Double.parseDouble(fields[1]);
            } else if (Pattern.matches("^PRIMER_LEFT_[0-9]+_GC_PERCENT=.*", line)){
                this.leftGC = Double.parseDouble(fields[1]);
            } else if (Pattern.matches("^PRIMER_RIGHT_[0-9]+_GC_PERCENT=.*", line)){
                this.rightGC = Double.parseDouble(fields[1]);
            } else if (Pattern.matches("^PRIMER_LEFT_[0-9]+_SELF_ANY_TH=.*", line)){
                this.leftSelfAnyTh = Double.parseDouble(fields[1]);
            } else if (Pattern.matches("^PRIMER_RIGHT_[0-9]+_SELF_ANY_TH=.*", line)){
                this.rightSelfAnyTh = Double.parseDouble(fields[1]);
            } else if (Pattern.matches("^PRIMER_LEFT_[0-9]+_SELF_END_TH=.*", line)){
                this.leftSelfEnd = Double.parseDouble(fields[1]);
            } else if (Pattern.matches("^PRIMER_RIGHT_[0-9]+_SELF_END_TH=.*", line)){
                this.rightSelfEnd = Double.parseDouble(fields[1]);
            } else if (Pattern.matches("^PRIMER_LEFT_[0-9]+_HAIRPIN_TH=.*", line)){
                this.leftHairpin = Double.parseDouble(fields[1]);
            } else if (Pattern.matches("^PRIMER_RIGHT_[0-9]+_HAIRPIN_TH=.*", line)){
                this.rightHairpin = Double.parseDouble(fields[1]);
            } else if (Pattern.matches("^PRIMER_LEFT_[0-9]+_END_STABILITY=.*", line)){
                this.leftEndStability = Double.parseDouble(fields[1]);
            } else if (Pattern.matches("^PRIMER_RIGHT_[0-9]+_END_STABILITY=.*", line)){
                this.rightEndStability = Double.parseDouble(fields[1]);
            } else if (Pattern.matches("^PRIMER_PAIR_[0-9]+_COMPL_ANY_TH=.*", line)){
                this.complAny = Double.parseDouble(fields[1]);
            } else if (Pattern.matches("^PRIMER_PAIR_[0-9]+_COMPL_END_TH=.*", line)){
                this.complEnd = Double.parseDouble(fields[1]);
            } else if (Pattern.matches("^PRIMER_PAIR_[0-9]+_PRODUCT_SIZE=.*", line)){
                this.productSize = Integer.parseInt(fields[1]);
            }

        }

        //check if the minimum required fields are present
        if (leftSequence == null || rightSequence == null || leftPosition == null || rightPosition == null){
            throw new MissingFormatArgumentException("Primer pair was malformed");
        }

    }

    public double getPairPenalty() {
        return pairPenalty;
    }
    public double getLeftPenalty() {
        return leftPenalty;
    }
    public double getRightPenalty() {
        return rightPenalty;
    }
    public double getLeftTm() {
        return leftTm;
    }
    public double getRightTm() {
        return rightTm;
    }
    public double getLeftGC() {
        return leftGC;
    }
    public double getRightGC() {
        return rightGC;
    }
    public double getLeftSelfAnyTh() {
        return leftSelfAnyTh;
    }
    public double getRightSelfAnyTh() {
        return rightSelfAnyTh;
    }
    public double getLeftSelfEnd() {
        return leftSelfEnd;
    }
    public double getRightSelfEnd() {
        return rightSelfEnd;
    }
    public double getLeftHairpin() {
        return leftHairpin;
    }
    public double getRightHairpin() {
        return rightHairpin;
    }
    public double getLeftEndStability() {
        return leftEndStability;
    }
    public double getRightEndStability() {
        return rightEndStability;
    }
    public double getComplAny() {
        return complAny;
    }
    public double getComplEnd() {
        return complEnd;
    }
    public String getLeftSequence() {
        return leftSequence;
    }
    public String getRightSequence() {
        return rightSequence;
    }
    public String getLeftPosition() {
        return leftPosition;
    }
    public String getRightPosition() {
        return rightPosition;
    }
    public int getProductSize() {
        return productSize;
    }
    public GenomicLocation getAmplifiableRegion() {
        return amplifiableRegion;
    }

    public void setAmplifiableRegion(GenomicLocation amplifiableRegion) {
        this.amplifiableRegion = amplifiableRegion;
    }
}
