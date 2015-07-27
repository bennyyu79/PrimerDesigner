package nhs.genetics.cardiff;

/**
 * Created by ml on 22/07/2015.
 */
public class Output {

    private String chromosome, leftSequence, rightSequence;
    private int startPosition, endPosition;
    private double leftTm,rightTm;

    public Output(){

    }

    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }
    public void setLeftSequence(String leftSequence) {
        this.leftSequence = leftSequence;
    }
    public void setRightSequence(String rightSequence) {
        this.rightSequence = rightSequence;
    }
    public void setStartPosition(int startPosition) {
        this.startPosition = startPosition;
    }
    public void setEndPosition(int endPosition) {
        this.endPosition = endPosition;
    }
    public void setLeftTm(double leftTm) {
        this.leftTm = leftTm;
    }
    public void setRightTm(double rightTm) {
        this.rightTm = rightTm;
    }

}
