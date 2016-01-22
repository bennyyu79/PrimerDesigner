package nhs.genetics.cardiff;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by ml on 20/04/15.
 */
public class MutationSurveyorReference {

    private static final Logger log = Logger.getLogger(MutationSurveyorReference.class.getName());
    private GenomicLocation roi;
    private ReferenceSequence paddedSequence;
    private int paddding;

    public MutationSurveyorReference(GenomicLocation roi, ReferenceSequence paddedSequence, int padding){
        this.roi = roi;
        this.paddedSequence = paddedSequence;
        this.paddding = padding;
    }

    public void writeMutationSurveyorSeqFile(){

        try(PrintWriter printWriter = new PrintWriter(roi.getContig() + "_" + roi.getStartPosition() + "_" + roi.getEndPosition() + ".seq")){

            printWriter.println("/Gene = \"GRCh37:" + roi.getContig() + ":" + roi.getStartPosition()  + "-" + roi.getEndPosition() + "\";");
            printWriter.println("/Exon_And_Note = \"\";");
            printWriter.println("/Reading Frame (1,2,3) = 1;");
            printWriter.println("/transl_table = ;");
            printWriter.println("/Remainder_Bases_of_the_Last_Exon = \"\";");
            printWriter.println("/Remainder_Bases_of_the_Next_Exon = \"\";");
            printWriter.println("/Number_of_the_First_Base = " + (roi.getStartPosition() - paddding) + ";");
            printWriter.println("/CDS = 0..0;");
            printWriter.println("/mCDSIndex = 1;");
            printWriter.println("/isLastmRNA = 0;");
            printWriter.println("/Exon_Base_Number = ;");
            printWriter.println("/mRNAIndex = ;");
            printWriter.println("/mRNARegion = ..;");
            printWriter.println("/Region of Interest = " + roi.getStartPosition() + ".." + roi.getEndPosition() + ";");
            printWriter.println("/Amplicon Id = \"\";");
            printWriter.println("/Amino Acid Sequence = ..;");
            printWriter.println("/Starting_vector_sequence = \"\";");
            printWriter.println("/Ending_vector_sequence = \"\";");
            printWriter.println("/NewVariation = \"\";");
            printWriter.println("/Translation = \"\";");
            printWriter.println((roi.getStartPosition() - paddding) + " " + paddedSequence.getReferenceSequence());

            printWriter.close();
        } catch (IOException e){
            log.log(Level.SEVERE, e.toString());
        }
    }

}
