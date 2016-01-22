package nhs.genetics.cardiff;

import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by ml on 21/04/15.
 */
public class ReferenceSequence {

    private static final Logger log = Logger.getLogger(ReferenceSequence.class.getName());

    private String referenceSequence;
    private GenomicLocation location;
    private File fastaFilePath, indexFilePath;
    private int padding = 0;

    public ReferenceSequence(GenomicLocation location, File fastaFilePath, File indexFilePath){
        this.location = location;
        this.fastaFilePath = fastaFilePath;
        this.indexFilePath = indexFilePath;
    }
    public ReferenceSequence(GenomicLocation location, File fastaFilePath, File indexFilePath, int padding){
        this.location = location;
        this.fastaFilePath = fastaFilePath;
        this.indexFilePath = indexFilePath;
        this.padding = padding;
    }

    public void populateReferenceSequence(){ //1-based

        //read fasta index
        FastaSequenceIndex refGenomeIndex = new FastaSequenceIndex(indexFilePath);

        //get fasta sequence
        try(IndexedFastaSequenceFile refGenomeFasta = new IndexedFastaSequenceFile(fastaFilePath, refGenomeIndex)) {

            //get sequence
            this.referenceSequence = new String(refGenomeFasta.getSubsequenceAt(location.getContig(),
                    location.getStartPosition() - padding, location.getEndPosition() + padding).getBases(), "UTF-8");

            refGenomeFasta.close();

        } catch (UnsupportedEncodingException e){
            log.log(Level.SEVERE, "Problem converting nucleotide sequence: " + e.toString());
        } catch(IOException e){
            log.log(Level.SEVERE, "Problem reading reference genome: " + e.toString());
        }

    }

    public boolean isRefAllNSites(){

        for (char base : referenceSequence.toCharArray()){
            if (base != 'N'){
                return false;
            }
        }

        return true;
    }

    public String getReferenceSequence() {
        return referenceSequence;
    }
}
