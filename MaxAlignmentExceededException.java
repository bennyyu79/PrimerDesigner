package nhs.genetics.cardiff;

/**
 * Exception for primers with too many genome alignments
 *
 * @author  Matt Lyon
 * @version 1.0
 * @since   2015-06-15
 */

class MaxAlignmentExceededException extends Exception
{
    public MaxAlignmentExceededException(String message)
    {
        super(message);
    }
}
