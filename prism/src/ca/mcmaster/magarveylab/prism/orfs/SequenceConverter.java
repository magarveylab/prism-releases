package ca.mcmaster.magarveylab.prism.orfs;

import org.biojava3.core.exceptions.CompoundNotFoundError;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.RNASequence;

/**
 * Utility class for biological sequence conversion and manipulation
 * @author skinnider
 *
 */
public class SequenceConverter {
	
	/**
	 * Converts DNA sequence to AA sequence
	 * @param dna	The DNA sequence to be converted
	 */
	public static String convertDNAToAA(String dna) throws CompoundNotFoundError {
		DNASequence dnaSeq = new DNASequence(dna);
		RNASequence rnaSeq = dnaSeq.getRNASequence();
		ProteinSequence aaSeq = rnaSeq.getProteinSequence();
		return aaSeq.toString();
	}
	
	/**
	 * Converts a DNA sequence to its complement
	 * @param dna	The DNA sequence to be converted
	 * @return		The complementary DNA sequence as a string
	 */
	public static String getComplementaryDNASequence(String dna) throws CompoundNotFoundError {
		DNASequence dnaSeq = new DNASequence(dna);
		return dnaSeq.getComplement().getSequenceAsString();
	}
	
	/**
	 * Converts a DNA sequence to its complement and reverses the resulting string
	 * @param dna	The DNA sequence to be converted
	 * @return		The complementary DNA sequence as a string
	 */
	public static String getReverseComplementaryDNASequence(String dna) throws CompoundNotFoundError {
		DNASequence dnaSeq = new DNASequence(dna);
		return dnaSeq.getReverseComplement().getSequenceAsString();
	}
}
