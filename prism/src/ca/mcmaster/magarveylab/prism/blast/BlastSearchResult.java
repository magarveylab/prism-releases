package ca.mcmaster.magarveylab.prism.blast;

import java.io.Serializable;

/**
 * Data package for a single row in a BLAST search result.
 * 
 * @author skinnider
 *
 */
public class BlastSearchResult implements Serializable {

	private static final long serialVersionUID = 7913519614233976975L;
	private double score;
	private double eValue;
	private double identity;
	private double positives;
	private int length;
	private int coverage;
	private String query;
	private String subject;

	/**
	 * Construct a new BLAST search result package, representing the
	 * intersection of a query and subject sequence.
	 * 
	 * @param query
	 *            the name of the query sequence
	 * @param subject
	 *            the name of the subject (database) sequence
	 * @param length
	 *            the length of the subject open reading frame
	 * @param score
	 *            the bitscore of the BLAST search result
	 * @param eValue
	 *            the E-value of the BLAST search result
	 * @param identity
	 *            the identity, in percent, of the BLAST search result
	 * @param positives
	 *            the positives, in percent, of the BLAST search result
	 * @param coverage
	 *            the coverage, in bases, of this BLAST search result
	 */
	public BlastSearchResult(String query, String subject, int length,
			double score, double eValue, double identity, double positives,
			int coverage) {
		this.query = query;
		this.subject = subject;
		this.length = length;
		this.score = score;
		this.eValue = eValue;
		this.identity = identity;
		this.positives = positives;
		this.coverage = coverage;
	}

	/**
	 * Get the name of the query sequence for this BLAST search result.
	 * 
	 * @return query sequence name
	 */
	public String query() {
		return query;
	}

	/**
	 * Get the name of the subject (database) sequence for this BLAST search
	 * result.
	 * 
	 * @return subject sequence name
	 */
	public String subject() {
		return subject;
	}

	/**
	 * Get the length of the subject (database) sequence for this BLAST search
	 * result.
	 * 
	 * @return subject sequence length
	 */
	public int length() {
		return length;
	}

	/**
	 * Get the bitscore of this BLAST search result.
	 * 
	 * @return bitscore
	 */
	public double score() {
		return score;
	}

	/**
	 * Set the bitscore of this BLAST search result.
	 * 
	 * @param score
	 *            bitscore
	 */
	public void setScore(double score) {
		this.score = score;
	}

	/**
	 * Get the E-value of this BLAST search result
	 * 
	 * @return e-value
	 */
	public double eValue() {
		return eValue;
	}

	/**
	 * Get the percent identity of this BLAST search result.
	 * 
	 * @return total identical bases divided by coverage
	 */
	public double identity() {
		return identity;
	}

	/**
	 * Get the percent positives of this BLAST search result.
	 * 
	 * @return total positive bases divided by coverage
	 */
	public double positives() {
		return positives;
	}

	/**
	 * Get the coverage of this BLAST search result.
	 * 
	 * @return the number of bases covered by the alignment
	 */
	public int coverage() {
		return coverage;
	}

	/**
	 * Set the subject of this BLAST search result.
	 * 
	 * @param subject
	 *            new subject to set
	 */
	public void setSubject(String subject) {
		this.subject = subject;
	}

}
