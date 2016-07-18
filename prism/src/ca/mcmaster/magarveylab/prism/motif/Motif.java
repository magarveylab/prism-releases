package ca.mcmaster.magarveylab.prism.motif;

import ca.mcmaster.magarveylab.enums.RibosomalPrecursorMotifs;
import ca.mcmaster.magarveylab.prism.util.exception.MotifException;

/**
 * A motif identified within a sequence using FIMO.
 * 
 * @author prees, skinnider 
 *
 */
public class Motif {

	private RibosomalPrecursorMotifs type;
	private String query;
	private String match;
	private int start;
	private int end;
	private double score;
	private double pValue;

	/**
	 * Instantiate a new motif from a FIMO search result package.
	 * 
	 * @param type
	 *            the motif type
	 * @param query
	 *            the open reading frame which was searched for the motif
	 * @param result
	 *            the FIMO search result
	 * @throws MotifException 
	 */
	public Motif(RibosomalPrecursorMotifs type, String query,
			FimoSearchResult result) throws MotifException {
		if (result.getStart() < 0)
			throw new MotifException(
					"Error: start cannot be less than 0! Start = " + start);
		if (result.getEnd() > query.length())
			throw new MotifException(
					"Error: end cannot be greater than length! End = " + end
							+ ", length = " + query.length());

		this.query = query;
		this.type = type;
		this.match = result.getMatch();
		this.start = result.getStart();
		this.end = result.getEnd();
		this.score = result.getScore();
		this.pValue = result.getPValue();
	}

	/**
	 * Instantiate a new motif.
	 * 
	 * @param type
	 *            the motif type
	 * @param query
	 *            the open reading frame which was searched for the motif
	 * @param match
	 *            the sequence which was matched by the motif
	 * @param start
	 *            the start position of the motif occurrence
	 * @param end
	 *            the end position of the motif occurrence
	 * @param score
	 *            the score for the motif occurrence
	 * @param pValue
	 *            the p-value of the motif occurrence
	 * @throws MotifException
	 */
	public Motif(RibosomalPrecursorMotifs type, String query, String match, int start,
			int end, double score, double pValue) throws MotifException {
		this.query = query;
		this.type = type;
		this.match = match;
		this.start = start;
		this.end = end;
		this.score = score;
		this.pValue = pValue;
	}

	/**
	 * Get the type of this motif.
	 * 
	 * @return the motif type
	 */
	public RibosomalPrecursorMotifs getType() {
		return type;
	}

	/**
	 * Get the sequence of the open reading frame which was searched for the
	 * motif.
	 * 
	 * @return the sequence of the open reading frame which was searched for the
	 *         motif
	 */
	public String getQuery() {
		return query;
	}

	/**
	 * Get the target sequence matched to this motif.
	 * 
	 * @return the target sequence matched to this motif 
	 */
	public String getMatch() {
		return match;
	}
	
	/**
	 * Get the start of the motif occurrence within the query sequence.
	 * 
	 * @return the start position of the motif occurrence
	 */
	public int getStart() {
		return start;
	}

	/**
	 * Get the end of the motif occurrence within the query sequence.
	 * 
	 * @return the end position of the motif occurrence
	 */
	public int getEnd() {
		return end;
	}

	/**
	 * Get the score assigned to the motif within the query sequence.
	 * 
	 * @return the score for the motif occurrence
	 */
	public double getScore() {
		return score;
	}

	/**
	 * Get the p-value assigned to the motif within the query sequence.
	 * 
	 * @return the p-value for the motif occurrence
	 */
	public double getPValue() {
		return pValue;
	}

	/**
	 * Get the substring of the query sequence corresponding to the motif.
	 * 
	 * @return the occurrence of the motif within the query sequence
	 */
	public String getMotifHit() {
		String motif = query.substring(start - 1, end);
		return motif;
	}

	@Override
	public String toString() {
		return getType() + " " + getPValue() + " " + getStart() + " " + getEnd() + " " + getMotifHit();
	}

	/**
	 * Check whether another motif overlaps with this one.
	 * 
	 * @param d2
	 *            the second motif
	 * @return true if their start-end ranges overlap
	 */
	public boolean overlaps(final Motif m2) {
		final int x1 = start;
		final int x2 = end;
		final int y1 = m2.getStart();
		final int y2 = m2.getEnd();
		return (x1 <= y2 && y1 <= x2);
	}

}