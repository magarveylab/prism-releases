package ca.mcmaster.magarveylab.prism.motif;

/**
 * Package for output from the FIMO program.
 * 
 * @author skinnider
 *
 */
public class FimoSearchResult {

	private String match;
	private int start;
	private int end;
	private double score;
	private double pValue;

	/**
	 * Instantiate a new FIMO results package.
	 * 
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
	 */
	public FimoSearchResult(String match, int start, int end, double score,
			double pValue) {
		this.match = match;
		this.start = start;
		this.end = end;
		this.score = score;
		this.pValue = pValue;
	}

	/**
	 * Get the sequence which was matched by the motif.
	 * 
	 * @return the sequence which was matched by the motif
	 */
	public String getMatch() {
		return match;
	}

	/**
	 * Get the start position of the motif occurrence.
	 * 
	 * @return the start position of the motif occurrence
	 */
	public int getStart() {
		return start;
	}

	/**
	 * Get the end position of the motif occurrence.
	 * 
	 * @return the end position of the motif occurrence
	 */
	public int getEnd() {
		return end;
	}

	/**
	 * Get the score for the motif occurrence.
	 * 
	 * @return the score for the motif occurrence
	 */
	public double getScore() {
		return score;
	}

	/**
	 * Get the p-value of the motif occurrence.
	 * 
	 * @return the p-value of the motif occurrence
	 */
	public double getPValue() {
		return pValue;
	}

}
