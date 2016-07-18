package ca.mcmaster.magarveylab.prism.data;

import ca.mcmaster.magarveylab.prism.motif.Motif;
import ca.mcmaster.magarveylab.prism.motif.MotifList;

/**
 * A cleaved leader peptide. In ribosomally synthesized and post-translationally
 * modified natural product (RiPP) biosynthesis, a precursor peptide undergoes
 * cleavage to produce a propeptide, which is subsequently enzymatically
 * tailored.
 * 
 * @author skinnider
 *
 */
public class Propeptide {

	private String sequence;
	private int start;
	private int end;
	private MotifList motifs = new MotifList();

	/**
	 * Instantiate a new propeptide.
	 * 
	 * @param sequence
	 *            the parent open reading frame sequence
	 * @param start
	 *            the start of the propeptide
	 * @param end
	 *            the end of the propeptide
	 */
	public Propeptide(String sequence, int start, int end) {
		this.sequence = sequence;
		this.start = start;
		this.end = end;
	}

	/**
	 * Get the sequence of the propeptide.
	 * 
	 * @return the propeptide sequence
	 */
	public String getSequence() {
		return sequence;
	}

	/**
	 * Set the sequence of the propeptide.
	 * 
	 * @param sequence
	 *            the propeptide sequence
	 */
	public void setSequence(String sequence) {
		this.sequence = sequence;
	}

	/**
	 * Get the point on the parent open reading frame at which the propeptide
	 * starts.
	 * 
	 * @return the start of the propeptide on the parent orf
	 */
	public int getStart() {
		return start;
	}

	/**
	 * Set the point on the parent open reading frame at which the propeptide
	 * starts.
	 * 
	 * @param start
	 *            the start of the propeptide on the parent orf
	 */
	public void setStart(int start) {
		this.start = start;
	}

	/**
	 * Get the point on the parent open reading frame at which the propeptide
	 * ends.
	 * 
	 * @return the end of the propeptide on the parent orf
	 */
	public int getEnd() {
		return end;
	}

	/**
	 * Set the point on the parent open reading frame at which the propeptide
	 * ends.
	 * 
	 * @param start
	 *            the end of the propeptide on the parent orf
	 */
	public void setEnd(int end) {
		this.end = end;
	}

	/**
	 * Get the length of the propeptide.
	 * 
	 * @return the absolute value of the difference between the start and end
	 *         points
	 */
	public int getLength() {
		return Math.abs(end - start);
	}

	/**
	 * Get the motifs used to identify this propeptide.
	 * 
	 * @return the motifs used to identify this propeptide
	 */
	public MotifList getMotifs() {
		return motifs;
	}

	/**
	 * Set the motifs used to identify this propeptide.
	 * 
	 * @param motifs
	 *            the motifs used to identify this propeptide
	 */
	public void setMotifs(MotifList motifs) {
		this.motifs = motifs;
	}

	/**
	 * Associate a new motif with the identification of this propeptide.
	 * 
	 * @param motif
	 *            motif to associate with this propeptide
	 */
	public void addMotif(Motif motif) {
		this.motifs.addMotif(motif);
	}

}
