package ca.mcmaster.magarveylab.prism.data;

import ca.mcmaster.magarveylab.enums.substrates.SubstrateType;
import ca.mcmaster.magarveylab.prism.genome.data.HmmSearchResultAnnotation;


/**
 * A substrate-specific region of a biosynthetic domain. 
 * @author skinnider
 *
 */
public class Substrate {
	
	protected int start;
	protected int end;
	protected double score;
	protected SubstrateType type;
		
	/**
	 * Instantiate a new substrate from a hmmsearch annotation.
	 * @param annotation	the annotation corresponding to this substrate
	 * @param type			the substrate type
	 */
	public Substrate(HmmSearchResultAnnotation annotation, SubstrateType type) {
		this.start = annotation.start();
		this.end = annotation.end();
		this.score = annotation.score();
		this.type = type;
	}

	/**
	 * Instantiate a new substrate.
	 * 
	 * @param start
	 *            residue on the parent domain at which this substrate starts
	 * @param end
	 *            residue on the parent domain at which this substrate ends
	 * @param score
	 *            score assigned to this subbstrate
	 * @param type
	 *            the substrate type
	 */
	public Substrate(int start, int end, double score, SubstrateType type) {
		this.start = start;
		this.end = end;
		this.score = score;
		this.type = type;
	}
	
	/**
	 * Instantiate a new substrate without any hidden Markov model results,
	 * setting start and end to -1 and score to 0.
	 * 
	 * @param type
	 *            the substrate type
	 */
	public Substrate(SubstrateType type) {
		this.type = type;
		this.start = -1;
		this.end = -1;
		this.score = 0.0d;
	}
	
	/**
	 * Get the SMILES of this substrate type.
	 * @return	substrate SMILES--which may be tagged with atoms indicating attachment sites!
	 */
	public String smiles() {
		return type.smiles();
	}
	
	/**
	 * Get the score of the hidden Markov model search for this substrate.
	 * 
	 * @return hidden Markov model search score
	 */
	public double score() {
		return score;
	}

	/**
	 * Get the start of this substrate region within the parent domain.
	 * 
	 * @return substrate region start, relative to parent domain
	 */
	public int start() {
		return start;
	}

	/**
	 * Get the end of this substrate region within the parent domain.
	 * 
	 * @return substrate region end, relative to parent domain
	 */
	public int end() {
		return end;
	}

	/**
	 * Get the domain-specific type of this substrate.
	 * 
	 * @return the substrate type
	 */
	public SubstrateType type() {
		return type;
	}

	/**
	 * Set the domain-specific type of this substrate.
	 * 
	 * @param type
	 *            the substrate type
	 */
	public void setType(SubstrateType type) {
		this.type = type;
	}

}