package ca.mcmaster.magarveylab.prism.data;

import java.io.Serializable;

/**
 * The package of combinatorial data evaluated for a combinatorial scaffold
 * library.
 * 
 * @author skinnider
 *
 */
public class CombinatorialData implements Serializable {

	private static final long serialVersionUID = 1L;
	private int numOrfPermutations = 0;  
	private int numPropeptides = 0; 
	private int numSugars = 0;  
	private int numCyclizations = 0;   
	private int numReactions = 0; 
	private int numCombinatorialPlans = 0; 

	/**
	 * Get the number of open reading frame permutations evaluated in
	 * thiotemplated combinatorial scaffold library generation.
	 * 
	 * @return the number of open reading frame permutations
	 */
	public int getNumOrfPermutations() {
		return numOrfPermutations;
	}

	/**
	 * Set the number of open reading frame permutations evaluated in
	 * thiotemplated combinatorial scaffold library generation.
	 * 
	 * @param numOrfPermutations
	 *            the number of open reading frame permutations
	 */
	public void setNumOrfPermutations(int numOrfPermutations) {
		this.numOrfPermutations = numOrfPermutations;
	}

	/**
	 * Get the number of cyclizations evaluated in combinatorial scaffold
	 * library generation.
	 * 
	 * @return the number of cyclizations
	 */
	public int getNumCyclizations() {
		return numCyclizations;
	}

	/**
	 * Get the number of cyclizations evaluated in combinatorial scaffold
	 * library generation.
	 * 
	 * @param numCyclizations
	 *            the number of cyclizations
	 */
	public void setNumCyclizations(int numCyclizations) {
		this.numCyclizations = numCyclizations;
	}

	/**
	 * Get the number of sugars evaluated in combinatorial scaffold library
	 * generation.
	 * 
	 * @return the number of sugars
	 */
	public int getNumSugars() {
		return numSugars;
	}

	/**
	 * Set the number of sugars evaluated in combinatorial scaffold library
	 * generation.
	 * 
	 * @param sugars
	 *            the number of sugars
	 */
	public void setNumSugars(int sugars) {
		this.numSugars = sugars;
	}

	/**
	 * Get the number of tailoring reaction plans evaluated in combinatorial
	 * scaffold library generation.
	 * 
	 * @return the number of tailoring reaction plans
	 */
	public int getNumReactions() {
		return numReactions;
	}

	/**
	 * Set the number of tailoring reaction plans evaluated in combinatorial
	 * scaffold library generation.
	 * 
	 * @param numReactions
	 *            the number of tailoring reaction plans
	 */
	public void setNumReactions(int numReactions) {
		this.numReactions = numReactions;
	}

	/**
	 * Get the total number of combinatorial plans evaluated in combinatorial
	 * scaffold library generation.
	 * 
	 * @return the total number of combinatorial plans
	 */
	public int getNumCombinatorialPlans() {
		return numCombinatorialPlans;
	}

	/**
	 * Set the total number of combinatorial plans evaluated in combinatorial
	 * scaffold library generation.
	 * 
	 * @param numCombinatorialPlans
	 *            the total number of combinatorial plans
	 */
	public void setNumCombinatorialPlans(int numCombinatorialPlans) {
		this.numCombinatorialPlans = numCombinatorialPlans;
	}

	/**
	 * Get the total number of propeptides evaluated in ribosomal combinatorial
	 * scaffold library generation.
	 * 
	 * @return the total number of propeptides
	 */
	public int getNumPropeptides() {
		return numPropeptides;
	}

	/**
	 * Set the total number of propeptides evaluated in ribosomal combinatorial
	 * scaffold library generation.
	 * 
	 * @param numPropeptides
	 *            the total number of propeptides
	 */
	public void setNumPropeptides(int numPropeptides) {
		this.numPropeptides = numPropeptides;
	}

}
