package ca.mcmaster.magarveylab.prism.data;

import ca.mcmaster.magarveylab.enums.CyclizationPatterns;

/**
 * A potential macrocyclization pattern for a predicted natural product
 * scaffold.
 * 
 * @author skinnider
 *
 */
public class Cyclization {
	
	private CyclizationPatterns type;
	private Module nTerminus;

	/**
	 * Instantiate a new macrocyclization pattern.
	 * 
	 * @param nTerminus
	 *            the N-terminal module
	 * @param type
	 *            the type of cyclization
	 */
	public Cyclization(Module nTerminus, CyclizationPatterns type) {
		this.nTerminus = nTerminus;
		this.type = type;
	}
	
	/**
	 * Get the type of cyclization pattern (e.g., cyclic imine or macrolactam)
	 * associated with this macrocyclization.
	 * 
	 * @return the cyclization pattern type
	 */
	public CyclizationPatterns type() {
		return type;
	}

	/**
	 * Get the N-terminal module for macrocyclization.
	 * 
	 * @return the N-terminal cyclization module
	 */
	public Module terminus() {
		return nTerminus;
	}

}
