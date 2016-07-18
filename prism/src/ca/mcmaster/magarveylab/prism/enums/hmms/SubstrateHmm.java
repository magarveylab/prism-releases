package ca.mcmaster.magarveylab.prism.enums.hmms;

import ca.mcmaster.magarveylab.enums.interfaces.SubstrateType;

/**
 * A substrate activated by a biosynthetic domain which is detected in PRISM by
 * a substrate-specific HMM.
 * 
 * @author skinnider
 *
 */
public interface SubstrateHmm extends SubstrateType {

	/**
	 * Get the name of the hidden Markov model (.hmm) file associated with this
	 * domain substrate.
	 * 
	 * @return the name of the .hmm file
	 */
	public String hmm();

}
