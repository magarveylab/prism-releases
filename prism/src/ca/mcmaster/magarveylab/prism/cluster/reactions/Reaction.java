package ca.mcmaster.magarveylab.prism.cluster.reactions;

import org.openscience.cdk.exception.CDKException;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.prism.util.exception.BondFormationException;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

/**
 * A reaction which modifies a natural product scaffold.
 * 
 * @author skinnider
 *
 */
public interface Reaction {

	/**
	 * Get the enzymatic domains which catalyze this reaction.
	 * 
	 * @return enzyme domains associated with this reaction
	 */
	public DomainType[] domains();

	/**
	 * Execute this reaction.
	 * 
	 * @throws NoResidueException
	 * @throws TailoringSubstrateException
	 * @throws ScaffoldGenerationException
	 * @throws BondFormationException 
	 * @throws CDKException 
	 */
	public void execute() throws NoResidueException,
			TailoringSubstrateException, ScaffoldGenerationException,
			BondFormationException, CDKException;

}
