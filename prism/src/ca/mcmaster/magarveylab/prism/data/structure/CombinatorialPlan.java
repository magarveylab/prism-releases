package ca.mcmaster.magarveylab.prism.data.structure;

import java.util.List;

import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.prism.data.Cyclization;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;

/**
 * A combinatorial plan for the generation of a natural product scaffold.
 * @author skinnider
 *
 */
public class CombinatorialPlan {
	
	private Cyclization cyclization;
	private List<ReactionPlan> reactions;
	private List<Module> permutation;
	
	/**
	 * Instantiate a new combinatorial scheme.
	 * @param reactions		list of reactions
	 * @param cyclization 	cyclization plan
	 * @param permutation	module permutation 
	 */
	public CombinatorialPlan(List<ReactionPlan> reactions, Cyclization cyclization, 
			List<Module> permutation) {
		this.reactions = reactions;
		this.cyclization = cyclization;
		this.permutation = permutation;
	}
	
	/**
	 * Get the cyclization pattern for this combinatorial plan.
	 * @return	combinatorial cyclization pattern
	 */
	public Cyclization cyclization() {
		return cyclization;
	}
	
	/**
	 * Get the tailoring reactions for this combinatorial scheme.
	 * @return	all tailoring reactions in this combinatorial scheme
	 */
	public List<ReactionPlan> reactions() {
		return reactions;
	}
	
	/**
	 * Get the module permutation associated with this combinatorial scheme.
	 * @return	cluster module permutation
	 */
	public List<Module> permutation() {
		return permutation;
	}
	
	/**
	 * Determine whether this combinatorial scheme includes a tailoring reaction domain.
	 * @param type	tailoring reaction to check for
	 * @return		true if this scheme contains that domain
	 */
	public boolean contains(ThiotemplatedDomains type) {
		boolean flag = false;
		for (ReactionPlan t : reactions) 
			if (t.type() == type)
				flag = true;
		return flag;
	}
	
}
