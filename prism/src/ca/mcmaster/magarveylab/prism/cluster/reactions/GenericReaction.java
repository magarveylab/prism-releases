package ca.mcmaster.magarveylab.prism.cluster.reactions;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;

/**
 * A generic chemical reaction on a natural product scaffold. Must be extended.
 * 
 * @author skinnider
 *
 */
public abstract class GenericReaction implements Reaction {
	
	protected DomainType[] domains;
	protected ReactionPlan plan;
	protected Scaffold scaffold;
	protected Cluster cluster;
	
	/**
	 * Instantiate a new generic reaction.	
	 * @param plan		the reaction plan for this reaction
	 * @param scaffold	the scaffold
	 * @param cluster	the parent cluster
	 */
	public GenericReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		this.plan = plan;
		this.scaffold = scaffold;
		this.cluster = cluster;
	}

	public DomainType[] domains() {
		return domains;
	}
	
}
