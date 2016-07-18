package ca.mcmaster.magarveylab.prism.data.reactions;

import ca.mcmaster.magarveylab.enums.ReactionPriorities;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

/**
 * Data package for a modification domain and its combinatorial substrate. 
 * @author skinnider
 *
 */
public class ReactionPlan {
	
	private ReactionPriorities rp;
	private Domain domain;
	private SubstrateSet modules;
	private String smiles;

	/**
	 * Instantiate a new tailoring domain reaction plan. 
	 * @param domain		tailoring domain
	 * @param modules		biosynthetic substrate modules
	 */
	public ReactionPlan(Domain domain, SubstrateSet modules, ReactionPriorities priority) {
		this.rp = priority;
		this.domain = domain;
		this.modules = modules;
	}
	
	/**
	 * Deep-copy a reaction plan. 
	 * @param reactionPlan	the plan to deep copy 
	 */
	public ReactionPlan(ReactionPlan reactionPlan) {
		this.rp = reactionPlan.reaction();
		this.domain = reactionPlan.domain();
		this.modules = reactionPlan.modules();
		this.smiles = reactionPlan.getSmiles();
	}
	
	public ReactionPriorities reaction() {
		return rp;
	}
	
	public int priority() {
		return rp.priority();
	}
	
	public Domain domain() {
		return domain;
	}
	
	public DomainType type() {
		return domain.type();
	}
	
	public SubstrateSet modules() {
		return modules;
	}
	
	public int size() {
		return modules.size();
	}
	
	public Module get(int n) throws TailoringSubstrateException {
		return modules.get(n);
	}

	public String getSmiles() {
		return smiles;
	}

	public void setSmiles(String smiles) {
		this.smiles = smiles;
	}
		
}
