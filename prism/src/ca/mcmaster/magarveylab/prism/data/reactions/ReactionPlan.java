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
	
	public Module get(int n) throws TailoringSubstrateException {
		return modules.get(n);
	}
		
}
