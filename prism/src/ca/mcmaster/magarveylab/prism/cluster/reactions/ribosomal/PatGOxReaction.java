package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

/**
 * Execute cyanobactin macrocyclization and thiazoline/oxazoline oxidation, as
 * catalyzed by the protease PatG with an integrated oxidase.
 * 
 * @author skinnider
 *
 */
public class PatGOxReaction extends GenericReaction implements Reaction {

	public PatGOxReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.PatG_ox };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		SubstrateSet s1 = new SubstrateSet();
		s1.add(plan.get(0));
		s1.add(plan.get(1));
		ReactionPlan p1 = new ReactionPlan(plan.domain(), s1, plan.reaction());
		try {
			MacrocyclizationReaction r1 = new MacrocyclizationReaction(p1, scaffold, cluster);
			r1.execute();
		} catch (ScaffoldGenerationException e) {
			System.out.println("Couldn't macrocyclize cyanobactin");
		}
		
		for (int i = 2; i < plan.modules().size() - 1; i += 2) {
			Module m1 = plan.get(i);
			Module m2 = plan.get(i+1);
			SubstrateSet s2 = new SubstrateSet(m1, m2);
			ReactionPlan p2 = new ReactionPlan(plan.domain(), s2, plan.reaction());
			AzoleReaction r2 = new AzoleReaction(p2, scaffold, cluster);
			r2.execute();
		}
	}
	
}
