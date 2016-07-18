package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.enums.interfaces.SubstrateType;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;
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
 * Execute the combination of serine/threonine dehydration and lantithione bond
 * formation catalyzed by the fused class 2 lantibiotic cyclase/dehydratase.
 * This reaction's execute method just calls the LanB and LanC reactions, except
 * in cases .where one or more additional, uncyclized serine or threonine
 * reactions are annotated as potential reaction sizes, in which case they are
 * dehydrated by the LanB reaction.
 * 
 * @author skinnider
 *
 */
public class LanMReaction extends GenericReaction implements Reaction {
	
	public LanMReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.LanM, 
				RibosomalDomains.ProcA };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		int size = plan.modules().size();
		for (int i = 0; i < size - 1; i+= 2) {
			Module m1 = plan.get(i);
			SubstrateSet s1 = new SubstrateSet(m1);
			ReactionPlan lanbPlan = new ReactionPlan(plan.domain(), s1, plan.reaction());
			LanBReaction lanb = new LanBReaction(lanbPlan, scaffold, cluster);
			lanb.execute();

			Module m2 = plan.get(i+1);
			if (m2.scaffold() != null
					&& m2.scaffold().topSubstrate() != null
					&& m2.scaffold().topSubstrate().type() != null) {
				SubstrateType type = m2.scaffold().topSubstrate().type();
				if (type == ProteinogenicAminoAcids.CYSTEINE) {
					// create lanthionine 
					SubstrateSet s2 = new SubstrateSet(m1, m2);
					ReactionPlan lancPlan = new ReactionPlan(plan.domain(), s2, plan.reaction());
					LanCReaction lanc = new LanCReaction(lancPlan, scaffold, cluster);
					lanc.execute();
				} else if (type == ProteinogenicAminoAcids.SERINE || type == ProteinogenicAminoAcids.THREONINE) {
					// dehydrate 
					SubstrateSet s2 = new SubstrateSet(m2);
					ReactionPlan lanbPlan2 = new ReactionPlan(plan.domain(), s2, plan.reaction());
					LanBReaction lanb2 = new LanBReaction(lanbPlan2, scaffold, cluster);
					lanb2.execute();
				}
			}
		}

		// get a potential last residue
		if (size % 2 != 0) {
			Module m1 = plan.get(size - 1);
			SubstrateSet s1 = new SubstrateSet(m1);
			ReactionPlan lanbPlan = new ReactionPlan(plan.domain(), s1,
					plan.reaction());
			LanBReaction lanb = new LanBReaction(lanbPlan, scaffold, cluster);
			lanb.execute();
		}
	}
	
}
