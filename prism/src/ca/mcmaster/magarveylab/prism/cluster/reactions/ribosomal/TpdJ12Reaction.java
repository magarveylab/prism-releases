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
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

/**
 * Execute isoleucine dihydroxylation or phenylalanine beta-hydroxylation, as
 * catalyzed by the TpdJ1 and TpdJ2 enzymes.
 * 
 * @author skinnider
 *
 */
public class TpdJ12Reaction extends GenericReaction implements Reaction {
	
	public TpdJ12Reaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.TpdJ1,
				RibosomalDomains.TpdJ2 };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		for (Module module : plan.modules().getAllModules()) {
			if (module.scaffold() == null
					|| module.scaffold().topSubstrate() == null
					|| module.scaffold().topSubstrate().type() == null)
				throw new TailoringSubstrateException("Error: could not get residue for TpdJ1/TpdJ2 reaction!");

			try {
				SubstrateType aa = module.scaffold().topSubstrate().type();
				if (aa == ProteinogenicAminoAcids.ISOLEUCINE) {
					EpoxidaseReaction reaction = new EpoxidaseReaction(plan, scaffold, cluster);
					reaction.execute();
				} else if (aa == ProteinogenicAminoAcids.PHENYLALANINE) {
					BetaHydroxylaseReaction reaction = new BetaHydroxylaseReaction(plan, scaffold, cluster);
					reaction.execute();
				}
			} catch (NullPointerException e) {
				throw new TailoringSubstrateException(e.getMessage());
			}
		}
	}

}
