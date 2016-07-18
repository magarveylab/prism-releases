package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.impl.NMethyltransferaseReaction;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

/**
 * Execute N,N-dimethylation catalyzed by the linardin enzyme CypM or the
 * LAP enzyme PznL.
 * 
 * @author skinnider
 *
 */
public class NNDimethyltransferaseReaction extends GenericReaction implements Reaction {

	public NNDimethyltransferaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.CypM, 
				RibosomalDomains.PznL };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		NMethyltransferaseReaction r1 = new NMethyltransferaseReaction(plan, scaffold, cluster);
		r1.execute();
		r1.execute();
	}
	
}
