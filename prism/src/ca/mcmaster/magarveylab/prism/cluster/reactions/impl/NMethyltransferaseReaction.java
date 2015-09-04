package ca.mcmaster.magarveylab.prism.cluster.reactions.impl;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.BondFormationException;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

/**
 * Execute N-methylation as catalyzed by thiotemplated sytem
 * N-methyltransferases.
 * 
 * @author skinnider
 *
 */
public class NMethyltransferaseReaction extends GenericReaction implements Reaction {

	public NMethyltransferaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { ThiotemplatedDomains.N_METHYLTRANSFERASE };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, 
			ScaffoldGenerationException, BondFormationException {
		IAtomContainer molecule = scaffold.molecule();
		Residue residue = scaffold.residue(plan.get(0));
		IAtom nitrogen = residue.nitrogen();
		
		String smiles = "CI";
		UtilityReactions.functionalize(smiles, nitrogen, molecule);
	}

}
