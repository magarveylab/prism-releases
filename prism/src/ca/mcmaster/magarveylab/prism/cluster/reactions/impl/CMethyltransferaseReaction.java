package ca.mcmaster.magarveylab.prism.cluster.reactions.impl;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.enums.domains.TypeIIPolyketideDomains;
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
 * Execute the C-methylation of a ketone alpha carbon as catalyzed by
 * thiotemplated system C-methyltransferases.
 * 
 * @author skinnider
 *
 */
public class CMethyltransferaseReaction extends GenericReaction implements Reaction {

	public CMethyltransferaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] {
				ThiotemplatedDomains.C_METHYLTRANSFERASE,
				TypeIIPolyketideDomains.C6CMT, 
				TypeIIPolyketideDomains.C8CMT,
				TypeIIPolyketideDomains.C10CMT };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException, BondFormationException {
		IAtomContainer molecule = scaffold.molecule();
		Residue residue = scaffold.residue(plan.get(0));
		IAtom alphaCarbon = residue.alphaCarbon();

		String smiles = "CI";
		UtilityReactions.functionalize(smiles, alphaCarbon, molecule);
	}

}
