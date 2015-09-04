package ca.mcmaster.magarveylab.prism.cluster.reactions.impl;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.enums.domains.TypeIIPolyketideDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Atoms;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.BondFormationException;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

/**
 * Execute O-methylation as catalyzed by thiotemplated system
 * O-methyltransferases.
 * 
 * @author skinnider
 *
 */
public class OMethyltransferaseReaction extends GenericReaction implements Reaction {

	public OMethyltransferaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] {
				ThiotemplatedDomains.O_METHYLTRANSFERASE,
				TypeIIPolyketideDomains.CARBOXY_OMT,
				TypeIIPolyketideDomains.C11OMT_1,
				TypeIIPolyketideDomains.C11OMT_2 };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, 
			ScaffoldGenerationException, BondFormationException {
		IAtomContainer molecule = scaffold.molecule();
		Residue residue = scaffold.residue(plan.get(0));
		
		IAtom oxygen = Atoms.getHydroxyl(residue.structure(), molecule);
		if (oxygen == null)
			throw new ScaffoldGenerationException("Error: could not find oxygen atom for O-methyltransferase!");
		
		String smiles = "CI";
		UtilityReactions.functionalize(smiles, oxygen, molecule);
	}

}
