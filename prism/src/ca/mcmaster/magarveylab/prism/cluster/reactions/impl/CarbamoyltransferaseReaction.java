package ca.mcmaster.magarveylab.prism.cluster.reactions.impl;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Atoms;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

public class CarbamoyltransferaseReaction extends GenericReaction implements Reaction {

	public CarbamoyltransferaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { TailoringDomains.CARBAMOYLTRANSFERASE };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		Module module = plan.get(0);
		
		Residue residue = scaffold.residue(module);
		IAtomContainer structure = residue.structure();
		IAtom hydroxyl = Atoms.getHydroxyl(structure, molecule);
		if (hydroxyl == null) 
			throw new TailoringSubstrateException("Error: carbamoylation hydroxyl is null!");
		
		String smiles = "IC(N)=O";
		UtilityReactions.functionalize(smiles, hydroxyl, molecule);
	}
	
}
