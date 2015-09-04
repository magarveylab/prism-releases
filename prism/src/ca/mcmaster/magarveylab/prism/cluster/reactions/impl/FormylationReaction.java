package ca.mcmaster.magarveylab.prism.cluster.reactions.impl;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.enums.substrates.AdenylationSubstrates;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

/**
 * Execute N-formylation of a natural product scaffold.
 * 
 * @author skinnider
 *
 */
public class FormylationReaction extends GenericReaction implements Reaction {

	public FormylationReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { TailoringDomains.FORMYLTRANSFERASE };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		
		Module module = plan.get(0);
		Residue residue = scaffold.residue(module);
		IAtomContainer structure = residue.structure();

		// get formylation atom
		IAtom formylation = null;
		if (module.scaffold() != null 
				&& module.scaffold().topSubstrate() != null 
				&& module.scaffold().topSubstrate().type() != null)
			if (module.scaffold().topSubstrate().type() == AdenylationSubstrates.ORNITHINE 
					|| module.scaffold().topSubstrate().type() == AdenylationSubstrates.N5_HYDROXYORNITHINE) {
				IAtom nitrogen = residue.nitrogen();
				for (IAtom atom : structure.atoms())
					if (atom.getSymbol().equals("N") && atom != nitrogen)
						formylation = atom;
			} else {
				formylation = residue.nitrogen();
			}
		
		if (formylation == null)
			throw new TailoringSubstrateException("Could not execute formylation: "
					+ "could not get formylation nitrogen!");
				
		String smiles = "IC=O";
		UtilityReactions.functionalize(smiles, formylation, molecule);
	}

}
