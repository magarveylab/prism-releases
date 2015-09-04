package ca.mcmaster.magarveylab.prism.cluster.reactions.impl;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Atoms;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.Substrate;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;
import ca.mcmaster.magarveylab.prism.util.SmilesIO;

/**
 * The reaction catalyzed by an acyl-adenylating ligase on a natural product
 * scaffold where the ligase is known not to contribute to scaffold biosynthesis
 * (e.g. type II polyketides).
 * 
 * @author skinnider
 *
 */
public class AcylAdenylateLigaseReaction extends GenericReaction implements Reaction {

	public AcylAdenylateLigaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { ThiotemplatedDomains.ACYL_ADENYLATING };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		Module module = plan.get(0);
		
		Residue residue = scaffold.residue(module);
		if (residue == null) 
			throw new TailoringSubstrateException("Error: could not get residue for acyl-adenylating ligase hydroxyl group!");
		IAtomContainer structure = residue.structure();
		IAtom hydroxyl = Atoms.getHydroxyl(structure, molecule);
		if (hydroxyl == null) 
			throw new TailoringSubstrateException("Error: acyl-adenylating hydroxyl is null!");

		Substrate topSubstrate = plan.domain().topSubstrate();
		if (topSubstrate == null)
			throw new TailoringSubstrateException("Error: acyl-adenylating substrate is null!");
			
		IAtomContainer substrate = null;
		try {
			substrate = SmilesIO.molecule(topSubstrate.smiles());
		} catch (Exception e) {
			throw new ScaffoldGenerationException("Error encountered building acyl-adenylating substrate SMILES");
		}
		IAtom iodine = Atoms.getIodine(substrate);
		IAtom fluorine = Atoms.getFluorine(substrate);
		IAtom ketone = Atoms.getConnectedCarbon(iodine, substrate);	
		UtilityReactions.removeIodine(substrate);
		if (fluorine != null)
			UtilityReactions.removeFluorine(substrate);
		
		// add to scaffold & add bond
		try {
			molecule.add(substrate);
			structure.add(substrate);
			UtilityReactions.addBond(hydroxyl, ketone, molecule);
		} catch (ArrayIndexOutOfBoundsException e) {
			throw new ScaffoldGenerationException("Error forming bond between hydroxyl atom " + molecule.getAtomNumber(hydroxyl)
					+ " and ketone " + molecule.getAtomNumber(ketone));
		}
	}
	
}
