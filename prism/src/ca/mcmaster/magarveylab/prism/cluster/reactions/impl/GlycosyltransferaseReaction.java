package ca.mcmaster.magarveylab.prism.cluster.reactions.impl;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
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
import ca.mcmaster.magarveylab.prism.util.SmilesIO;

/**
 * Execute glycosylation of a natural product scaffold.
 * 
 * @author skinnider
 *
 */
public class GlycosyltransferaseReaction extends GenericReaction implements Reaction {

	public GlycosyltransferaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { TailoringDomains.GLYCOSYLTRANSFERASE };
	}

	public void execute() throws NoResidueException,
			TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		Module module = plan.get(0);
		
		String smiles = plan.getSmiles();
		if (smiles == null)
			throw new TailoringSubstrateException("Error: could not execute glycosyltransferase without sugar!");
		
		Residue residue = scaffold.residue(module);
		IAtomContainer structure = residue.structure();
		IAtom hydroxyl = Atoms.getHydroxyl(structure, molecule);
		if (hydroxyl == null) 
			throw new TailoringSubstrateException("Error: "
					+ "no hydroxyl in glycosyltransferase residue!");
		
		IAtomContainer sugar = null;
		try {
			sugar = SmilesIO.molecule(smiles);
		} catch (Exception e) {
			throw new ScaffoldGenerationException("Error: "
					+ "could not build sugar SMILES");
		}
		IAtom iodine = Atoms.getIodine(sugar);
		IAtom connectionAtom = Atoms.getConnectedCarbon(iodine, sugar);	
		UtilityReactions.removeIodine(sugar);
		
		// add to scaffold & add bond
		try {
			molecule.add(sugar);
			structure.add(sugar);
			UtilityReactions.addBond(hydroxyl, connectionAtom, molecule);
		} catch (ArrayIndexOutOfBoundsException e) {
			throw new ScaffoldGenerationException("Error: "
					+ "could not add bond between hydroxyl atom "
					+ molecule.getAtomNumber(hydroxyl) + " and sugar atom "
					+ molecule.getAtomNumber(connectionAtom));
		}
	}

}
