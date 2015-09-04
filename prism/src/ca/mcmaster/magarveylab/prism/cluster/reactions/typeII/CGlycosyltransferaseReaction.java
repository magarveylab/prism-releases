package ca.mcmaster.magarveylab.prism.cluster.reactions.typeII;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Atoms;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;
import ca.mcmaster.magarveylab.prism.util.SmilesIO;

public class CGlycosyltransferaseReaction  extends GenericReaction implements Reaction {

	public CGlycosyltransferaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { TailoringDomains.C_GLYCOSYLTRANSFERASE };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		Module module = plan.get(0);
		
		Domain domain = plan.domain();
		if (domain.sugar() == null)
			throw new TailoringSubstrateException("Error: could not execute C-glycosyltransferase without sugar!");
		
		Residue residue = scaffold.residue(module);
		IAtomContainer structure = residue.structure();
		
		// alpha carbon w/ only two bonds, bond order sum = 3
		IAtom alphaCarbon = residue.alphaCarbon();
		if (molecule.getBondOrderSum(alphaCarbon) != 3 || molecule.getConnectedBondsCount(alphaCarbon) != 2)
			throw new TailoringSubstrateException("Error: could not execute C-glycosyltransferase on alpha carbon"
					+ "with incorrect bonding!");

		IAtomContainer sugar = null;
		try {
			sugar = SmilesIO.molecule(domain.sugar().smiles());
		} catch (Exception e) {
			throw new ScaffoldGenerationException("Error: could not build sugar SMILES");
		}
		IAtom iodine = Atoms.getIodine(sugar);
		IAtom connectionAtom = Atoms.getConnectedCarbon(iodine, sugar);	
		UtilityReactions.removeIodine(sugar);
		
		// add to scaffold & add bond
		try {
			molecule.add(sugar);
			structure.add(sugar);
			UtilityReactions.addBond(alphaCarbon, connectionAtom, molecule);
		} catch (ArrayIndexOutOfBoundsException e) {
			throw new ScaffoldGenerationException("Error: could not add bond between alpha carbon " 
					+ molecule.getAtomNumber(alphaCarbon) + " and sugar atom " + molecule.getAtomNumber(connectionAtom));
		}
	}

}
