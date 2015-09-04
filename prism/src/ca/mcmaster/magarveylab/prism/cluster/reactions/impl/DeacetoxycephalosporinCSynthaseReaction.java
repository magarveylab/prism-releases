package ca.mcmaster.magarveylab.prism.cluster.reactions.impl;

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import ca.mcmaster.magarveylab.enums.domains.BetaLactamDomains;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Atoms;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Bonds;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

public class DeacetoxycephalosporinCSynthaseReaction extends GenericReaction implements Reaction {
	
	public DeacetoxycephalosporinCSynthaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { BetaLactamDomains.DOACS };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		
		Residue valineResidue = scaffold.residue(plan.get(1));
		IAtomContainer valine = valineResidue.structure();
		
		// remove methyl
		IAtom methyl = Atoms.getMethylCarbon(valine);
		IAtom carbon = Atoms.getConnectedCarbon(methyl, valine);
		IBond bond = Bonds.getConnectedCarbonBond(valine, methyl);
		molecule.removeBond(bond);
		molecule.removeAtom(methyl);
		
		// add double bond
		IAtom c2 = null;
		for (IAtom atom : molecule.getConnectedAtomsList(carbon)) 
			if (atom.getSymbol().equals("C") && molecule.getConnectedBondsCount(atom) == 3)
				c2 = atom;
		IBond bond2 = molecule.getBond(carbon, c2);
		bond2.setOrder(IBond.Order.DOUBLE);
		
		// remove sulfur-carbon bond
		IAtom sulfur = Atoms.getConnectedSulfur(carbon, molecule);
		IBond bond3 = molecule.getBond(sulfur, carbon);
		molecule.removeBond(bond3);
		
		// add another methylene group
		IAtomContainer methane = new AtomContainer();
		IAtom methylene = new Atom("C");
		methane.addAtom(methylene);
		molecule.add(methane);
		valine.add(methane);
		
		UtilityReactions.addBond(sulfur, methylene, molecule);
		UtilityReactions.addBond(carbon, methylene, molecule);

		// create hydroxyl group
		IAtomContainer hydroxyl = new AtomContainer();
		IAtom oxygen = new Atom("O");
		hydroxyl.addAtom(oxygen);
		molecule.add(hydroxyl);
		valine.add(hydroxyl);
		
		// hydroxylate
		for (IAtom atom : molecule.getConnectedAtomsList(carbon))
			if (atom.getSymbol().equals("C") && molecule.getConnectedBondsCount(atom) == 1)
				UtilityReactions.addBond(atom, oxygen, molecule);
	}
	
}
