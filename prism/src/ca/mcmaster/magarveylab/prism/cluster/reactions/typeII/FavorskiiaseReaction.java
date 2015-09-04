package ca.mcmaster.magarveylab.prism.cluster.reactions.typeII;

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.TypeIIPolyketideDomains;
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

public class FavorskiiaseReaction extends GenericReaction implements Reaction {

	public FavorskiiaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { TypeIIPolyketideDomains.FAVORSKIIASE };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		
		Module m1 = plan.get(0);
		Module m2 = plan.get(1);
		Module m3 = plan.get(2);
		Module m4 = plan.get(3);
		Module m5 = plan.get(4);
		Module m6 = plan.get(5);
		Module m7 = plan.get(6);
		Module m8 = plan.get(7);
		
		Residue r1 = scaffold.residue(m1);
		Residue r2 = scaffold.residue(m2);
		Residue r3 = scaffold.residue(m3);
		Residue r4 = scaffold.residue(m4);
		Residue r5 = scaffold.residue(m5);
		Residue r6 = scaffold.residue(m6);
		Residue r7 = scaffold.residue(m7);
		Residue r8 = scaffold.residue(m8);
		
		IAtom c1 = r1.ketone();
		IAtom c2 = r1.alphaCarbon();
		IAtom c3 = r2.ketone();
		IAtom c4 = r2.alphaCarbon();
		IAtom c5 = r3.ketone();
		IAtom c6 = r3.alphaCarbon();
		IAtom c7 = r4.ketone();
		IAtom c9 = r5.ketone();
		IAtom c11 = r6.ketone();
		IAtom c12 = r6.alphaCarbon();
		IAtom c13 = r7.ketone();
		IAtom c14 = r7.alphaCarbon();
		IAtom c15 = r8.ketone();
		
		// remove the C1 alcohol 
		UtilityReactions.removeKetone(c1, molecule);
		UtilityReactions.oxidizeAlcohol(c1, molecule);
		
		// reduce C3, C5, C7, C11, C15 ketone
		UtilityReactions.reduceKetone(c3, molecule);
		UtilityReactions.reduceKetone(c5, molecule);
		UtilityReactions.reduceKetone(c7, molecule);
		UtilityReactions.reduceKetone(c11, molecule);
		UtilityReactions.reduceKetone(c15, molecule);
		
		//  remove C9 ketone/hydroxyl
		UtilityReactions.removeOxygen(c9, molecule);
		
		// remove bond between C13-C14 
		UtilityReactions.removeBond(c13, c14, molecule);

		// add bond between C12-C14
		UtilityReactions.addBond(c12, c14, molecule);
		
		// add bond between C5 oxygen and C1
		IAtom c5o = Atoms.getConnectedOxygen(c5, molecule);
		UtilityReactions.addBond(c1, c5o, molecule);
		
		// add bond betweeen C6-C11, C7-C14
		UtilityReactions.addBond(c6, c11, molecule);
		UtilityReactions.addBond(c7, c14, molecule);
		
		// add hydroxyl to C12
		IAtomContainer hydroxyl = new AtomContainer();
		IAtom c12oxygen = new Atom("O");
		hydroxyl.addAtom(c12oxygen);
		molecule.add(hydroxyl);
		UtilityReactions.addBond(c12oxygen, c12, molecule);
		
		// add hydroxyl to C9
		IAtomContainer hydroxyl2 = new AtomContainer();
		IAtom c9oxygen = new Atom("O");
		hydroxyl2.addAtom(c9oxygen);
		molecule.add(hydroxyl2);
		UtilityReactions.addBond(c9oxygen, c9, molecule);
				
		// add bond between C9 oxygen and C13
		UtilityReactions.addBond(c13, c9oxygen, molecule);
		
		// set double bond between C2-C3, C4-C5, C15-C16
		UtilityReactions.setBondOrder(c2, c3, molecule, IBond.Order.DOUBLE);
		UtilityReactions.setBondOrder(c4, c5, molecule, IBond.Order.DOUBLE);
		UtilityReactions.setBondOrder(c14, c15, molecule, IBond.Order.DOUBLE);
	}

}
