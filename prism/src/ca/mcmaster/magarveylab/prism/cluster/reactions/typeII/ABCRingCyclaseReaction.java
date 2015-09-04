package ca.mcmaster.magarveylab.prism.cluster.reactions.typeII;

import java.util.List;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.TypeIIPolyketideDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

/**
 * Type II polyketide cyclase which catalyzes the linear cyclization of the A,
 * B, and C rings of the aromatic polyketide scaffold.
 * 
 * @author skinnider
 *
 */
public class ABCRingCyclaseReaction extends GenericReaction implements Reaction {

	public ABCRingCyclaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] {
				TypeIIPolyketideDomains.CYCLASE_CLADE_7,
				TypeIIPolyketideDomains.CYCLASE_CLADE_8_9 };
	}

	public void execute() throws NoResidueException,
			TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();

		Residue r1 = scaffold.residue(plan.get(0));
		Residue r2 = scaffold.residue(plan.get(1));
		Residue r3 = scaffold.residue(plan.get(2));
		Residue r4 = scaffold.residue(plan.get(3));
		Residue r5 = scaffold.residue(plan.get(4));
		Residue r6 = scaffold.residue(plan.get(5));
		Residue r7 = scaffold.residue(plan.get(6));

		IAtom c1 = r1.ketone();
		IAtom c2 = r1.alphaCarbon();
		IAtom c3 = r2.ketone();
		IAtom c4 = r2.alphaCarbon();
		IAtom c5 = r3.ketone();
		IAtom c6 = r3.alphaCarbon();
		IAtom c7 = r4.ketone();
		IAtom c8 = r4.alphaCarbon();
		IAtom c9 = r5.ketone();
		IAtom c10 = r5.alphaCarbon();
		IAtom c11 = r6.ketone();
		IAtom c12 = r6.alphaCarbon();
		IAtom c13 = r7.ketone();
		IAtom c14 = r7.alphaCarbon();

		// add bond between C5-C10, C3-C12, C1-C14
		UtilityReactions.addBond(c5, c10, molecule);
		UtilityReactions.addBond(c3, c12, molecule);
		UtilityReactions.addBond(c1, c14, molecule);
		
		// remove ketone from C1, C3, C5
		UtilityReactions.removeOxygen(c1, molecule);
		UtilityReactions.removeOxygen(c3, molecule);
		UtilityReactions.removeOxygen(c5, molecule);

		// set double bond betweeen C1-C2, C3-C4, C5-6, C7-C8, C9-C10, C11-C12, C13-C14
		UtilityReactions.setBondOrder(c1, c2, molecule, IBond.Order.DOUBLE);
		UtilityReactions.setBondOrder(c3, c4, molecule, IBond.Order.DOUBLE);
		UtilityReactions.setBondOrder(c5, c6, molecule, IBond.Order.DOUBLE);
		UtilityReactions.setBondOrder(c7, c8, molecule, IBond.Order.DOUBLE);
		UtilityReactions.setBondOrder(c9, c10, molecule, IBond.Order.DOUBLE);
		UtilityReactions.setBondOrder(c11, c12, molecule, IBond.Order.DOUBLE);
		UtilityReactions.setBondOrder(c13, c14, molecule, IBond.Order.DOUBLE);

		// reduce C9, C11, C13 ketones
		UtilityReactions.reduceKetone(c9, molecule);
		UtilityReactions.reduceKetone(c11, molecule);
		UtilityReactions.reduceKetone(c13, molecule);

		// remove or reduce oxygen from C7
		List<IBond> c7Bonds = molecule.getConnectedBondsList(c7);
		for (IBond bond : c7Bonds) {
			IAtom atom = bond.getConnectedAtom(c7);
			if (atom.getSymbol().equals("O")) {
				if (bond.getOrder() == IBond.Order.SINGLE) {
					molecule.removeBond(bond);
					molecule.removeAtom(atom);
					break;
				} else if (bond.getOrder() == IBond.Order.DOUBLE) {
					bond.setOrder(IBond.Order.SINGLE);
					break;
				}
			}
		}
		
	}

}
