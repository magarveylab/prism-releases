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
 * Type II polyketide cyclase which catalyzes the linear cyclization of the A
 * and B rings of the aromatic polyketide scaffold.
 * 
 * @author skinnider
 *
 */
public class ABRingCyclaseReaction extends GenericReaction implements Reaction {

	public ABRingCyclaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] {
				TypeIIPolyketideDomains.CYCLASE_CLADE_10,
				TypeIIPolyketideDomains.CYCLASE_CLADE_8a };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		
		Residue r1 = scaffold.residue(plan.get(0));
		Residue r2 = scaffold.residue(plan.get(1));
		Residue r3 = scaffold.residue(plan.get(2));
		Residue r4 = scaffold.residue(plan.get(3));
		Residue r5 = scaffold.residue(plan.get(4));
		
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

		// add bond between C3-C8, C1-C10
		UtilityReactions.addBond(c3, c8, molecule);
		UtilityReactions.addBond(c1, c10, molecule);

		// remove ketone from C3, C1
		UtilityReactions.removeKetone(c3, molecule);
		UtilityReactions.removeKetone(c1, molecule);

		// set double bond betweeen C1-C2, C3-C4, C5-C6, C7-C8, C9-C10
		UtilityReactions.setBondOrder(c1, c2, molecule, IBond.Order.DOUBLE);
		UtilityReactions.setBondOrder(c3, c4, molecule, IBond.Order.DOUBLE);
		UtilityReactions.setBondOrder(c5, c6, molecule, IBond.Order.DOUBLE);
		UtilityReactions.setBondOrder(c7, c8, molecule, IBond.Order.DOUBLE);
		UtilityReactions.setBondOrder(c9, c10, molecule, IBond.Order.DOUBLE);

		// reduce C7, C9 ketones
		UtilityReactions.reduceKetone(c7, molecule);
		UtilityReactions.reduceKetone(c9, molecule);
		
		// remove or reduce oxygen from C5
		List<IBond> c5Bonds = molecule.getConnectedBondsList(c5);
		for (IBond bond : c5Bonds) {
			IAtom atom = bond.getConnectedAtom(c5);
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
