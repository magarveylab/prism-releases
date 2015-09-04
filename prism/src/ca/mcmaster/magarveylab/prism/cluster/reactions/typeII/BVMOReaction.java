package ca.mcmaster.magarveylab.prism.cluster.reactions.typeII;

import org.openscience.cdk.Atom;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.TypeIIPolyketideDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

public class BVMOReaction extends GenericReaction implements Reaction {

	public BVMOReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] {
				TypeIIPolyketideDomains.BVMO };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		
		Residue c1_2 = scaffold.residue(plan.get(0));
		Residue c3_4 = scaffold.residue(plan.get(1));
		Residue c15_16 = scaffold.residue(plan.get(2));
		Residue c17_18 = scaffold.residue(plan.get(3));
		
		IAtom c1 = c1_2.ketone();
		IAtom c2 = c1_2.alphaCarbon();
		IAtom c3 = c3_4.ketone();
		IAtom c4 = c3_4.alphaCarbon();
		IAtom c15 = c15_16.ketone();
		IAtom c16 = c15_16.alphaCarbon();
		IAtom c17 = c17_18.ketone();
		IAtom c18 = c17_18.alphaCarbon();
		
		// remove C16-C17, C17-C18 bonds; remove C17 and attached oxygen
		for (IBond bond : molecule.getConnectedBondsList(c17)) {
			IAtom connectedAtom = bond.getConnectedAtom(c17);
			if (connectedAtom.getSymbol().equals("O")) 
				molecule.removeAtom(connectedAtom);
			molecule.removeBond(bond);
		}
		molecule.removeAtom(c17);
		
		// add hydroxyl to C18
		IAtom oxygen = new Atom("O");
		molecule.addAtom(oxygen);
		c17_18.structure().addAtom(oxygen);
		UtilityReactions.addBond(c18, oxygen, molecule);
		
		// set single bond between C1-C2
		UtilityReactions.setBondOrder(c1, c2, molecule, IBond.Order.SINGLE);
		UtilityReactions.setBondOrder(c3, c4, molecule, IBond.Order.SINGLE);
		UtilityReactions.setBondOrder(c15, c16, molecule, IBond.Order.SINGLE);
		
		// oxidize C1, C15
		UtilityReactions.oxidizeAlcohol(c1, molecule);
		UtilityReactions.oxidizeAlcohol(c15, molecule);
	}

}
