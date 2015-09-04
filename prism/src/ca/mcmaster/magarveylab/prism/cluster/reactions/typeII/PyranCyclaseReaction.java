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

public class PyranCyclaseReaction extends GenericReaction implements Reaction {
	
	public PyranCyclaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { TypeIIPolyketideDomains.CYCLASE_CLADE_3 };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		
		Residue r1 = scaffold.residue(plan.get(0));
		Residue r2 = scaffold.residue(plan.get(1));
		
		IAtom c1 = r1.ketone();
		IAtom c13 = r2.ketone();
	
		// get C1 oxygen
		IAtom c1Oxygen = null;
		List<IBond> c1Bonds = molecule.getConnectedBondsList(c1);
		for (IBond bond : c1Bonds) {
			IAtom atom = bond.getConnectedAtom(c1);
			if (atom.getSymbol().equals("O") && bond.getOrder() == IBond.Order.DOUBLE)
				c1Oxygen = atom;
		}
		if (c1Oxygen == null)
			throw new ScaffoldGenerationException("Error: could not find C1 oxygen for pyran ring formation!");
		

		// add bond between C1 oxygen and C13
		UtilityReactions.addBond(c1Oxygen, c13, molecule);

		// reduce C1 oxygen
		UtilityReactions.reduceKetone(c1, molecule);

		// reduce C13 oxygen
		UtilityReactions.reduceKetone(c13, molecule);
	}

}
