package ca.mcmaster.magarveylab.prism.cluster.reactions.typeII;

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

public class TetracenomycinFourthRingCyclaseReaction extends GenericReaction implements Reaction {
	
	public TetracenomycinFourthRingCyclaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { TypeIIPolyketideDomains.CYCLASE_CLADE_5b };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		
		Residue r1 = scaffold.residue(plan.get(0));
		Residue r2 = scaffold.residue(plan.get(1));
		Residue r3 = scaffold.residue(plan.get(2));
		
		IAtom c2 = r1.alphaCarbon();
		IAtom c3 = r2.ketone();
		IAtom c4 = r2.alphaCarbon();
		IAtom c19 = r3.ketone();
		
		// add double bond between C2/C19
		UtilityReactions.addBond(c2, c19, molecule, IBond.Order.DOUBLE);
		
		// remove C19 ketone
		UtilityReactions.removeKetone(c19, molecule);
		
		// reduce C3 ketone
		UtilityReactions.reduceKetone(c3, molecule);
		
		// set double bond between C3/C4
		UtilityReactions.setBondOrder(c3, c4, molecule, IBond.Order.DOUBLE);
	}

}
