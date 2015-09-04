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

public class TetracyclineFourthRingCyclaseReaction extends GenericReaction implements Reaction {
	
	public TetracyclineFourthRingCyclaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { TypeIIPolyketideDomains.CYCLASE_CLADE_6b_SUBTYPE_2 };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		
		Residue r1 = scaffold.residue(plan.get(0));
		Residue r2 = scaffold.residue(plan.get(1));
		
		IAtom c1 = r1.ketone();
		IAtom c2 = r1.alphaCarbon();
		IAtom c17 = r2.ketone();
		IAtom c18 = r2.alphaCarbon();
		
		// add bond between C1-C18
		UtilityReactions.addBond(c1, c18, molecule);
		
		// remove C1 carboxyl alcohol
		UtilityReactions.removeCarboxylAlcohol(c1, molecule);
		
		// reduce C1 ketone 
		UtilityReactions.reduceKetone(c1, molecule);
		
		// reduce C17 ketone 
		UtilityReactions.reduceKetone(c17, molecule);
		
		// set double bond between C17-C18
		UtilityReactions.setBondOrder(c17, c18, molecule, IBond.Order.DOUBLE);
		
		// set double bond between C1-C2
		UtilityReactions.setBondOrder(c1, c2, molecule, IBond.Order.DOUBLE);
		
	}

}
