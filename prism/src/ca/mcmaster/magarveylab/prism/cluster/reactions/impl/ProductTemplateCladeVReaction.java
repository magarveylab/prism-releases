package ca.mcmaster.magarveylab.prism.cluster.reactions.impl;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
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
 * Execute the reaction catalyzed by product template domains of clade V (C6-C11
 * and C4-C13 cyclization).
 * 
 * @author skinnider
 *
 */
public class ProductTemplateCladeVReaction extends GenericReaction implements Reaction {

	public ProductTemplateCladeVReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { ThiotemplatedDomains.PRODUCT_TEMPLATE_V };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		
		Residue r1 = scaffold.residue(plan.get(0));
		Residue r2 = scaffold.residue(plan.get(1));
		Residue r3 = scaffold.residue(plan.get(2));
		Residue r4 = scaffold.residue(plan.get(3));
		Residue r5 = scaffold.residue(plan.get(4));
		Residue r6 = scaffold.residue(plan.get(5));
		
		IAtom c4 = r1.alphaCarbon();
	//	IAtom c5 = r2.ketone();
		IAtom c6 = r2.alphaCarbon();
		IAtom c7 = r3.ketone();
		IAtom c8 = r3.alphaCarbon();
		IAtom c9 = r4.ketone();
		IAtom c10 = r4.alphaCarbon();
		IAtom c11 = r5.ketone();
	//	IAtom c12 = r5.alphaCarbon();
		IAtom c13 = r6.ketone();
	//	IAtom c14 = r6.alphaCarbon();
		
		// add double bonds between C4-C13, C6-C11
		UtilityReactions.addBond(c4, c13, molecule, IBond.Order.DOUBLE);
		UtilityReactions.addBond(c6, c11, molecule, IBond.Order.DOUBLE);
		
		// remove ketones from C11, C13
		UtilityReactions.removeKetone(c11, molecule);
		UtilityReactions.removeKetone(c13, molecule);

		// reduce ketones at C7, C9
		UtilityReactions.reduceKetone(c7, molecule);
		UtilityReactions.reduceKetone(c9, molecule);
		
		// set double bonds between C7-C8, C9-C10
		UtilityReactions.setBondOrder(c7, c8, molecule, IBond.Order.DOUBLE);
		UtilityReactions.setBondOrder(c9, c10, molecule, IBond.Order.DOUBLE);
	}
	
}
