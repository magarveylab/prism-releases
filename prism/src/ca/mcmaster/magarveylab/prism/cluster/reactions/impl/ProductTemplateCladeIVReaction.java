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
 * Execute the reaction catalyzed by product template domains of clade IV (C4-C9
 * and C2-C11 cyclization). <br>
 * <br>
 * See figure 1 of Crawford et al., Nature 461, 1139 (2009).
 * 
 * @author skinnider
 *
 */
public class ProductTemplateCladeIVReaction extends GenericReaction implements Reaction {

	public ProductTemplateCladeIVReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { ThiotemplatedDomains.PRODUCT_TEMPLATE_IV };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		
		Residue r1 = scaffold.residue(plan.get(0));
		Residue r2 = scaffold.residue(plan.get(1));
		Residue r3 = scaffold.residue(plan.get(2));
		Residue r4 = scaffold.residue(plan.get(3));
		Residue r5 = scaffold.residue(plan.get(4));
		Residue r6 = scaffold.residue(plan.get(5));
		
		IAtom c2 = r1.alphaCarbon();
		IAtom c4 = r2.alphaCarbon();
		IAtom c5 = r3.ketone();
		IAtom c6 = r3.alphaCarbon();
		IAtom c7 = r4.ketone();
		IAtom c8 = r4.alphaCarbon();
		IAtom c9 = r5.ketone();
		IAtom c11 = r6.ketone();

		// add double bonds between C2-C11, C4-C9
		UtilityReactions.addBond(c2, c11, molecule, IBond.Order.DOUBLE);
		UtilityReactions.addBond(c4, c9, molecule, IBond.Order.DOUBLE);
		
		// reduce ketones at C5, C7
		UtilityReactions.reduceKetone(c5, molecule);
		UtilityReactions.reduceKetone(c7, molecule);

		// remove ketones at C9, C11
		UtilityReactions.removeKetone(c9, molecule);
		UtilityReactions.removeKetone(c11, molecule);

		// set double bond between C5-C6, C7-C8
		UtilityReactions.setBondOrder(c5, c6, molecule, IBond.Order.DOUBLE);
		UtilityReactions.setBondOrder(c7, c8, molecule, IBond.Order.DOUBLE);
	}
	
}
