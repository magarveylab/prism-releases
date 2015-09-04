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
 * Execute the reaction catalyzed by product template domains of clade II
 * (tetrahydroxynaphthalene synthases).
 * 
 * @author skinnider
 *
 */
public class ProductTemplateCladeIIReaction extends GenericReaction implements Reaction {

	public ProductTemplateCladeIIReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { ThiotemplatedDomains.PRODUCT_TEMPLATE_II };
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
		
		// add bonds between C1-C10, C2-C7
		UtilityReactions.addBond(c1, c10, molecule, IBond.Order.DOUBLE);
		UtilityReactions.addBond(c2, c7, molecule);

		// reduce ketones at C1, C3, C5, C9
		UtilityReactions.reduceKetone(c1, molecule);
		UtilityReactions.reduceKetone(c3, molecule);
		UtilityReactions.reduceKetone(c5, molecule);
		UtilityReactions.reduceKetone(c9, molecule);

		// remove ketone at C7
		UtilityReactions.removeKetone(c7, molecule);
		
		// set double bonds between C2-C3, C4-C5, C6-C7, C8-C9
		UtilityReactions.setBondOrder(c2, c3, molecule, IBond.Order.DOUBLE);
		UtilityReactions.setBondOrder(c4, c5, molecule, IBond.Order.DOUBLE);
		UtilityReactions.setBondOrder(c6, c7, molecule, IBond.Order.DOUBLE);
		UtilityReactions.setBondOrder(c8, c9, molecule, IBond.Order.DOUBLE);
	}
	
}
