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
 * Execute the reaction catalyzed by product template domains of clade I (C2/C7
 * cyclization, as in orsellinic acid or zearalenone).
 * 
 * @author skinnider
 *
 */
public class ProductTemplateCladeIReaction extends GenericReaction implements Reaction {

	public ProductTemplateCladeIReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { ThiotemplatedDomains.PRODUCT_TEMPLATE_I };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		
		Residue r1 = scaffold.residue(plan.get(0));
		Residue r2 = scaffold.residue(plan.get(1));
		Residue r3 = scaffold.residue(plan.get(2));
		Residue r4 = scaffold.residue(plan.get(3));
		
		IAtom c2 = r1.alphaCarbon();
		IAtom c3 = r2.ketone();
		IAtom c4 = r2.alphaCarbon();
		IAtom c5 = r3.ketone();
		IAtom c6 = r3.alphaCarbon();
		IAtom c7 = r4.ketone();
		
		// add double bond between C2-C7
		UtilityReactions.addBond(c2, c7, molecule, IBond.Order.DOUBLE);
		
		// reduce C3, C5 ketones
		UtilityReactions.reduceKetone(c3, molecule);
		UtilityReactions.reduceKetone(c5, molecule);
		
		// remove C7 ketone
		UtilityReactions.removeKetone(c7, molecule);
		
		// set double bonds between C3-C4, C5-C6
		UtilityReactions.setBondOrder(c3, c4, molecule, IBond.Order.DOUBLE);
		UtilityReactions.setBondOrder(c5, c6, molecule, IBond.Order.DOUBLE);
	}
	
}
