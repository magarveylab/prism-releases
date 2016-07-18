package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.RibosomalUtil;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

/**
 * Execute lysinoalanine formation as catalyzed by cinnamycin orf 7 (cinorf7).
 * 
 * @author skinnider
 *
 */
public class Cinorf7Reaction extends GenericReaction implements Reaction {
	
	public Cinorf7Reaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.Cinorf7 };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module m1 = plan.get(0);
		Residue serine = scaffold.residue(m1);
		if (serine == null) 
			throw new TailoringSubstrateException("Error: could not get Dha residue for Cinorf7!");
		IAtomContainer s1 = serine.structure();
		
		Module m2 = plan.get(1);
		Residue lysine = scaffold.residue(m2);
		if (lysine == null) 
			throw new TailoringSubstrateException("Error: could not get Lys residue for Cinorf7!");
		IAtomContainer s2 = lysine.structure();
		
		if (!RibosomalUtil.isDehydrated(serine, s1, molecule)) 
			throw new TailoringSubstrateException("Error: could not execute lysinoalanine "
					+ "formation on non-dehydrated serine!");
		
		// get serine alpha,beta carbons
		IAtom betaCarbon = RibosomalUtil.getBetaCarbon(serine, s1);
		IAtom alphaCarbon = serine.alphaCarbon();
			
		// also can't already be reaction site
		for (IAtom atom : molecule.getConnectedAtomsList(betaCarbon))
			if (atom.getSymbol().equals("S"))
				throw new TailoringSubstrateException(
						"Error: could not form lanthionine: residue is already the site of lanthionine!");
		
		// get lysine terminal nitrogen 
		IAtom connectionNitrogen = null;
		IAtom nitrogen = lysine.nitrogen();
		for (IAtom atom : s2.atoms())
			if (atom.getSymbol().equals("N") && atom != nitrogen)
				connectionNitrogen = atom;
		
		// oxidize double bond
		UtilityReactions.setBondOrder(betaCarbon, alphaCarbon, s1, IBond.Order.SINGLE);
		
		// add bond 
		UtilityReactions.addBond(betaCarbon, connectionNitrogen, molecule);
	}
	
}
