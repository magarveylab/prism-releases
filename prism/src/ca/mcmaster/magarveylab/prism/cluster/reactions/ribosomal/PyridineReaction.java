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
import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Atoms;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

/**
 * Execute pyridine ring formation, as catalyzed by the thiopeptide
 * cycloaddition enzyme LazC.
 * 
 * @author skinnider
 *
 */
public class PyridineReaction extends GenericReaction implements Reaction {
	
	public PyridineReaction(ReactionPlan plan, Scaffold scaffold,
			Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.LazC,
				RibosomalDomains.LazC_b };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module m1 = plan.get(0);
		Residue r1 = scaffold.residue(m1);
		if (r1 == null) 
			throw new TailoringSubstrateException("Error: could not get first Dha residue for LazC!");
		IAtomContainer s1 = r1.structure();
		
		Module m2 = plan.get(1);
		Residue r2 = scaffold.residue(m2);
		if (r2 == null) 
			throw new TailoringSubstrateException("Error: could not get second Dha residue for LazC!");
		IAtomContainer s2 = r2.structure();
		
		Module m3 = plan.get(2);
		Residue r3 = scaffold.residue(m3);
		if (r3 == null) 
			throw new TailoringSubstrateException("Error: could not get final residue for LazC!");
		IAtomContainer s3 = r3.structure();

		if (!RibosomalUtil.isDehydrated(r1, s1, molecule)) 
			throw new TailoringSubstrateException("Error: could not execute pyridine "
					+ "formation with first non-dehydrated serine!");
		if (!RibosomalUtil.isDehydrated(r2, s2, molecule)) 
			throw new TailoringSubstrateException("Error: could not execute pyridine "
					+ "formation with second non-dehydrated serine!");
		if (!RibosomalUtil.hasKetoneOxygen(r3, molecule)) 
			throw new TailoringSubstrateException("Error: could not execute pyridine "
					+ "formation: third residue lacks ketone oxygen!");
		
		IAtom betaCarbon1 = RibosomalUtil.getBetaCarbon(r1, s1);
		IAtom betaCarbon2 = RibosomalUtil.getBetaCarbon(r2, s2);
		if (molecule.getConnectedBondsCount(betaCarbon1) != 1)
			throw new TailoringSubstrateException("Error: could not execute pyridine "
					+ "formation: first Dha has > 1 bond!");
		if (molecule.getConnectedBondsCount(betaCarbon2) != 1)
			throw new TailoringSubstrateException("Error: could not execute pyridine "
					+ "formation: second Dha has > 1 bond!");

		// third residue can't be a thiazole
		IAtom sulfur = Atoms.getSulfur(s3);
		if (sulfur != null && molecule.getConnectedBondsCount(sulfur) != 1)
			throw new TailoringSubstrateException("Error: could not execute pyridine "
					+ "formation: third residue is a cyclized cysteine!");
		
		if (molecule.getBond(betaCarbon2, betaCarbon1) != null)
			throw new TailoringSubstrateException("Error: could not execute pyridine "
					+ "formation: bond already exists!");

		// add bond between two Dha beta carbons
		UtilityReactions.addBond(betaCarbon1, betaCarbon2, molecule);
		
		if (scaffold.indexOf(m1) > 0) {
			// set bond order between first Dha alpha/beta carbons for non-pyridines
			UtilityReactions.setBondOrder(betaCarbon1, r1.alphaCarbon(), molecule, IBond.Order.SINGLE);
		} else {
			// remove the ketone from tautomerized Dha
			for (IAtom atom : molecule.getConnectedAtomsList(r1.alphaCarbon()))
				if (atom.getSymbol().equals("O")
						&& molecule.getBond(r1.alphaCarbon(), atom).getOrder() == IBond.Order.DOUBLE)
					UtilityReactions.removeAtom(atom, molecule);
			
			// set bond order between first Dha alpha/beta carbons for pyridines 
			UtilityReactions.setBondOrder(betaCarbon1, r1.alphaCarbon(), molecule, IBond.Order.DOUBLE);
		}
		
		// add bond between first Dha alpha carbon and 3rd residue ketone carbon
		UtilityReactions.addBond(r1.alphaCarbon(), r3.ketone(), molecule);

		// remove ketone oxygen
		UtilityReactions.removeOxygen(r3.ketone(), molecule);

		// if tetrahydropyridine-forming...
		if (plan.domain().type() == RibosomalDomains.LazC_b) {
			// reduce second Dha
			UtilityReactions.setBondOrder(r2.alphaCarbon(), betaCarbon2,
					molecule, IBond.Order.SINGLE);
			// set bond between second Dha alpha carbon and nitrogen
			UtilityReactions.setBondOrder(r2.alphaCarbon(), r2.nitrogen(),
					molecule, IBond.Order.DOUBLE);
		} else {
			// if dihydropyridine forming, set bond order between 3rd residue
			// ketone carbon and 2nd residue nitrogen
			UtilityReactions.setBondOrder(r3.ketone(), r2.nitrogen(), molecule,
					IBond.Order.DOUBLE);
		}
	}
	
}
