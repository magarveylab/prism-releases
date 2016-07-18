package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

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
 * Execute lysine-tryptophan crosslinking in streptide biosynthesis.
 * 
 * @author skinnider
 *
 */
public class StrBReaction extends GenericReaction implements Reaction {

	public StrBReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.StrB };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module m1 = plan.get(0);
		Module m2 = plan.get(1);
		Residue r1 = scaffold.residue(m1);
		Residue r2 = scaffold.residue(m2);
		if (r1 == null)
			throw new TailoringSubstrateException("Could not get lysine residue for StrB reaction!");
		if (r2 == null)
			throw new TailoringSubstrateException("Could not get tryptophan residue for StrB reaction!");
		IAtomContainer s1 = r1.structure();
		IAtomContainer s2 = r2.structure();
		
		// get lysine beta carbon
		IAtom betaCarbon = RibosomalUtil.getBetaCarbon(r1, s1);
		
		// get tryptophan atom 
		IAtom n1 = r2.nitrogen();
		IAtom n2 = null;
		for (IAtom atom : s2.atoms())
			if (atom.getSymbol().equals("N") && atom != n1
					&& molecule.contains(atom))
				n2 = atom;
		if (n2 == null)
			throw new TailoringSubstrateException("Could not get "
					+ "tryptophan indole nitrogen for StrB reaction!");
		IAtom c1 = null;
		for (IAtom atom : molecule.getConnectedAtomsList(n2))
			if (atom.getSymbol().equals("C")
					&& molecule.getConnectedBondsCount(atom) == 3)
				c1 = atom;
		if (c1 == null)
			throw new TailoringSubstrateException("Could not get "
					+ "tryptophan indole bridge carbon for StrB reaction!");
		IAtom c2 = null;
		for (IAtom atom : molecule.getConnectedAtomsList(c1))
			if (atom.getSymbol().equals("C")
					&& molecule.getConnectedBondsCount(atom) == 2)
				c2 = atom;
		if (c2 == null)
			throw new TailoringSubstrateException("Could not get "
					+ "tryptophan cross-linking carbon for StrB reaction!");
		
		// add bond
		UtilityReactions.addBond(betaCarbon, c2, molecule);
	}
	
}
