package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;
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
 * Execute aminovinylcysteine formation between a C-terminal cysteine and
 * another cysteine.
 * 
 * @author skinnider
 *
 */
public class AminovinylcysteineReaction extends GenericReaction implements Reaction {
	
	public AminovinylcysteineReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.LanD,
				RibosomalDomains.TvaF };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		
		Module m1 = plan.get(0);
		Module m2 = plan.get(1);
		Residue r1 = scaffold.residue(m1);
		Residue r2 = scaffold.residue(m2);
		if (r1 == null) 
			throw new TailoringSubstrateException("Error: could not get first residue for LanD!");
		if (r2 == null) 
			throw new TailoringSubstrateException("Error: could not get second residue for LanD!");
		
		IAtomContainer s1 = r1.structure();
		IAtomContainer s2 = r2.structure();
		IAtom sulfur2 = Atoms.getSulfur(s2);
		
		if (m1.scaffold().topSubstrate().type() == ProteinogenicAminoAcids.CYSTEINE) {
			// remove residue #1 sulfur
			IAtom sulfur1 = Atoms.getSulfur(s1);
			if (sulfur1 == null)
				throw new TailoringSubstrateException("Error: could not get cysteine sulfur!");
			IAtom betaCarbon1 = Atoms.getConnectedCarbon(sulfur1, molecule);
			UtilityReactions.removeAtom(sulfur1, molecule);
			
			// add bond between sulfur #2, beta carbon #1
			UtilityReactions.addBond(sulfur2, betaCarbon1, molecule);
		} else {
			System.out.println("[AminovinylcysteineReaction] " + m1.scaffold().topSubstrate().toString());
			// must be Dha/Dhb
			if (!RibosomalUtil.isDehydrated(r1, s1, molecule))
				throw new TailoringSubstrateException(
						"Error: could not form AviCys from non-dehydrated serine/threonine!");
			
			IAtom ketone = r1.ketone();
			IAtom alphaCarbon = r1.alphaCarbon();
			IAtom betaCarbon = null;
			for (IAtom atom : s1.getConnectedAtomsList(alphaCarbon))
				if (atom.getSymbol().equals("C") && atom != ketone)
					betaCarbon = atom;
			
			// also can't already be reaction site
			for (IAtom atom : molecule.getConnectedAtomsList(betaCarbon))
				if (atom.getSymbol().equals("S"))
					throw new TailoringSubstrateException(
							"Error: could not form AviCys: residue is already the site of lanthionine!");

			// oxidize double bond
			UtilityReactions.setBondOrder(alphaCarbon, betaCarbon, molecule, IBond.Order.SINGLE);
			
			// add bond
			UtilityReactions.addBond(sulfur2, betaCarbon, molecule);
		}
		
		// decarboxylate C-terminal Cys
		IAtom ketone2 = r2.ketone();
		IAtom alphaCarbon2 = r2.alphaCarbon();
		UtilityReactions.decarboxylate(ketone2, alphaCarbon2, molecule);
		
		// form double bond
		IAtom betaCarbon2 = Atoms.getConnectedCarbon(sulfur2, s2);
		UtilityReactions.setBondOrder(alphaCarbon2, betaCarbon2, s2, IBond.Order.DOUBLE);
	}

}
