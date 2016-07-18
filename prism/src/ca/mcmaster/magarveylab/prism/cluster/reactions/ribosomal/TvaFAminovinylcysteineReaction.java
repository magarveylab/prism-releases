package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.exception.CDKException;
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
 * Execute aminovinylcysteine formation between a C-terminal cysteine and a
 * non-dehydrated serine or threonine residue in thioviridamide.
 * 
 * @author skinnider
 *
 */
public class TvaFAminovinylcysteineReaction extends GenericReaction implements
		Reaction {

	public TvaFAminovinylcysteineReaction(ReactionPlan plan, Scaffold scaffold,
			Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.TvaF };
	}

	public void execute() throws NoResidueException,
			TailoringSubstrateException, ScaffoldGenerationException,
			CDKException {
		IAtomContainer molecule = scaffold.molecule();

		Module m1 = plan.get(0);
		Module last = plan.get(1);
		Residue r1 = scaffold.residue(m1);
		Residue r2 = scaffold.residue(last);
		if (r1 == null)
			throw new TailoringSubstrateException(
					"Error: could not get first residue for TvaF!");
		if (r2 == null)
			throw new TailoringSubstrateException(
					"Error: could not get second residue for TvaF!");
		IAtomContainer s1 = r1.structure();
		IAtomContainer s2 = r2.structure();

		// get sulfur
		IAtom sulfur = Atoms.getSulfur(s2);
		if (sulfur == null)
			throw new TailoringSubstrateException(
					"Error: could not get cysteine sulfur atom for TvaF!");

		// must not be cyclized
		if (molecule.getConnectedBondsCount(sulfur) != 1)
			throw new TailoringSubstrateException(
					"Error: cysteine sulfur is already cyclized!");

		// decarboxylate C-terminal Cys
		IAtom ketone2 = r2.ketone();
		IAtom alphaCarbon2 = r2.alphaCarbon();
		UtilityReactions.decarboxylate(ketone2, alphaCarbon2, molecule);

		// form double bond
		IAtom betaCarbon2 = Atoms.getConnectedCarbon(sulfur, s2);
		UtilityReactions.setBondOrder(alphaCarbon2, betaCarbon2, s2,
				IBond.Order.DOUBLE);

		// get beta carbon
		IAtom betaCarbon = RibosomalUtil.getBetaCarbon(r1, s1);

		// remove serine/threonine hydroxyl
		IAtom hydroxyl = RibosomalUtil.getSerineOrThreonineHydroxyl(r1, s1,
				molecule);
		if (hydroxyl == null)
			throw new TailoringSubstrateException(
					"Error: could not get serine hydroxyl!");
		UtilityReactions.removeAtom(hydroxyl, molecule);

		// add bond
		UtilityReactions.addBond(betaCarbon, sulfur, molecule);
	}

}
