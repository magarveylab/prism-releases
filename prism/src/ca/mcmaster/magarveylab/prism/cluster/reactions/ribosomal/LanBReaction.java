package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.Atom;
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
 * Execute the reaction catalyzed by lantibiotic LanB enzmye, which converts
 * serine and threonine to dehydroalanine and dehydroaminobutyric acid,
 * respectively. For N-terminal Dha/Dhb residues, spontaneous tautomerization to
 * form the alpha-keto acid (i.e., pyruvate or 2-oxobutyrate) will also be
 * executed. Also called for the linardin serine/threonine-specific enzymes CypH
 * and CypL, and the thiopeptide serine dehydrogenase LazB.
 * 
 * @author skinnider
 *
 */
public class LanBReaction extends GenericReaction implements Reaction {

	public LanBReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.LanB,
				RibosomalDomains.CypH, RibosomalDomains.LegH,
				RibosomalDomains.CypL, RibosomalDomains.LazB,
				RibosomalDomains.GodF };
	}

	public void execute() throws NoResidueException,
			TailoringSubstrateException, ScaffoldGenerationException {
		for (int i = 0; i < plan.modules().size(); i++) {
			IAtomContainer molecule = scaffold.molecule();
			if (molecule == null)
				throw new TailoringSubstrateException(
						"Error: could not get molecule for LanB!");

			Module module = plan.get(i);
			Residue residue = scaffold.residue(module);
			if (residue == null)
				throw new TailoringSubstrateException(
						"Error: could not get residue for LanB!");
			IAtomContainer structure = residue.structure();
			if (structure == null)
				throw new TailoringSubstrateException(
						"Error: could not get residue structure for LanB!");

			IAtom ketone = residue.ketone();
			if (ketone == null)
				throw new TailoringSubstrateException(
						"Error: could not get ketone carbon for LanB!");
			IAtom ketoneOxygen = Atoms.getConnectedOxygen(ketone, structure);
			if (ketoneOxygen == null)
				throw new TailoringSubstrateException(
						"Error: could not get ketone oxygen for LanB!");
			IAtom alphaCarbon = residue.alphaCarbon();
			if (alphaCarbon == null)
				throw new TailoringSubstrateException(
						"Error: could not get alpha carbon for LanB!");

			IAtom betaCarbon = RibosomalUtil.getBetaCarbon(residue, structure);
			if (betaCarbon == null)
				throw new TailoringSubstrateException(
						"Error: could not get beta carbon for LanB!");

			IAtom oxygen = RibosomalUtil.getSerineOrThreonineHydroxyl(residue,
					structure, molecule);
			if (oxygen == null)
				throw new TailoringSubstrateException(
						"Error: could not get serine/threonine oxygen for LanB!");
			if (molecule.getBond(oxygen, betaCarbon) == null)
				throw new TailoringSubstrateException(
						"Error: no bond between serine/threonine oxygen and beta carbon!");

			// O can't have 2 bonds (e.g. methylation, esterification)
			if (molecule.getConnectedBondsCount(oxygen) == 2)
				throw new TailoringSubstrateException(
						"Error: could not dehydrate serine/threonine: oxygen has 2 bonds!");

			// set double bond
			IBond bond = molecule.getBond(alphaCarbon, betaCarbon);
			if (bond == null)
				throw new TailoringSubstrateException(
						"Error: no alpha/beta carbon bond for LanB!");
			UtilityReactions.setBondOrder(alphaCarbon, betaCarbon, molecule,
					IBond.Order.DOUBLE);

			// remove oxygen
			UtilityReactions.removeAtom(oxygen, molecule);
			UtilityReactions.removeAtom(oxygen, structure);

			// tautomerize N-terminal Dha/Dhb to pyruvate/2-oxobutyrate
			if (scaffold.indexOf(module) == 0) {
				// set single bond
				UtilityReactions.setBondOrder(alphaCarbon, betaCarbon,
						molecule, IBond.Order.SINGLE);

				// remove N
				IAtom nitrogen = Atoms.getNitrogen(structure);
				UtilityReactions.removeAtom(nitrogen, molecule);

				// add new ketone
				IAtom oxygen2 = new Atom("O");
				molecule.addAtom(oxygen2);
				UtilityReactions.addBond(alphaCarbon, oxygen2, molecule,
						IBond.Order.DOUBLE);
			}
		}
	}
}
