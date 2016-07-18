package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.Atom;
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
 * Execute beta-hydroxylation of an amino acid substrate, as catalyzed by the
 * YmE, TvaJ, TclD, and CinX enzymes, or the iterative proteusin enzyme PoyI.
 * 
 * @author skinnider
 *
 */
public class BetaHydroxylaseReaction extends GenericReaction implements
		Reaction {

	public BetaHydroxylaseReaction(ReactionPlan plan, Scaffold scaffold,
			Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.CinX,
				RibosomalDomains.YmE, RibosomalDomains.TvaJ,
				RibosomalDomains.TclD, RibosomalDomains.BerH,
				RibosomalDomains.TpdQ, RibosomalDomains.PoyI };
	}

	public void execute() throws NoResidueException,
			TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();

		for (int i = 0; i < plan.modules().size(); i++) {
			Module module = plan.get(i);
			Residue residue = scaffold.residue(module);
			if (residue == null)
				throw new TailoringSubstrateException(
						"Error: could not get residue for beta-hydroxylation!");
			IAtomContainer structure = residue.structure();

			IAtom betaCarbon = RibosomalUtil.getBetaCarbon(residue, structure);

			// make sure there isn't already a beta-hydroxylation
			boolean hasBetaHydroxylation = false;
			for (IAtom atom : molecule.getConnectedAtomsList(betaCarbon))
				if (atom.getSymbol().equals("O"))
					hasBetaHydroxylation = true;
			if (hasBetaHydroxylation)
				throw new TailoringSubstrateException(
						"Error: could not beta-hydroxylate a residue "
						+ "which is already beta-hydroxylated!");
			
			// make sure a hydroxyl group can be added 
			if (molecule.getConnectedBondsCount(betaCarbon) >= 4)
				throw new TailoringSubstrateException("Error: could not "
						+ "beta-hydroxylate beta carbon with 4 or more bonds!");

			IAtom oxygen = new Atom("O");
			structure.addAtom(oxygen);
			molecule.addAtom(oxygen);

			UtilityReactions.addBond(betaCarbon, oxygen, molecule);
			oxygen.setImplicitHydrogenCount(1);
		}
	}

}
