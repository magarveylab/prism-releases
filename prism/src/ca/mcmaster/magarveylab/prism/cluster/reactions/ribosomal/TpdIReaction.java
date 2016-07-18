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
 * Execute thiazole C-methylation, as catalyzed by the TpdI and TpdL enzymes.
 * 
 * @author skinnider
 *
 */
public class TpdIReaction extends GenericReaction implements Reaction {

	public TpdIReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.TpdI, 
				RibosomalDomains.TpdL };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module module = plan.get(0);
		Residue residue = scaffold.residue(module);
		if (residue == null)
			throw new TailoringSubstrateException("Could not get cysteine residue for C-methylation!");
		IAtomContainer structure = residue.structure();
		
		// must be cyclized
		IAtom sulfur = Atoms.getSulfur(structure);
		if (molecule.getConnectedBondsCount(sulfur) != 2) 
			throw new TailoringSubstrateException("Could not execute thiazole C-methylation: "
					+ "sulfur does not have two bonds!");
		
		// must be azole
		IAtom betaCarbon = RibosomalUtil.getBetaCarbon(residue, structure);
		IBond bond = molecule.getBond(betaCarbon, residue.alphaCarbon());
		if (bond.getOrder() != IBond.Order.DOUBLE)
			throw new TailoringSubstrateException("Could not execute thiazole C-methylation: "
					+ "alpha carbon-beta carbon bond order is not double!");

		// cannot already be C-methylated
		if (molecule.getConnectedBondsCount(betaCarbon) > 2)
			throw new TailoringSubstrateException("Could not execute thiazole C-methylation: "
					+ "beta carbon already has > 2 bonds!");
		
		UtilityReactions.functionalize("CI", betaCarbon, molecule);
	}
	
}
