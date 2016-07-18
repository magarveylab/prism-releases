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
 * Execute tertiary thioether formation, as catalyzed by the cyclothiazomycin
 * enzyme CltM.
 * 
 * @author skinnider
 *
 */
public class CltMReaction extends GenericReaction implements Reaction {

	public CltMReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.CltM };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module m1 = plan.get(0);
		Module m2 = plan.get(1);
		Residue r1 = scaffold.residue(m1);
		Residue r2 = scaffold.residue(m2);
		IAtomContainer s1 = r1.structure();
		IAtomContainer s2 = r2.structure();
		
		// cysteine cannot be cyclized
		IAtom sulfur = Atoms.getSulfur(s1);
		if (molecule.getConnectedBondsCount(sulfur) != 1)
			throw new ScaffoldGenerationException("Could not execute tertiary thioether formation: "
					+ "cysteine has already been cyclized!");
		
		// ser residue must be dehydrated
		if (!RibosomalUtil.isDehydrated(r2, s2, molecule))
			throw new ScaffoldGenerationException("Could not execute tertiary thioether formation: "
					+ "serine residue has not been dehydrated!");

		// get ser alpha and beta carbons 
		IAtom alphaCarbon = r2.alphaCarbon();
		IAtom betaCarbon = RibosomalUtil.getBetaCarbon(r2, s2);
		
		// ser cannot be in a pyridine
		if (molecule.getConnectedAtomsCount(betaCarbon) > 1)
			throw new ScaffoldGenerationException("Could not execute tertiary thioether formation: "
					+ "serine residue is in a pyridine!");
		
		// reduce alpha/beta carbon bond
		UtilityReactions.setBondOrder(alphaCarbon, betaCarbon, s2, IBond.Order.SINGLE);

		// add bond between alpha carbon and sulfur
		UtilityReactions.addBond(alphaCarbon, sulfur, molecule);
	}
	
}
