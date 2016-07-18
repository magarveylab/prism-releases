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
 * Execute the biosynthesis of the unique thioviridamide N-terminal serine-derived unit.
 * 
 * @author skinnider
 *
 */
public class TvaDReaction extends GenericReaction implements Reaction {

	public TvaDReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.TvaD };
	}

	public void execute() throws NoResidueException,
			TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();

		Module module = plan.get(0);
		Residue residue = scaffold.residue(module);
		IAtomContainer structure = residue.structure();
		
		IAtom nitrogen = residue.nitrogen();
		IAtom alphaCarbon = residue.alphaCarbon();
		
		// get beta carbon and hydroxyl
		IAtom betaCarbon = RibosomalUtil.getBetaCarbon(residue, structure);
		if (betaCarbon == null)
			throw new ScaffoldGenerationException("Couldn't get beta carbon for N-terminal serine!");
		IAtom hydroxyl = RibosomalUtil.getSerineOrThreonineHydroxyl(residue, structure, molecule);
		if (hydroxyl == null)
			throw new ScaffoldGenerationException("Couldn't get hydroxyl group for N-terminal serine!");
		
		// remove nitrogen
		UtilityReactions.removeAtom(nitrogen, molecule);
		
		// hydroxylate and methylate alpha carbon
		String oh = "OI";
		String ch3 = "CI";
		UtilityReactions.functionalize(oh, alphaCarbon, molecule);
		UtilityReactions.functionalize(ch3, alphaCarbon, molecule);
		
		// remove hydroxyl
		UtilityReactions.removeAtom(hydroxyl, molecule);
		
		// acetylate 
		String acetyl = "IC(C)=O";
		UtilityReactions.functionalize(acetyl, betaCarbon, molecule);
	}
	
}
