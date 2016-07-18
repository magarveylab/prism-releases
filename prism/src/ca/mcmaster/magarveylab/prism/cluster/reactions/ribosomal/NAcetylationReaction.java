package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
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
 * Execute N-acetylation, as catalyzed by the microviridin enzyme MdnD, the
 * goadsporin enzmye GodH, and the paenibacillin enzyme PaeN.
 * 
 * @author skinnider
 *
 */
public class NAcetylationReaction extends GenericReaction implements Reaction {

	public NAcetylationReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.MdnD,
				RibosomalDomains.GodH, RibosomalDomains.PaeN };
	}

	public void execute() throws NoResidueException,
			TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();

		Module module = plan.get(0);
		Residue residue = scaffold.residue(module);
		IAtom nitrogen = residue.nitrogen();

		// can't have >2 bonds
		if (molecule.getConnectedAtomsCount(nitrogen) > 2)
			throw new ScaffoldGenerationException("Couldn't acetylate "
					+ "N-terminus nitrogen with 3 or more bonds!");
		
		// acetylate 
		String acetyl = "IC(C)=O";
		UtilityReactions.functionalize(acetyl, nitrogen, molecule);
	}
	
}
