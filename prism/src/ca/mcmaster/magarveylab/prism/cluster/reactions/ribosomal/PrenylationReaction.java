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
 * Execute O- or N-prenylation as catalyzed by the cyanobactin prenyltransferase
 * PatF.
 * 
 * @author skinnider
 *
 */
public class PrenylationReaction extends GenericReaction implements Reaction {

	public PrenylationReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.PatF };
	}

	public void execute() throws NoResidueException,
			TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();

		for (int i = 0; i < plan.modules().size(); i++) {
			Module module = plan.get(i);
			Residue residue = scaffold.residue(module);
			IAtomContainer structure = residue.structure();
			
			// get prenylation atom
			IAtom nitrogen = residue.nitrogen();
			IAtom prenylation = RibosomalUtil.getSerineOrThreonineHydroxyl(
					residue, structure, molecule);
			if (prenylation == null
					&& molecule.getConnectedBondsCount(nitrogen) < 2)
				prenylation = nitrogen;
			
			if (prenylation == null)
				throw new ScaffoldGenerationException("Could not "
						+ "prenylate residue: could not get prenylation atom");
			
			// prenylate 
			String prenyl = "CC(I)(C=C)C";
			UtilityReactions.functionalize(prenyl, prenylation, molecule);
		}
	}
	
}
