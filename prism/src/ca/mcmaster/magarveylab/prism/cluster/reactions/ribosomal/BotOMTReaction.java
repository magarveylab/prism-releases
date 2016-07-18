package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import java.util.List;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
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
 * Execute aspartate side chain O-methylation, as catalyzed by the bottromycin
 * enzyme BotOMT.
 * 
 * @author skinnider
 *
 */
public class BotOMTReaction extends GenericReaction implements Reaction {
	
	public BotOMTReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.BotOMT };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module module = plan.get(0);
		Residue residue = scaffold.residue(module);
		IAtomContainer structure = residue.structure();
		
		// get carboxyl C
		IAtom carboxyl = null;
		IAtom ketone = residue.ketone();
		List<IAtom> carboxyls = Atoms.getAllCarboxyls(structure);
		if (carboxyls.size() > 2)
			throw new TailoringSubstrateException("Couldn't execute BotOMT reaction: "
					+ "aspartate residue has >2 carboxyls!");
		for (IAtom atom : carboxyls)
			if (atom != ketone)
				carboxyl = atom;
		
		// get carboxyl alcohol
		IAtom carboxylAlcohol = Atoms.getCarboxylAlcohol(carboxyl, structure);

		String smiles = "CI";
		UtilityReactions.functionalize(smiles, carboxylAlcohol, molecule);
	}
	
}
