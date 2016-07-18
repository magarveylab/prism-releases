package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

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
 * Execute C-terminal amide formation via methylesterification, as catalyzed by
 * TsrC. Requires TsrB (methylesterase) reaction to have already been executed.
 * 
 * @author skinnider
 *
 */
public class TsrCReaction extends GenericReaction implements Reaction {

	public TsrCReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.TsrC };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module module = plan.get(0);
		Residue residue = scaffold.residue(module);
		if (residue == null)
			throw new TailoringSubstrateException("Could not get C-terminal residue for amide formation!");
		
		// get carboxyl alcohol
		IAtom ketone = residue.ketone();
		IAtom carboxylAlcohol = Atoms.getCarboxylAlcohol(ketone, molecule);
		if (carboxylAlcohol == null)
			throw new TailoringSubstrateException("Could not get C-terminal carboxyl alcohol for amide formation!");
		if (molecule.getConnectedBondsCount(carboxylAlcohol) != 2)
			throw new TailoringSubstrateException("Carboxyl alcohol is not methylesterified!");
		
		// get methyl group
		IAtom methyl = null;
		for (IAtom atom : molecule.getConnectedAtomsList(carboxylAlcohol))
			if (atom.getSymbol().equals("C") && atom != ketone)
				methyl = atom;
		if (methyl == null)
			throw new TailoringSubstrateException("Could not get methyl ester for amide formation!");

		// remove both
		UtilityReactions.removeAtom(methyl, molecule);
		UtilityReactions.removeAtom(carboxylAlcohol, molecule);
		
		// add NH2
		String smiles = "NI";
		UtilityReactions.functionalize(smiles, ketone, molecule);
	}
	
}
