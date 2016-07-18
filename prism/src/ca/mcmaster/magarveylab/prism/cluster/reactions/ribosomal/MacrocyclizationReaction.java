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
 * Execute the macrocyclization of a ribosomal peptide.
 * 
 * @author skinnider
 *
 */
public class MacrocyclizationReaction extends GenericReaction implements Reaction {

	public MacrocyclizationReaction(ReactionPlan plan, Scaffold scaffold,
			Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.PatG,
				RibosomalDomains.DUF95, RibosomalDomains.YmF,
				RibosomalDomains.AlbE, RibosomalDomains.SkfC };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module m1 = plan.get(0);
		Residue r1 = scaffold.residue(m1);
		if (r1 == null)
			throw new TailoringSubstrateException("Error: could not get first "
					+ "residue for macrocyclization!");

		Module m2 = plan.get(1);
		Residue r2 = scaffold.residue(m2);
		if (r2 == null)
			throw new TailoringSubstrateException("Error: could not get last "
					+ "residue for macrocyclization!");

		IAtom ketone = r2.ketone();
		IAtom nitrogen = r1.nitrogen();
		
		// check cyclization atoms
		if (molecule.getConnectedBondsCount(nitrogen) > 2)
			throw new ScaffoldGenerationException("Could not macrocyclize: nitrogen atom has >2 bonds");
		IAtom carboxylAlcohol = Atoms.getCarboxylAlcohol(ketone, molecule);
		if (carboxylAlcohol == null)
			throw new ScaffoldGenerationException("Could not macrocyclize: carboxyl alcohol atom is null");
		if (molecule.getConnectedBondsCount(carboxylAlcohol) > 1)
			throw new ScaffoldGenerationException("Could not macrocyclize: carboxyl alcohol atom has >1 bond");
		
		// remove carboxylic acid
		UtilityReactions.removeAlcohol(ketone, molecule);
		
		// add bond
		UtilityReactions.addBond(ketone, nitrogen, molecule);
	}
	
}
