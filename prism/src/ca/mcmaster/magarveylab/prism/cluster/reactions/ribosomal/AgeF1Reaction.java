package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.Atom;
import org.openscience.cdk.exception.CDKException;
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
import ca.mcmaster.magarveylab.prism.util.SmilesIO;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

/**
 * Execute tandem C-terminal O-methylation and N-terminal N-prenylation, as
 * catalyzed by the bifunctional cyanobactin enzyme AgeF1.
 * 
 * @author skinnider
 *
 */
public class AgeF1Reaction extends GenericReaction implements Reaction {

	public AgeF1Reaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.AgeF1 };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException, CDKException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module m1 = plan.get(0);
		Module m2 = plan.get(1);
		Residue r1 = scaffold.residue(m1);
		Residue r2 = scaffold.residue(m2);
		
		// prenylate the N-terminal nitrogen
		IAtom nitrogen = r1.nitrogen();
		if (molecule.getConnectedBondsCount(nitrogen) > 2)
			throw new TailoringSubstrateException("Error: could not prenylate N-terminal "
					+ "nitrogen with more than 2 bonds!");
		String prenyl = "CC(I)(C=C)C";
		UtilityReactions.functionalize(prenyl, nitrogen, molecule);
		
		System.out.println("SIMLES: " + SmilesIO.smiles(molecule));
		
		// O-methylate the carboxy alcohol 
		IAtom carboxylCarbon = r2.ketone();
		IAtom carboxylAlcohol = Atoms.getCarboxylAlcohol(carboxylCarbon, molecule);
		if (molecule.getConnectedBondsCount(carboxylAlcohol) > 1)
			throw new TailoringSubstrateException("Error: could not methylate C-terminal "
					+ "carboxy alcohol with more than 1 bond!");
		if (carboxylAlcohol == null)
			throw new TailoringSubstrateException("Error: could not methylate C-terminal "
					+ "carboxy alcohol: could not get alcohol!");
		
		IAtom methyl = new Atom("C");
		molecule.addAtom(methyl);
		System.out.println("SIMLES: " + SmilesIO.smiles(molecule));
		UtilityReactions.addBond(methyl, carboxylAlcohol, molecule);
		System.out.println("SIMLES: " + SmilesIO.smiles(molecule));
	}
	
}
