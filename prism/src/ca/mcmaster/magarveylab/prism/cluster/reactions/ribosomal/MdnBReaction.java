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
 * Execute omega-amide bond formation between the epsilon nitrogen of a lysine
 * residue and the omega carboxyl group of a glutamate or aspartate residue, as
 * catalyzed by the microviridin ATP-grasp enzyme MdnB.
 * 
 * @author skinnider
 *
 */
public class MdnBReaction extends GenericReaction implements Reaction {

	public MdnBReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.MdnB };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module m1 = plan.get(0);
		Module m2 = plan.get(1);
		Residue r1 = scaffold.residue(m1);
		Residue r2 = scaffold.residue(m2);
		IAtomContainer s1 = r1.structure();
		IAtomContainer s2 = r2.structure();
		
		// get terminal lysine nitrogen
		IAtom nitrogen = null;
		IAtom n1 = r1.nitrogen();
		for (IAtom atom : s1.atoms()) 
			if (atom.getSymbol().equals("N") && atom != n1)
				nitrogen = atom;
		
		// must have only one bond
		if (nitrogen == null)
			throw new ScaffoldGenerationException("Could not execute omega-amide formation: "
					+ "could not get lysine nitrogen!");
		if (molecule.getConnectedBondsCount(nitrogen) != 1)
			throw new ScaffoldGenerationException("Could not execute omega-amide formation: "
					+ "lysine nitrogen has >1 bond!");
		
		// get omega carboxy group
		IAtom carboxyl = null;
		List<IAtom> carboxyls = Atoms.getAllCarboxyls(s2);
		if (carboxyls.size() != 1)
			throw new ScaffoldGenerationException("Could not execute omega-amide formation: "
					+ "more than one carboxylic acid in glutamate/aspartate residue!");
		carboxyl = carboxyls.get(0);
		
		// must not already be involved in MdnB/MdnC reaction
		for (IAtom atom : molecule.getConnectedAtomsList(carboxyl)) {
			if (atom.getSymbol().equals("N"))
				throw new ScaffoldGenerationException("Could not execute omega-amide formation: "
						+ " side chain carboxylic acid already has an amide bond!");
			if (atom.getSymbol().equals("O")
					&& molecule.getConnectedBondsCount(atom) > 1)
				throw new ScaffoldGenerationException("Could not execute omega-amide formation: "
						+ " side chain carboxylic acid already has an ester bond!");
		}

		// remove carboxyl alcohol
		UtilityReactions.removeCarboxylAlcohol(carboxyl, molecule);

		// add lactam bond
		UtilityReactions.addBond(carboxyl, nitrogen, molecule);
	}
	
}
