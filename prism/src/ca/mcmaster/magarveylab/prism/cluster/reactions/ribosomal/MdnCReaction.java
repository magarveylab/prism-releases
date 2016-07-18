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
 * Execute omega-ester bond formation between the omega-carboxyl group of a
 * glutamate or aspartate residue, and the hydroxyl group of a serine or
 * threonine residue, as catalyzed by the microviridin ATP-grasp enzyme MdnC.
 * 
 * @author skinnider
 *
 */
public class MdnCReaction extends GenericReaction implements Reaction {

	public MdnCReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.MdnC };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		for (int i = 0; i < plan.modules().size() - 1; i += 2) {
			Module m1 = plan.get(i); // serine/threonine
			Module m2 = plan.get(i+1); // glutamate/aspartate
			Residue r1 = scaffold.residue(m1);
			Residue r2 = scaffold.residue(m2);
			IAtomContainer s1 = r1.structure();
			IAtomContainer s2 = r2.structure();
			
			// get serine/threonine hydroxyl
			IAtom hydroxyl = null;
			IAtom ketone = r1.ketone();
			IAtom ketoneOxygen = Atoms.getConnectedOxygen(ketone, s1);
			for (IAtom atom : s1.atoms()) 
				if (atom.getSymbol().equals("O") && atom != ketoneOxygen
						&& molecule.getBond(atom, ketone) == null) // can't be carboxyl
					hydroxyl = atom;
			
			// must have only one bond
			if (hydroxyl == null)
				throw new ScaffoldGenerationException("Could not execute omega-ester formation: "
						+ "could not get serine/threonine hydroxyl!");
			if (molecule.getConnectedBondsCount(hydroxyl) != 1)
				throw new ScaffoldGenerationException("Could not execute omega-ester formation: "
						+ "serine/threonine hydroxyl has >1 bond!");

			// get omega carboxy group
			IAtom carboxyl = null;
			List<IAtom> carboxyls = Atoms.getAllCarboxyls(s2);
			if (carboxyls.size() != 1)
				throw new ScaffoldGenerationException("Could not execute omega-ester formation: "
						+ "more than one carboxylic acid in glutamate/aspartate residue!");
			carboxyl = carboxyls.get(0);

			// must not already be involved in MdnB/MdnC reaction
			for (IAtom atom : molecule.getConnectedAtomsList(carboxyl)) {
				if (atom.getSymbol().equals("N"))
					throw new ScaffoldGenerationException("Could not execute omega-ester formation: "
							+ " side chain carboxylic acid already has an amide bond!");
				if (atom.getSymbol().equals("O")
						&& molecule.getConnectedBondsCount(atom) > 1)
					throw new ScaffoldGenerationException("Could not execute omega-ester formation: "
							+ " side chain carboxylic acid already has an ester bond!");
			}
			
			// remove carboxyl alcohol
			UtilityReactions.removeCarboxylAlcohol(carboxyl, molecule);

			// add lactam bond
			UtilityReactions.addBond(carboxyl, hydroxyl, molecule);
		}

	}
	
}
