package ca.mcmaster.magarveylab.prism.cluster.scaffold;

import java.io.IOException;
import java.util.List;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import ca.mcmaster.magarveylab.enums.CyclizationPatterns;
import ca.mcmaster.magarveylab.prism.cluster.analysis.TypeIIPolyketideAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.ReactionUtil;
import ca.mcmaster.magarveylab.prism.cluster.reactions.ScaffoldReactions;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Atoms;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Cyclization;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.CombinatorialPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.Sorter;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ClassInstantiationException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;
import ca.mcmaster.magarveylab.prism.util.SmilesIO;

/**
 * Generates predicted chemical structures.
 * 
 * @author skinnider
 *
 */
public class StructureGenerator {

	/**
	 * Generate one member of a combinatorial library of scaffolds given a
	 * cluster, a module permutation to consider, and a combinatorial scheme of
	 * cyclizations and tailoring reactions.
	 * 
	 * @param cluster
	 *            cluster in question
	 * @param scheme
	 *            combinatorial scheme
	 * @return generated scaffold corresponding to cluster, modules, and scheme
	 * @throws InvalidSmilesException
	 * @throws IOException
	 * @throws NoResidueException
	 * @throws ScaffoldGenerationException
	 * @throws TailoringSubstrateException
	 * @throws ClassInstantiationException
	 */
	public static Scaffold generateScaffold(Cluster cluster, CombinatorialPlan scheme) 
			throws CDKException, IOException, NoResidueException, ScaffoldGenerationException, 
			ClassInstantiationException, TailoringSubstrateException {
		Scaffold scaffold = null;
		List<Module> modules = scheme.permutation();
		Cyclization cyclization = scheme.cyclization();

		if (TypeIIPolyketideAnalyzer.isTypeIIPolyketideCluster(cluster)) {
			scaffold = TypeIIPolyketideGenerator.generateTypeIIPolyketideScaffold(modules, cluster);
		} else {
			scaffold = generateLinearScaffold(modules);
		}

		scaffold = executeReactionPlan(cluster, scaffold, scheme);
		scaffold = cyclize(scaffold, cyclization, cluster);

		return scaffold;
	}

	/**
	 * Generate a natural product scaffold from a list of biosynthetic modules.	
	 * @param modules	list of modules from which to construct scaffold
	 * @return			the scaffold 
	 * @throws IOException
	 * @throws NoResidueException 
	 * @throws ScaffoldGenerationException 
	 * @throws CDKException 
	 */
	public static Scaffold generateLinearScaffold(List<Module> modules) 
			throws IOException, NoResidueException, ScaffoldGenerationException, CDKException {
		Scaffold scaffold = null;
		int size = modules.size();
		System.out.println("[StructureGenerator] Generating new scaffold with " + size + " modules");

		int i = 0;
		while (i < size) {
			Module module = modules.get(i);
			if (scaffold == null) {
				scaffold = ScaffoldReactions.createScaffold(module);
			} else {
				ScaffoldReactions.extendScaffold(module, scaffold);
			}
			i++;
		}
		ScaffoldReactions.finishScaffold(scaffold);

		if (scaffold == null)
			throw new ScaffoldGenerationException("Could not generate linear scaffold");

		return scaffold;
	}

	/**
	 * Cyclize a generated natural product scaffold.	
	 * @param scaffold		the constructed natural product scaffold
	 * @param cyclization	the cyclization module	
	 * @return				the cyclized module, or null if cyclization failed 
	 * @throws NoResidueException
	 * @throws ScaffoldGenerationException
	 * @throws CDKException 
	 */
	public static Scaffold cyclize(Scaffold scaffold, Cyclization cyclization, Cluster cluster) 
			throws NoResidueException, ScaffoldGenerationException, CDKException {
		if (cyclization == null)
			return scaffold;
		if (cyclization.type() == CyclizationPatterns.LINEAR)
			return scaffold;
		
		// remove -COOH hydroxyl 
		IAtomContainer molecule = scaffold.molecule();
		Residue last = scaffold.getLastResidue();
		IAtom ketone = last.ketone();
		if (ketone == null)
			System.out.println("Error: could not get ketone for last residue!");
		UtilityReactions.removeCarboxylAlcohol(ketone, molecule);
		
		if (cyclization.type() == CyclizationPatterns.LINEAR_ALDEHYDE)
			return scaffold;

		// remove =O for imines
		if (cyclization.type() == CyclizationPatterns.IMINE) {
			IAtom oxygen = Atoms.getConnectedOxygen(ketone, molecule);
			if (molecule.getConnectedBondsCount(oxygen) > 1)
				throw new ScaffoldGenerationException("Could not generate cylic imine: "
						+ "ketone oxygen has > 1 bond!");
			UtilityReactions.removeAtom(oxygen, molecule);
		}
		
		// get residue
		Module module = cyclization.terminus();
		if (module == null)
			throw new ScaffoldGenerationException("Could not cyclize scaffold: "
					+ "no cyclization module!");
		Residue residue = scaffold.residue(module);
		if (residue == null)
			throw new ScaffoldGenerationException("Could not cyclize scaffold: "
					+ "could not get residue for cyclization module!");

		// get cyclization atom
		IAtom atom = null;
		if (cyclization.type() == CyclizationPatterns.IMINE 
				|| cyclization.type() == CyclizationPatterns.LACTAM) {
			atom = residue.nitrogen();
			if (atom == null)
				atom = Atoms.getPrimaryAmine(residue.structure(), molecule);
		} else if (cyclization.type() == CyclizationPatterns.LACTONE) {
			IAtom cyclicKetoneCarbon = residue.ketone();
			atom = Atoms.getConnectedHydroxyl(cyclicKetoneCarbon, molecule);
			if (atom == null) 
				atom = Atoms.getHydroxyl(residue.structure(), molecule);
		}
		if (atom == null)
			throw new ScaffoldGenerationException("Could not cyclize scaffold: "
					+ "could not get macrocyclization atom!");

		// form bond 
		if (cyclization.type() == CyclizationPatterns.IMINE) {
			UtilityReactions.addBond(atom, ketone, molecule, IBond.Order.DOUBLE);
			UtilityReactions.removeAlcohol(ketone, molecule);
		} else if (cyclization.type() == CyclizationPatterns.LACTAM
				|| cyclization.type() == CyclizationPatterns.LACTONE) {
			UtilityReactions.addBond(atom, ketone, molecule);
		}

		System.out.println("Generated cyclic scaffold with SMILES " + SmilesIO.smiles(molecule));
		return scaffold;
	}

	/**
	 * Execute tailoring reactions on a cluster scaffold.
	 * @param cluster	cluster in question
	 * @param scaffold	the generated cluster scaffold
	 * @param scheme	combinatorial tailoring reaction scheme
	 * @throws IOException 
	 * @throws ScaffoldGenerationException 
	 * @throws ClassInstantiationException 
	 * @throws TailoringSubstrateException 
	 * @throws NoResidueException 
	 * @throws CDKException 
	 */
	public static Scaffold executeReactionPlan(Cluster cluster, Scaffold scaffold, CombinatorialPlan scheme) 
			throws IOException, ScaffoldGenerationException, ClassInstantiationException, 
			NoResidueException, TailoringSubstrateException, CDKException {
		List<ReactionPlan> plans = scheme.reactions();
		if (plans != null) {
			// sort plans by priority
			Sorter.sortPlansByPriority(plans);
			// execute reaction plans
			for (ReactionPlan plan : plans) {
				Reaction reaction = ReactionUtil.getReaction(plan, scaffold, cluster);
				System.out.println("[StructureGenerator] Executing " + reaction.toString()
						+ "; current SMILES: " + SmilesIO.smiles(scaffold.molecule()));
				try {
					reaction.execute();
				} catch (ScaffoldGenerationException e) {
					System.out.println("In ScaffoldGenerationException");
					e.printStackTrace();
					continue;
				} catch (TailoringSubstrateException e) {
					System.out.println("In TailoringSubstrateException");
					e.printStackTrace();
					continue;
				}
				scaffold.incrementReactionCount();
			}

		}
		return scaffold;
	}	

}
