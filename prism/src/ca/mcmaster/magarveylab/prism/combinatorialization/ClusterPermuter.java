package ca.mcmaster.magarveylab.prism.combinatorialization;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;

import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.prism.cluster.analysis.ClusterAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.analysis.TransATAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.analysis.TypeIIPolyketideAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.TypeIIPolyketideGenerator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.CombinatorialData;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.data.Propeptide;
import ca.mcmaster.magarveylab.enums.ClusterFamilies;
import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.prism.util.exception.FimoSearchException;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.RibosomalScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.cluster.analysis.RibosomalClusterAnalyzer;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Generate permutations for a cluster.
 * @author skinnider, Robyn Edgar
 *
 */
public class ClusterPermuter {
	
	/**
	 * Generate all permutations of biosynthetic modules within a cluster.
	 * 
	 * @param cluster
	 * @param session
	 * @return
	 * @throws IOException
	 * @throws ClassNotFoundException
	 * @throws NoResidueException
	 * @throws FimoSearchException
	 * @throws RibosomalScaffoldGenerationException
	 * @throws InterruptedException 
	 */
	public static List<List<Module>> getModulePermutations(Cluster cluster,
			Session session) throws IOException, NoResidueException,
			FimoSearchException, RibosomalScaffoldGenerationException, InterruptedException {
		List<List<Module>> modulePermutations = new ArrayList<List<Module>>();
		CombinatorialData cd = cluster.combinatorialData();

		if (ClusterAnalyzer.isTransAT(cluster)) {
			modulePermutations = TransATAnalyzer
					.getTransATModulePermutations(cluster);
		} else if (TypeIIPolyketideAnalyzer.isTypeIIPolyketideCluster(cluster)) {
			modulePermutations = TypeIIPolyketideGenerator
					.generateTypeIIPolyketideModulePermutation(cluster);
			cd.setNumOrfPermutations(1);
		} else if (cluster.families().contains(ClusterFamilies.RIBOSOMAL)) {
			// cleave leader sequence from precursor for a final propeptide
			System.out.println("[ClusterPermuter] Getting propeptide module permutations");
			for (Orf orf : cluster.orfs()) {
				/*
				 * TODO: In ribosomal PRISM, "domain" refers to an individual
				 * residue within a sequence (i.e. propeptide), in order to
				 * force existing data structures to work with ribosomal peptide
				 * prediction. The method called here in
				 * RibosomalClusterAnalyzer converts the sequence into the
				 * appropriate data structure. This is an awkward compromise.
				 * Eventually, this should be refactored to create a graph-based
				 * prediction system.
				 */

				// generate module permutations from propeptides
				List<Propeptide> propeptides = orf.propeptides();
				for (Propeptide p : propeptides) {
					modulePermutations.addAll(RibosomalClusterAnalyzer
							.generateRibosomalModulePermutations(
									p.getSequence(), orf));
				}
				cd.setNumPropeptides(modulePermutations.size());
			}
		} else {
			List<List<Orf>> orfPermutations = OrfPermuter.permuteOrfs(cluster);
			modulePermutations = ModulePermuter.generateModulePermutations(orfPermutations);
		}
		
		checkAcylAdenylateModules(modulePermutations, cluster);
		checkExtenderModules(modulePermutations);
		
		return modulePermutations;
	}

	/**
	 * If a module which cannot extend a scaffold is located in the middle of
	 * this permutation, inactivate the module.
	 * 
	 * @param permutations
	 *            module permutations to analyze
	 */
	public static void checkExtenderModules(List<List<Module>> permutations) {
		for (List<Module> permutation : permutations) {
			int i = 0;
			Iterator<Module> itr = permutation.iterator();
			while (itr.hasNext()) {
				Module next = itr.next();
				if (!next.canExtend() && i > 0)
					itr.remove();
				i++;
			}
		}
	}

	/**
	 * If a cluster has multiple non-extender acyl-adenylating domains,
	 * permutations which do not start with an acyl-adenylating module can be
	 * ignored.
	 * 
	 * @param permutations
	 *            module permutations to analyze
	 * @param cluster
	 *            parent cluster
	 */
	public static void checkAcylAdenylateModules(
			List<List<Module>> permutations, Cluster cluster) {
		int start = 0;
		for (Domain domain : cluster
				.domains(ThiotemplatedDomains.ACYL_ADENYLATING))
			if (domain.topSubstrate().smiles().indexOf("F") == -1)
				start++;

		if (start > 1) {
			Iterator<List<Module>> itr = permutations.iterator();
			while (itr.hasNext()) {
				List<Module> permutation = itr.next();
				if (permutation.size() > 0
						&& permutation.get(0).type() != ModuleTypes.ACYL_ADENYLATE)
					itr.remove();
			}
		}
	}
	
}