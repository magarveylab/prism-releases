package ca.mcmaster.magarveylab.prism.combinatorialization;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.paukov.combinatorics.Factory;
import org.paukov.combinatorics.Generator;
import org.paukov.combinatorics.ICombinatoricsVector;

import ca.mcmaster.magarveylab.enums.Frames;
import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.prism.cluster.analysis.ClusterAnalyzer;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.CombinatorialData;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.Orf;

/**
 * Generates permutations of biosynthetic orfs.
 * 
 * @author skinnider
 *
 */
public class OrfPermuter {
	
	/**
	 * Generate all biosynthetically conceivable permutations of the open
	 * reading frames in this cluster, to account for non-colinearity.
	 * 
	 * @return all orf permutations
	 * @throws ClassNotFoundException
	 * @throws IOException
	 */
	public static List<List<Orf>> permuteOrfs(Cluster cluster) throws IOException {
		Orf start = ClusterAnalyzer.start(cluster.orfs());
		Orf end = ClusterAnalyzer.end(cluster.orfs());
	
		// get trans A domains
		List<Module> transAdenylation = cluster.activeModules(ModuleTypes.TRANS_ADENYLATION);
		List<Module> transAdenylationInsertion = cluster.modules(ModuleTypes.TRANS_ADENYLATION_INSERTION);
		
		// permute orfs
		List<List<Orf>> permutations = new ArrayList<List<Orf>>();
		List<Orf> orfs = cluster.moduleOrfs();
		System.out.println(cluster.frame().toString());
		if (transAdenylation.size() > 0 && transAdenylationInsertion.size() > 0) {
			System.out.println("[OrfPermuter] Getting trans-A permutations for cluster " + cluster.index());
			permutations = permuteOrfsWithTransAdenylationModules(start, end, cluster);
		} else {
			if (cluster.frame() == Frames.POSITIVE || cluster.frame() == Frames.NEGATIVE) {
				System.out.println("[OrfPermuter] Getting colinear permutation for cluster " + cluster.index());
				permutations = getColinearPermutation(start, end, orfs, cluster);
			} else {
				System.out.println("[OrfPermuter] Getting all permutations for cluster " + cluster.index());
				permutations = permuteOrfs(start, end, orfs, cluster);
			}
		}
		
		CombinatorialData cd = cluster.combinatorialData();
		cd.setNumOrfPermutations(permutations.size());
		
		System.out.println("[OrfPermuter] Generated " + permutations.size() + " permutations");
		return permutations;
	}
	
	/**
	 * Get permutations of orfs assuming colinearity.
	 * 
	 * @param start
	 *            start orf
	 * @param end
	 *            end orfs
	 * @param orfs
	 *            orfs to permute
	 * @param cluster
	 *            the parent cluster 
	 * @return all biosynthetically conceivable orf permutations
	 */
	public static List<List<Orf>> getColinearPermutation(Orf start, Orf end,
			List<Orf> orfs, Cluster cluster) {
		// get orfs in right order
		Frames frame = cluster.frame();
		if (frame == Frames.NEGATIVE)
			Collections.reverse(orfs);
		
		// if orfs don't start with start/don't end with end: permute
		if ((start != null && orfs.indexOf(start) != 0)
				|| (end != null && orfs.indexOf(end) != orfs.size() - 1)) {
			System.out.println("[OrfPermuter] Error: orfs are not colinear. Permuting orfs:");
			return permuteOrfs(start, end, orfs, cluster);
		}
		
		// else, generate a single permutation
		List<List<Orf>> permutations = new ArrayList<List<Orf>>();
		permutations.add(orfs);
		
		return permutations;
	}
	
	/**
	 * Get permutations of orfs from a cluster which contains at least one
	 * trans-adenylation module and one trans-adenylation insertion module.
	 * 
	 * @param start
	 *            start orf
	 * @param end
	 *            end orf
	 * @param cluster
	 *            the parent cluster 
	 * @return all biosnythetically conceivable orf permutations
	 */
	public static List<List<Orf>> permuteOrfsWithTransAdenylationModules(Orf start, Orf end, Cluster cluster) {
		List<List<Orf>> orfPermutations = new ArrayList<List<Orf>>();

		List<Module> transA = cluster.activeModules(ModuleTypes.TRANS_ADENYLATION);
		List<Module> transAInsertion = cluster.modules(ModuleTypes.TRANS_ADENYLATION_INSERTION);
		List<Orf> orfs = cluster.moduleOrfs();
		
		// get all permutations of trans A
		List<List<Module>> transAPermutations = new ArrayList<List<Module>>();
		ICombinatoricsVector<Module> transAVector = Factory.createVector(transA);
		Generator<Module> transAGenerator = Factory.createPermutationGenerator(transAVector);
		for (ICombinatoricsVector<Module> vector : transAGenerator) {
			List<Module> permutation = vector.getVector();
			transAPermutations.add(permutation);
		}
		
		for (List<Module> transAPermutation : transAPermutations) {
			// get combinations of two lists 
			List<List<Module>> toPermute = new ArrayList<List<Module>>();
			toPermute.add(transAInsertion);
			toPermute.add(transAPermutation);
			List<List<Module>> results = Combinations.getCombinations(toPermute);
			System.out.println("[OrfPermuter] Generated " + results.size() + " combinations of two lists");
		
			// create orfs copy
			List<Orf> orfsToPermute = new ArrayList<Orf>();
			for (Orf orf : orfs) {
				Orf copy = new Orf(orf);
				if (orf == end)
					end = copy;
				if (orf == start)
					start = copy;
				orfsToPermute.add(copy);
			}
				
			// insert trans A permutations into insertion modules 
			for (List<Module> result : results) {
				Module insertion = result.get(0);
				Module trans = result.get(1);

				for (Orf orf : orfs) {
					if (orf.contains(trans)) {
						Iterator<Orf> itr = orfsToPermute.iterator();
						while (itr.hasNext()) {
							Orf copy = itr.next();
							if (copy.name().equals(orf.name()))
								itr.remove();
						}
						continue;
					}
					for (Module module : orf.modules())
						if (module == insertion) {
							System.out.println("[OrfPermuter] Identified insertion module in " + orf.name());
							int idx = orf.modules().indexOf(insertion);
							for (Orf copy : orfsToPermute)
								if (copy.name().equals(orf.name())) 
									copy.modules().add(idx, new Module(trans));
						}
				}
			}
			
			List<List<Orf>> permuted;
			if (cluster.frame() == Frames.POSITIVE || cluster.frame() == Frames.NEGATIVE) {
				System.out.println("[OrfPermuter] Getting colinear permutation for cluster " + cluster.index());
				permuted = getColinearPermutation(start, end, orfsToPermute, cluster);
			} else {
				System.out.println("[OrfPermuter] FLAG: Detected non-colinear cluster ");
				System.out.println("[OrfPermuter] Getting all permutations for cluster " + cluster.index());
				permuted = permuteOrfs(start, end, orfsToPermute, cluster);
			}
			orfPermutations.addAll(permuted);
		}
		
		return orfPermutations;
	}

	/**
	 * Generate all permutations of a list of orfs, taking into account start
	 * and end orfs.
	 * 
	 * @param start
	 *            the start orf
	 * @param end
	 *            the end orf
	 * @param orfs
	 *            the list of orfs
	 * @param cluster
	 *            the parent cluster
	 * @return all permutations of the orfs
	 */
	public static List<List<Orf>> permuteOrfs(Orf start, Orf end, List<Orf> orfs, Cluster cluster) {
		List<List<Orf>> newPermutations = new ArrayList<List<Orf>>();
		System.out.println("[OrfPermuter] Permuting " + orfs.size() + " orfs");
		
		// remove start & end from permutable array
		List<Orf> permutable = new ArrayList<Orf>();
		for (Orf orf : orfs) {
			if (orf != start && orf != end)
				permutable.add(orf);
		}
	
		// if there are only start/end orfs, don't permute anything
		if (permutable.size() == 0) {
			List<Orf> permutation = new ArrayList<Orf>();
			if (start != null)
				permutation.add(start);
			if (end != null)
				permutation.add(end);
			newPermutations.add(permutation);
			return newPermutations;
		}
		
		List<List<Orf>> permutations = permuteOrfs(permutable);
		System.out.println("[OrfPermuter] Generated " + permutations.size()
				+ " permutations of " + permutable.size() + " orfs");
		
		// prepend start/append end to each orf order
		for (List<Orf> permutation : permutations) {
			if (start != null)
				permutation.add(0, start);
			if (end != null)
				permutation.add(end);
			newPermutations.add(permutation);
		}
		return newPermutations;
	}
	
	/**
	 * Generate all permutations of a list of orfs.
	 * 
	 * @param orfs
	 *            orfs to permute
	 * @return all permutations, or 500 random permutations if there are too
	 *         many to search the entire space 
	 */
	public static List<List<Orf>> permuteOrfs(List<Orf> orfs) {
		List<List<Orf>> permutedOrfs = new ArrayList<List<Orf>>();
		
		int size = orfs.size();
		System.out.println("[OrfPermuter] After start/end analysis, permuting " + size + " orfs");
		
		List<int[]> permutations = Permutations.permutations(size, size, 500);
		for (int[] permutation : permutations) {
			List<Orf> orfPermutation = new ArrayList<Orf>();
			for (int i = 0; i < permutation.length; i++)
				orfPermutation.add(i, orfs.get(permutation[i] - 1));
			permutedOrfs.add(orfPermutation);
		}
			
		return permutedOrfs;
	}

}
