package ca.mcmaster.magarveylab.prism.combinatorialization;

import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.Orf;

/**
 * Generates module permutations.
 * 
 * @author skinnider
 *
 */
public class ModulePermuter {

	/**
	 * Generate a permutation of biosynthetic modules from each generated orf
	 * permutation.
	 * 
	 * @param orfPermutations
	 *            list of orf permutations
	 * @return permutation of biosynthetic modules corresponding to each orf
	 *         permutation
	 */
	public static List<List<Module>> generateModulePermutations(
			List<List<Orf>> orfPermutations) {
		List<List<Module>> modulePermutations = new ArrayList<List<Module>>();
		for (List<Orf> orfPermutation : orfPermutations) {
			List<Module> permutedModule = new ArrayList<Module>();
			for (Orf orf : orfPermutation)
				for (Module module : orf.modules())
					if (module.isActive())
						permutedModule.add(module);
			modulePermutations.add(permutedModule);
		}
		System.out.println("[ModulePermuter] Generated "
				+ modulePermutations.size() + " total module permutations");
		return modulePermutations;
	}

}
