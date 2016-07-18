package ca.mcmaster.magarveylab.prism.combinatorialization;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Low-level class to generate integer permutations.
 * @author skinnider, prees
 *
 */
public class Permutations {

	/**
	 * Generate all permutations of k elements of a set of n elements (nPk).
	 * Works by first generating all combinations nCk, then permuting each
	 * combination.
	 * 
	 * @param n
	 *            size of the set of elements
	 * @param k
	 *            number of elements per permutation
	 * @return list of all permutations nPk
	 */
	public static List<int[]> permutations(int n, int k, int limit) {
		System.out.println("[Permutations] Generating all permutations " + n
				+ "P" + k);

		if (k > n)
			System.out.println("[Permutations] Error: "
					+ "cannot permute when k is greater than n!");

		List<int[]> permutations = new ArrayList<int[]>();
		//if there are less than 6 elements it will permute all combos to get all permutations
		if (n < 6) {
			List<int[]> combinations = Combinations.combinations(n, k);
			for (int[] combination : combinations) 
				permutations.addAll(permutations(combination));
		} else { // if there are 6 or more elements it will sample the
					// combinations to get a sample of permutations
			permutations.addAll(samplePermutations(n, k, limit));
		}

		System.out.println("[Permutations] Generated " + permutations.size()
				+ " permutations " + n + "P" + k);
		return permutations;
	}
	
	public static List<int[]> samplePermutations(int n, int k, int limit) {
		if (k > n)
			throw new ArrayIndexOutOfBoundsException("Couldn't sample "
					+ "combinations: k (" + k + ") > n (" + n
							+ ")");
		List<int[]> permutations = new ArrayList<int[]>();
		Random random = new Random(0);
		int i = 0, used = 0;
		while (i < limit) {
			// if 1,000 used permutations have been generated, stop trying 
			if (used > 10_000)
				break; 
			int[] permutation = new int[k];
			Arrays.fill(permutation, -1);
			for (int j = 0; j < k; j++) {
				int r = -1;
				while (r == -1 || Combinations.contains(permutation, r))
					r = random.nextInt(n) + 1;
				permutation[j] = r;
			}
			if (!Combinations.usedPermutation(permutation, permutations)
					|| permutations.size() == 0) {
				permutations.add(permutation);
				i++;
			} else {
				used++;
			}
		}
		return permutations;
	}

	
	/**
	 * Generate all permutations of an array of integers.
	 * 
	 * @param array
	 *            array to permute
	 * @return all permutations
	 */
	public static List<int[]> permutations(int[] array) {
		List<int[]> permutations = new ArrayList<int[]>();
		permute(permutations, array, 0);
		return permutations;
	}

	public static void permute(List<int[]> results, int[] array, int k) {
		for (int i = k; i < array.length; i++) {
			swap(array, i, k);
			permute(results, array, k+1);
			swap(array, k, i);
		}
		if (k == array.length - 1) {
			results.add(Arrays.copyOf(array, array.length));
		}
	}

	public static void swap(int[] array, int k, int i) {
		int temp = array[i];
		array[i] = array[k];
		array[k] = temp;
	}

}
