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
	 * Generate all permutations of k elements of a set of n elements (nPk). Works by first generating all combinations nCk,
	 * then permuting each combination. 
	 * @param n		size of the set of elements
	 * @param k		number of elements per permutation
	 * @return		list of all permutations nPk
	 */
	public static List<int[]> permutations(int n, int k) {
		System.out.println("[Permutations] Generating all permutations " + n + "P" + k);
		
		if (k > n)
			System.out.println("[Permutations] Error: "
					+ "cannot permute when k is greater than n!");

		List<int[]> permutations = new ArrayList<int[]>();
		//if there are less than 6 elements it will permute all combos to get all permutations
		if (n < 6) {
			List<int[]> combinations = Combinations.combinations(n,k);
			for (int[] combination : combinations) {
				permutations.addAll(permutations(combination));
			}
		} else { //if there are 6 or more elements it will sample the combinations to get a sample of permutations
			permutations.addAll(samplePermutations(n, 500)); //returns 500 randomly generated permutations
		}	
		return permutations;
	}

	/**
	 * Generate a pseudorandom sample of permutations of an array of integers.
	 * 
	 * @param n
	 *            number of elements to put in the array
	 * @param limit
	 *            number of arrays to generate
	 * @return all permutations
	 */
	public static List<int[]> samplePermutations(int n, int limit) {
		List<int[]> permutations = new ArrayList<int[]>();
		int[] newpermute = new int[n];
		for (int y = 0; y < n; y++) {
			newpermute[y] = y;
		}
		Random rand = new Random(0);
		for (int v = 0; v < limit; v++) {
			boolean[] used = new boolean[n];
			int z = 0;
			int[] permuted = new int[n];
			while (z<n) {
				int rnd = rand.nextInt(n) +1;
				if (!used[rnd-1]) {
					permuted[z] = rnd;
					z++;
					used[rnd-1] = true;
				}
				
			}
			permutations.add(permuted);
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
