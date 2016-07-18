package ca.mcmaster.magarveylab.prism.combinatorialization;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Low-level class to generate integer combinations.
 * @author skinnider
 *
 */
public class Combinations {


	/**
	 * Get all combinations of multiple lists, or when the number of combinations is computationally expensive (>1,000,000), 
	 * a random subset.
	 * 
	 * @param lists
	 *            list of lists to combine
	 * @param results
	 *            object which holds resulting combinations
	 * @param depth
	 *            current list depth
	 * @param current
	 *            current combination
	 */
	public static <T> List<List<T>> getCombinations(List<List<T>> lists) {
		List<List<T>> output = new ArrayList<List<T>>();
		
		// how many combinations will the complete set produce?
		int count = 1;
		for (List<T> list : lists)
			count *= list.size();
		
		// if >1 million, sample
		if (count <= 1_000_000) {
			combinations(lists, output, 0, new ArrayList<T>());
		} else {
			output = sampleListCombinations(lists, 1_000_000);
		}

		return output;
	}
	
	public static <T> List<List<T>> sampleListCombinations(List<List<T>> lists, int limit) {
		Random random = new Random(0);
		int i = 0, used = 0;
		List<int[]> combinations = new ArrayList<int[]>();
		while (i < limit) {
			if (used > 1_000)
				break;
			
			int[] combination = new int[lists.size()];
			for (int j = 0; j < lists.size(); j++) {
				List<T> list = lists.get(j);
				int r = random.nextInt(list.size());
				combination[j] = r;
			}
			if (!usedCombination(combination, combinations)
					|| combinations.size() == 0) {
				combinations.add(combination);
				i++;
			} else {
				used++;
			}
		}
		
		List<List<T>> output = new ArrayList<List<T>>();
		for (int[] combination : combinations) {
			List<T> list = new ArrayList<T>();
			for (int a = 0; a < combination.length; a++) {
				int c = combination[a];
				list.add(lists.get(a).get(c));
			}		
			output.add(list);
		}
		
		return output;
	}
	
	/**
	 * Recursive function to generate all combinations of multiple lists. For
	 * two lists { 1, 2 } and { A, B, C }, this function will generate { 1A, 1B,
	 * 1C, 2A, 2B, 2C }.
	 * 
	 * @param lists
	 *            list of lists to combine
	 * @param results
	 *            object which holds resulting combinations
	 * @param depth
	 *            current list depth
	 * @param current
	 *            current combination
	 */
	public static <T> void combinations(List<List<T>> lists,
			List<List<T>> results, int depth, List<T> current) {
		if (depth == lists.size()) {
			List<T> result = new ArrayList<T>();
			for (int i = 0; i < current.size(); i++)
				result.add(current.get(i));
			results.add(result);
			current.clear();
			return;
		}
		
		if (lists.get(depth).size() > 0) {
			for (int i = 0; i < lists.get(depth).size(); i++) {
				if (depth > 0 && results.size() > 0 && current.size() == 0) {
					int lastIdx = results.size() - 1;
					List<T> last = results.get(lastIdx);
					for (int j = 0; j < depth; j++) {
						current.add(last.get(j));
					}
				}
				current.add(lists.get(depth).get(i));
				combinations(lists, results, depth + 1, current);
			}
		} else {
			current.add(null);
			combinations(lists, results, depth + 1, current);
		}
	}

	/**
	 * Generate all combinations of k elements of a set of n elements.
	 * @param n		size of the set of elements
	 * @param k		number of elements to consider
	 * @return		all possible combinations nCk
	 */
	public static List<int[]> combinations(int n, int k) {
		List<int[]> combinations = new ArrayList<int[]>();
		
		if (n == 0) {
			int[] combination = new int[] { 1 };
			combinations.add(combination);
			return combinations;
		}
		
		// instantiate array of integers from 1 to n 
		int[] original = new int[n];
		for (int i = 0; i < n; i++)
			original[i] = i + 1;
		
		int[] branch = new int[k];
		Combinations.combine(combinations, original, k, 0, branch, 0);
		
		return combinations;
	}

	public static void combine(List<int[]> results, int[] arr, int k,
			int startIndex, int[] branch, int elementNumber) {
		if (elementNumber == k) {
			results.add(Arrays.copyOf(branch, branch.length));
			return;
		} 
		for (int i = startIndex; i < arr.length; ++i) {
			branch[elementNumber++] = arr[i];
			combine(results, arr, k, ++startIndex, branch, elementNumber);
			--elementNumber;
		}
	}

	/**
	 * Sample randomly from a set of permutations nPn.
	 * 
	 * @param n
	 *            the number of items to permute
	 * @param limit
	 *            the number of permutations to sample
	 * @return a random set of permutations from the set nPn
	 */
	public static List<int[]> samplePermutations(int n, int limit) {
		List<int[]> permutations = new ArrayList<int[]>();
		Random random = new Random(0);
		int i = 0;
		while (i < limit) {
			int[] permutation = new int[n];
			Arrays.fill(permutation, -1);
			for (int j = 0; j < n; j++) {
				int r = -1;
				while (r == -1 || contains(permutation, r))
					r = random.nextInt(n);
				permutation[j] = r;
			}
			if (!usedPermutation(permutation, permutations)
					|| permutations.size() == 0) {
				permutations.add(permutation);
				i++;
			}
		}
		return permutations;
	}

	/**
	 * Sample randomly from a set of combinations nCk.
	 * 
	 * @param n
	 *            the number of items that can be combined
	 * @param k
	 *            the size of each combination
	 * @param limit
	 *            the number of combinations to sample
	 * @return a random set of combinations from the set nCk
	 */
	public static List<int[]> sampleCombinations(int n, int k, int limit) {
		if (k > n)
			throw new ArrayIndexOutOfBoundsException(
					"Couldn't sample combinations: k (" + k + ") > n (" + n
							+ ")");
		List<int[]> combinations = new ArrayList<int[]>();
		Random random = new Random(0);
		int i = 0, used = 0;
		while (i < limit) {
			// if 1,000 used combinations have been generated, stop trying 
			if (used > 1_000)
				break; 
			int[] combination = new int[k];
			Arrays.fill(combination, -1);
			for (int j = 0; j < k; j++) {
				int r = -1;
				while (r == -1 || contains(combination, r))
					r = random.nextInt(n);
				combination[j] = r;
			}
			if (!usedCombination(combination, combinations)
					|| combinations.size() == 0) {
				combinations.add(combination);
				i++;
			} else {
				used++; 
			}
		}
		return combinations;
	}
	
	/**
	 * Sample randomly from sets of combinations nC{k}, where k is an integer
	 * whose bounds are specified by kMin and kMax.
	 * 
	 * @param n
	 *            the number of items that can be combined
	 * @param kMin
	 *            the minimum size of each combination
	 * @param kMax
	 *            the maximum size of each combination
	 * @param limit
	 *            the number of combinations to sample
	 * @return a random set of combinations from the sets nCkMin through nCkMax
	 */
	public static List<int[]> sampleCombinations(int n, int kMin, int kMax, int limit) {
		if (kMax > n)
			throw new ArrayIndexOutOfBoundsException(
					"Couldn't sample combinations: upper limit " + kMax
							+ " is greater than size " + n);
		List<int[]> combinations = new ArrayList<int[]>();
		Random random = new Random(0);
		int i = 0, used = 0;
		while (i < limit) {
			// if 1,000 used combinations have been generated, stop trying 
			if (used > 1_000)
				break; 
			int k = random.nextInt(kMax - kMin + 1) + kMin;
			int[] combination = new int[k];
			Arrays.fill(combination, -1);
			for (int j = 0; j < k; j++) {
				int r = -1;
				while (r == -1 || contains(combination, r))
					r = random.nextInt(n);
				combination[j] = r;
			}
			if (!usedCombination(combination, combinations)
					|| combinations.size() == 0) {
				combinations.add(combination);
				i++;
			} else {
				used++;
			}
		}
		return combinations;
	}
	
	public static boolean usedPermutation(int[] permutation, List<int[]> permutations) {
		boolean flag = false;
		for (int[] p : permutations) {
			boolean same = true;
			if (permutation.length != p.length)
				continue;
			for (int i = 0; i < p.length && i < permutation.length; i++)
				if (permutation[i] != p[i])
					same = false;
			if (same)
				flag = true;
		}
		return flag;
	}
	
	public static boolean usedCombination(int[] combination, List<int[]> combinations) {
		for (int[] c : combinations) {
			boolean same = true;
			if (combination.length != c.length)
				continue;
			for (int i : combination) {
				if (!contains(c, i))
					same = false;
			}
			if (same) 
				return true;
		}
		return false;
	}

	public static boolean contains(int[] array, int i) {
		boolean contains = false;
		for (int j : array)
			if (i == j)
				contains = true;
		return contains;
	}
	
}
