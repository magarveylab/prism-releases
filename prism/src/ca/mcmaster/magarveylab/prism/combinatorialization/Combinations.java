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
	public static <T> void combinations(List<List<T>> lists, List<List<T>> results, int depth, List<T> current) {
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

	public static void combine(List<int[]> results, int[] arr, int k, int startIndex, int[] branch, int elementNumber) {
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

	public static List<int[]> sampleCombinations(int k, int limit) {
		List<int[]> combinations = new ArrayList<int[]>();
		Random random = new Random(0);
		int i = 0;
		while (i < limit) {
			int[] combination = new int[k];
			for (int j = 0; j < k; j++) {
				int r = random.nextInt(k) + 1;
				combination[j] = r;
			}
			if (!used(combination, combinations) || combinations.size() == 0) {
				combinations.add(combination);
				i++;
			}
		}
		return combinations;
	}
	
	public static boolean used(int[] combination, List<int[]> combinations) {
		boolean flag = true;
		for (int[] c : combinations) 
			for (int i = 0; i < c.length && i < combination.length; i++) 
				if (combination[i] != c[i])
					flag = false;
		return flag;
	}

}
