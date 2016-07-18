package ca.mcmaster.magarveylab.prism.cluster.reactions;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.exception.MathArithmeticException;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.paukov.combinatorics.Factory;
import org.paukov.combinatorics.Generator;
import org.paukov.combinatorics.ICombinatoricsVector;

import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Atoms;
import ca.mcmaster.magarveylab.prism.combinatorialization.Combinations;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

/**
 * Utilities class for ribosomal peptide-specific chemoinformatic operations or
 * reactions.
 * 
 * @author skinnider
 *
 */
public class RibosomalUtil {
	
	private static final int PLANS_PER_REACTION = 500;
	
	/**
	 * Get the beta carbon of an amino acid, used for serine/threonine
	 * dehydration.
	 * 
	 * @param residue
	 *            amino acid residue to analyze
	 * @param structure
	 *            the structure of the residue
	 * @return the beta carbon
	 * @throws TailoringSubstrateException
	 */
	public static IAtom getBetaCarbon(Residue residue, IAtomContainer structure)
			throws TailoringSubstrateException {
		IAtom betaCarbon = null;
		IAtom ketone = residue.ketone();
		if (ketone == null)
			throw new TailoringSubstrateException("Error: could not get ketone!");
		IAtom alphaCarbon = residue.alphaCarbon();
		if (alphaCarbon == null) 
			throw new TailoringSubstrateException("Error: could not get alpha carbon!");
		for (IAtom atom : structure.getConnectedAtomsList(alphaCarbon))
			if (atom.getSymbol().equals("C") && atom != ketone)
				betaCarbon = atom;
		if (betaCarbon == null)
			throw new TailoringSubstrateException("Error: could not get beta carbon!");
		return betaCarbon;
	}

	/**
	 * Check whether a serine or threonine residue has been dehydrated to form
	 * dehydroalanine or dehydroaminobutyric acid.
	 * 
	 * @param residue
	 *            amino acid residue to analyze
	 * @param structure
	 *            the structure of the residue
	 * @return true if the residue has only a single oxygen (the ketone)
	 * @throws TailoringSubstrateException
	 */
	public static boolean isDehydrated(Residue residue,
			IAtomContainer structure, IAtomContainer molecule)
			throws TailoringSubstrateException {
		List<IAtom> carboxyls = Atoms.getAllCarboxyls(molecule);
		boolean flag = true;
		IAtom ketone = residue.ketone();
		if (ketone == null) 
			throw new TailoringSubstrateException("Error: could not get ketone carbon!");
		for (IAtom atom : structure.atoms())
			if (atom.getSymbol().equals("O")
					&& !carboxyls.contains(atom)
					&& molecule.getBond(atom, ketone) == null) {
				boolean connectedToKetone = false;
				for (IAtom atom2 : structure.getConnectedAtomsList(atom))
					if (atom2 == ketone)
						connectedToKetone = true;
				if (!connectedToKetone)
					flag = false;
			}
		return flag;
	}

	/**
	 * Check whether an amino acid has an oxygen double-bonded to its ketone
	 * carbon; this is required for some reactions, such as 4+2 pyridine
	 * formation in thiopeptides.
	 * 
	 * @param residue
	 *            amino acid residue to analyze
	 * @param structure
	 *            the structure of the residue
	 * @return true if the residue has an oxygen double-bonded to its ketone
	 *         carbon
	 * @throws TailoringSubstrateException
	 *             if there is no ketone carbon
	 */
	public static boolean hasKetoneOxygen(Residue residue,
			IAtomContainer structure) throws TailoringSubstrateException {
		boolean flag = true;
		IAtom ketone = residue.ketone();
		if (ketone == null)
			throw new TailoringSubstrateException(
					"Error: could not get ketone carbon!");
		for (IAtom atom : structure.getConnectedAtomsList(ketone)) {
			IBond bond = structure.getBond(atom, ketone);
			if (atom.getSymbol().equals("O")
					&& bond.getOrder() == IBond.Order.DOUBLE)
				flag = true;
		}
		return flag;
	}
	
	/**
	 * Get the R-group hydroxyl of a serine or threonine residue.
	 * 
	 * @param residue
	 *            amino acid (serine/threonine) residue to analyze
	 * @param structure
	 *            the structure of the residue
	 * @param molecule
	 *            the parent molecule
	 * @return the serine or threonine hydroxyl, or null if no hydroxyl is found
	 */
	public static IAtom getSerineOrThreonineHydroxyl(Residue residue,
			IAtomContainer structure, IAtomContainer molecule) {
		IAtom hydroxyl = null;
		IAtom ketone = residue.ketone();
		List<IAtom> carboxyls = Atoms.getAllCarboxyls(molecule);
		IAtom ketoneOxygen = Atoms.getConnectedOxygen(ketone, structure);
		for (IAtom atom : structure.atoms()) {
			if (atom.getSymbol().equals("O") && atom != ketoneOxygen
					&& !carboxyls.contains(atom) && molecule.contains(atom)
					&& molecule.getBond(ketone, atom) == null) {
				hydroxyl = atom;			
			}
		}
		return hydroxyl;
	}
	
	/**
	 * Get all subsets of a list. <br>
	 * In order to limit the amount of combinatorial information that is
	 * considered, when the input list is greater than 4, only all combinations
	 * where 0, 1, or 2 of the input items are inactive, and combinations where
	 * 1 or 2 of the input items are active, will be considered. 
	 * 
	 * @param input
	 *            list to get subsets of
	 * @return all subsets of the list (or a selection, for large lists--see
	 *         above)
	 */
	public static <T> List<List<T>> getSubsets(List<T> input) {
		List<List<T>> subsets = new ArrayList<List<T>>();
		ICombinatoricsVector<T> vector = Factory.createVector(input);
		if (input.size() <= 4) {
			// get all subsets of input list
			Generator<T> generator = Factory.createSubSetGenerator(vector);
			for (ICombinatoricsVector<T> subset : generator) 
				if (subset.getVector().size() > 0)
					subsets.add(subset.getVector());
		} else {
			// if nCk > PLANS_PER_REACTION, then get a random subset from the set { nC1 -> nCn }
			List<int[]> combinations = Combinations.sampleCombinations(
					input.size(), 1, input.size(), PLANS_PER_REACTION);
			for (int[] combination : combinations) {
				List<T> subset = new ArrayList<T>();
				for (int c : combination)
					subset.add(input.get(c));
				subsets.add(subset);
			}
		}
		return subsets;
	}
	
	/**
	 * Get all subsets of a list of a fixed size, or, if more than 500 subsets
	 * are possible, a random subset.
	 * 
	 * @param input
	 *            list to get subsets of
	 * @param size
	 *            size of the subsets (combinations) to retrieve
	 * @return all subsets of the given size list (or a selection, for large
	 *         lists--see above)
	 */
	public static <T> List<List<T>> getSubsets(List<T> input, int size) {
		List<List<T>> subsets = new ArrayList<List<T>>();
		
		// if nCk > PLANS_PER_REACTION, then get a random subset
		long numerator = CombinatoricsUtils.factorial(input.size());
		long denominator = CombinatoricsUtils.factorial(input.size() - size)
				* CombinatoricsUtils.factorial(size);
		if (numerator / denominator > PLANS_PER_REACTION) {
			List<int[]> combinations = Combinations.sampleCombinations(
					input.size(), size, PLANS_PER_REACTION);
			for (int[] combination : combinations) {
				List<T> subset = new ArrayList<T>();
				for (int idx : combination)
					subset.add(input.get(idx));
				subsets.add(subset);
			}
		} else {
			ICombinatoricsVector<T> vector = Factory.createVector(input);
			Generator<T> gen = Factory.createSimpleCombinationGenerator(vector,
					size);
			for (ICombinatoricsVector<T> subset : gen)
				if (subset.getVector().size() > 0)
					subsets.add(subset.getVector());
		}
		
		return subsets;
	}
	
	/**
	 * Get all subsets of a list with a minimum size, or, if more than 500 such
	 * subsets are possible, a random subset.
	 * 
	 * @param input
	 *            list to get subsets of
	 * @param minSize
	 *            minimum size of the subsets to retrieve
	 * @return all subsets of the list with the given minimum size (or a
	 *         selection, for large lists--see above)
	 */
	public static <T> List<List<T>> getSubsetsWithMinSize(List<T> input, int minSize) {
		List<List<T>> subsets = new ArrayList<List<T>>();
		ICombinatoricsVector<T> vector = Factory.createVector(input);
		
		// if this query will produce >PLANS_PER_REACTION subsets, get a random selection
		int count = 0;
		long numerator = CombinatoricsUtils.factorial(input.size());
		for (int i = input.size(); i >= minSize; i--) {
			long denominator = CombinatoricsUtils.factorial(input.size() - i)
					* CombinatoricsUtils.factorial(i);
			count += numerator/denominator;
		}
		if (count > PLANS_PER_REACTION) {
			System.out.println("[RibosomalUtil] "
					+ "Getting a random selection of " + PLANS_PER_REACTION + " subsets of a "
					+ input.size() + "-item list with minimum size " + minSize
					+ " (" + count + " possible)");
			List<int[]> combinations = Combinations.sampleCombinations(
					input.size(), minSize, input.size(), PLANS_PER_REACTION);
			for (int[] combination : combinations) {
				List<T> subset = new ArrayList<T>();
				for (int c : combination)
					subset.add(input.get(c));
				subsets.add(subset);
			}
		} else {
			// get all subsets of input list
			Generator<T> generator = Factory.createSubSetGenerator(vector);
			for (ICombinatoricsVector<T> subset : generator) {
				if (subset.getVector().size() >= minSize)
					subsets.add(subset.getVector());
			}
		}
		
		return subsets;
	}
	
	/**
	 * Get all subsets of a list with a specified minimum size and a specified
	 * maximum size, or, if more than 500 such subsets are possible, a random
	 * subset.
	 * 
	 * @param input
	 *            list to get subsets of
	 * @param minSize
	 *            minimum size of the subsets to retrieve
	 * @param maxSize
	 *            maximum size of the subsets to retrieve
	 * @return all subsets of the list with the given size range (or a
	 *         selection, for large lists--see above)
	 */
	public static <T> List<List<T>> getSubsetsWithSizeRange(List<T> input,
			int minSize, int maxSize) {
		List<List<T>> subsets = new ArrayList<List<T>>();
		ICombinatoricsVector<T> vector = Factory.createVector(input);

		// if this query will produce >PLANS_PER_REACTION subsets, get a random selection
		int count = 0;
		try {
			long numerator = CombinatoricsUtils.factorial(input.size());
			for (int i = maxSize; i >= minSize; i--) {
				long denominator = CombinatoricsUtils.factorial(input.size() - i)
						* CombinatoricsUtils.factorial(i);
				count += numerator/denominator;
			}
		} catch (MathArithmeticException e) {
			count = PLANS_PER_REACTION + 1;
		}
		if (count > PLANS_PER_REACTION) {
			System.out.println("[RibosomalUtil] "
					+ "Getting a random selection of " + PLANS_PER_REACTION + " subsets of a "
					+ input.size() + "-item list with minimum size " + minSize
					+ " and maximum size " + maxSize + " (" + count
					+ " possible)");
			List<int[]> combinations = Combinations.sampleCombinations(
					input.size(), minSize, maxSize, PLANS_PER_REACTION);
			for (int[] combination : combinations) {
				List<T> subset = new ArrayList<T>();
				for (int c : combination)
					subset.add(input.get(c));
				subsets.add(subset);
			}
		} else {
			// get all subsets of input list
			Generator<T> generator = Factory.createSubSetGenerator(vector);
			for (ICombinatoricsVector<T> subset : generator) {
				if (subset.getVector().size() >= minSize)
					subsets.add(subset.getVector());
			}
		}
		
		return subsets;
	}
	
}
