package ca.mcmaster.magarveylab.prism.cluster.reactions;

import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.paukov.combinatorics.Factory;
import org.paukov.combinatorics.Generator;
import org.paukov.combinatorics.ICombinatoricsVector;

import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Atoms;
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
	public static IAtom getBetaCarbon(Residue residue, IAtomContainer structure) throws TailoringSubstrateException {
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
	public static boolean isDehydrated(Residue residue, IAtomContainer structure) 
			throws TailoringSubstrateException {
		boolean flag = true;
		IAtom ketone = residue.ketone();
		IAtom ketoneOxygen = Atoms.getConnectedOxygen(ketone, structure);
		if (ketone == null) 
			throw new TailoringSubstrateException("Error: could not get ketone carbon!");
		if (ketoneOxygen == null)
			throw new TailoringSubstrateException("Error: could not get ketone oxygen!");
		for (IAtom atom : structure.atoms())
			if (atom.getSymbol().equals("O") && atom != ketoneOxygen) {
				System.out.println("Serine/threonine residue has not been dehydated");
				flag = false;
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
	 * @return the serine or threonine hydroxyl, or null if no hydroxyl is found
	 */
	public static IAtom getSerineOrThreonineHydroxyl(Residue residue, IAtomContainer structure) {
		IAtom hydroxyl = null;
		IAtom ketone = residue.ketone();
		IAtom ketoneOxygen = Atoms.getConnectedOxygen(ketone, structure);
		for (IAtom atom : structure.atoms())
			if (atom.getSymbol().equals("O") && atom != ketoneOxygen)
				hydroxyl = atom;
		return hydroxyl;
	}
	
	/**
	 * Get all subsets of a list. <br>
	 * In order to limit the amount of combinatorial information that is
	 * considered, when the input list is greater than 10, only all combinations
	 * where 0, 1, 2, or 3 of the input items are inactive will be considered.
	 * 
	 * @param input
	 *            list to get subsets of
	 * @return all subsets of the list (or a selection, for large lists--see
	 *         above)
	 */
	public static <T> List<List<T>> getSubsets(List<T> input) {
		List<List<T>> subsets = new ArrayList<List<T>>();
		ICombinatoricsVector<T> vector = Factory.createVector(input);
		if (input.size() <= 10) {
			// get all subsets of input list
			Generator<T> generator = Factory.createSubSetGenerator(vector);
			for (ICombinatoricsVector<T> subset : generator) {
				subsets.add(subset.getVector());
			}
		} else {
			// get all combinations with 0, 1, 2, or 3 domains inactive
			Generator<T> gen0 = Factory.createSimpleCombinationGenerator(vector, input.size());
			Generator<T> gen1 = Factory.createSimpleCombinationGenerator(vector, input.size() - 1);
			Generator<T> gen2 = Factory.createSimpleCombinationGenerator(vector, input.size() - 2);
			Generator<T> gen3 = Factory.createSimpleCombinationGenerator(vector, input.size() - 3);
			for (ICombinatoricsVector<T> subset : gen0) 
				subsets.add(subset.getVector());
			for (ICombinatoricsVector<T> subset : gen1) 
				subsets.add(subset.getVector());
			for (ICombinatoricsVector<T> subset : gen2) 
				subsets.add(subset.getVector());
			for (ICombinatoricsVector<T> subset : gen3) 
				subsets.add(subset.getVector());
		}
		return subsets;
	}
	
}
