package ca.mcmaster.magarveylab.prism.data;

import java.util.ArrayList;
import java.util.List;

/**
 * Represents a combinatorial library of possible natural products.
 * 
 * @author skinnider
 *
 */
public class Library {

	private int reactionCount = 0;
	private List<String> scaffolds = new ArrayList<String>();

	/**
	 * Add a new scaffold to this library.
	 * 
	 * @param scaffold
	 *            scaffold to add
	 */
	public void add(String scaffold) {
		scaffolds.add(scaffold);
	}

	/**
	 * Get all scaffolds within this library.
	 * 
	 * @return all scaffolds
	 */
	public List<String> scaffolds() {
		return scaffolds;
	}

	/**
	 * Get the size of this library--i.e., the total number of generated
	 * scaffolds.
	 * 
	 * @return library size
	 */
	public int size() {
		return scaffolds.size();
	}

	/**
	 * Wipe the library of all generated SMILES.
	 */
	public void clear() {
		this.scaffolds = new ArrayList<String>();
	}

	/**
	 * Test whether this library contains a given scaffold.
	 * 
	 * @param scaffold
	 *            scaffold to check
	 * @return true if this scaffold is already contained within this library
	 */
	public boolean contains(String scaffold) {
		boolean flag = false;
		for (String s : scaffolds)
			if (s.equals(scaffold))
				flag = true;
		return flag;
	}

	/**
	 * Get the total number of tailoring reactions executed in this scaffold
	 * library.
	 * 
	 * @return the total number of scaffold tailoring reactions executed on each
	 *         scaffold
	 */
	public int reactionCount() {
		return reactionCount;
	}

	/**
	 * Set the the total number of tailoring reactions executed in this scaffold
	 * library.
	 * 
	 * @param the
	 *            total number of scaffold tailoring reactions executed on each
	 *            scaffold
	 */
	public void setReactionCount(int reactions) {
		this.reactionCount = reactions;
	}

}
