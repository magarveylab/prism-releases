package ca.mcmaster.magarveylab.prism.data;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.smiles.SmilesParser;


import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

/**
 * Represents a combinatorial library of possible natural products.
 * 
 * @author skinnider
 *
 */
public class Library implements Serializable {

	private static final long serialVersionUID = -8235025422333740682L;
	
	private int reactionCount = 0;
	private List<String> scaffolds = new ArrayList<String>();

	private List<Object> moleculeMasses = new ArrayList<Object>();
	private List<Double> fragmentation = new ArrayList<Double>();	
	
	
	public List<Object> getMoleculeMasses(){
		return moleculeMasses;
	}
	
	/**
	 * Returns a list of the top 50 most common fragment masses
	 * @return
	 */
	public List<Double> getFragments(){
		return fragmentation;
	}
	
	public void setFragments(List<Double> fragments){
		this.fragmentation = fragments;
	}
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
	 * Gets all the smiles in the library sorted according to molecular mass.
	 * 
	 * @return Map where the keys represent the mass and the values are a list of SMILES associated with that mass
	 * 
	 */
	public void sortByMass(){
		if (scaffolds.size() > 0){
			// Instantiates all the parsers and CDK stuff
			SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
			Map<String, Set<String>> stringOutput = new HashMap<String, Set<String>>();
			Set<String> imfstrings = new HashSet<String>();
			
			// Sorts all the scaffolds according to their molecular formula
			for(String s : scaffolds){
				IMolecule mol = null;
				try {
					mol = sp.parseSmiles(s);
				} catch (InvalidSmilesException e) {
					System.out.println("SMILES failed to parse");
					continue;
				}
				IMolecularFormula imf = MolecularFormulaManipulator.getMolecularFormula(mol);
				String imfstr = MolecularFormulaManipulator.getString(imf);
				if(imfstrings.contains(imfstr)){
					stringOutput.get(imfstr).add(s);
				}
				else{
					Set<String> toAdd = new HashSet<String>();
					toAdd.add(s);
					stringOutput.put(imfstr, toAdd);
					imfstrings.add(imfstr);
				}
			}
			// Converts molecular formula into mass and outputs SMILES
			for(String molform : stringOutput.keySet()){
				Map<String, Object> eachEntry = new HashMap<String,Object>();
				
				//Adding the mass section
				IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
				IMolecularFormula imfnew = MolecularFormulaManipulator.getMolecularFormula(molform, builder);
				double mass = MolecularFormulaManipulator.getNaturalExactMass(imfnew);
				eachEntry.put("monoisotopic_mass", mass);
				
				// Add the smiles section
				List<String> smiles = new ArrayList<String>();
				for(String smile : stringOutput.get(molform)){
					smiles.add(smile);
				}
				eachEntry.put("smiles", smiles);
				
				// Add it all to the main list
				moleculeMasses.add(eachEntry);
			}
		}
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
