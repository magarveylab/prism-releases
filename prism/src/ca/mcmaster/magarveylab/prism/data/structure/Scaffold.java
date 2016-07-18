package ca.mcmaster.magarveylab.prism.data.structure;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;

/**
 * A natural product scaffold.
 * @author skinnider
 *
 */
public class Scaffold {
	
	private int reactionCount = 0;
	private IAtomContainer molecule = new AtomContainer();
	private Map<Module,Residue> residues  = new LinkedHashMap<Module,Residue>();
	
	/**
	 * Instantiate a new scaffold.
	 */
	public Scaffold() { 
	}
	
	/**
	 * Get the number of reactions that were successfully executed in the construction of this scaffold.
	 * @return	the total number of successfully executed reactions
	 */
	public int reactionCount() {
		return reactionCount;
	}
	
	/**
	 * Increment the number of reactions that have been successfully executed in the construction of this scaffold. 
	 */
	public void incrementReactionCount() {
		reactionCount++;
	}
	
	/**
	 * Get the molecule associated with this scaffold.
	 * @return	the molecule associated with this scaffold
	 */
	public IAtomContainer molecule() {
		return molecule;
	}
	
	/**
	 * Set the molecule associated with this scaffold.
	 * @param container	the molecule associated with this scaffold
	 */
	public void setMolecule(IAtomContainer container) {
		this.molecule = container;
	}
	
	/**
	 * Get the map of modules to residues associated with this scaffold.
	 * @return	the map of modules to residues
	 */
	public Map<Module,Residue> residues() {
		return residues;
	}
	
	/**
	 * Set the map of modules to residues associated with this scaffold.
	 * @param residues	the map of modules to residues 
	 */
	public void setResidues(Map<Module,Residue> residues) {
		this.residues = residues;
	}
	
	/**
	 * Get the residue corresponding to a given module in this scaffold.
	 * @param module	the module to check	
	 * @return			the corresponding residue
	 * @throws NoResidueException 
	 */
	public Residue residue(Module module) throws NoResidueException {
		try {
			Residue residue = residues.get(module);
			return residue;
		} catch (ArrayIndexOutOfBoundsException e) {
			String name = (module.domains().size() > 0) ? module.domains().get(0).name() : "empty module";
			throw new NoResidueException("Error: could not get scaffold residue for " + module.domains().size() + 
					"-domain module with first domain: " + name, e);
		}
	}
	
	/**
	 * Get the last residue in this scaffold.
	 * @return	the last residue
	 */
	public Residue getLastResidue() {
		List<Residue> values = getAllResidues();
		int last = residues.size() - 1;
		return values.get(last);
	}
	
	/**
	 * Get all residues associated with this scaffold.
	 * @return	all residues
	 */
	public List<Residue> getAllResidues() {
		return new ArrayList<Residue>(residues.values());
	}
	
	/**
	 * Add a new module-residue pair to this scaffold.	
	 * @param module	the module to add
	 * @param residue	the residue to add
	 */
	public void addResidue(Module module, Residue residue) {
		residues.put(module, residue);
	}
	
	/**
	 * Get the index of a module within this scaffold. 
	 * @param module	module to search 
	 * @return			the module's index 
	 */
	public int indexOf(Module module) {
		Object[] modules = residues.keySet().toArray();
		return (Arrays.asList(modules).indexOf(module));
	}
	
}
