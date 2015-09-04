package ca.mcmaster.magarveylab.prism.data.reactions;

import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

/**
 * A package for one or more substrates used in combinatorial library
 * generation.
 * 
 * @author skinnider
 *
 */
public class SubstrateSet {

	private List<Module> modules = new ArrayList<Module>();

	/**
	 * Instantiate a new, empty substrate package. 
	 */
	public SubstrateSet() {
	}

	/**
	 * Instantiate a new deep copy of an existing substrate package.
	 * 
	 * @param substrate
	 *            substrate package to deep copy
	 */
	public SubstrateSet(SubstrateSet substrate) {
		for (Module module : substrate.getAllModules())
			this.modules.add(module);
	}

	/**
	 * Instantiate a new substrate package from a list of modules.
	 * 
	 * @param modules
	 *            list of substrate modules
	 */
	public SubstrateSet(List<Module> modules) {
		this.modules = modules;
	}

	/**
	 * Instantiate a new substrate package from a single module.
	 * 
	 * @param module
	 *            substrate module
	 */
	public SubstrateSet(Module module) {
		this.modules.add(module);
	}

	/**
	 * Instantiate a new substrate package from two modules.
	 * 
	 * @param first
	 *            first substrate module
	 * @param second
	 *            second substrate module
	 */
	public SubstrateSet(Module first, Module second) {
		this.modules.add(first);
		this.modules.add(second);
	}

	/**
	 * Get the number of substrate modules contained in this substrate package.
	 * 
	 * @return the number of substrate modules
	 */
	public int size() {
		return modules.size();
	}

	/**
	 * Add a new substrate module to this substrate package.
	 * 
	 * @param module
	 *            module to add
	 */
	public void add(Module module) {
		this.modules.add(module);
	}
	
	/**
	 * Add a list of new substrate modules to this substrate package.
	 * 
	 * @param modules
	 *            modules to add
	 */
	public void addAll(List<Module> modules) {
		this.modules.addAll(modules);
	}

	/**
	 * Get a substrate module from a substrate package.
	 * 
	 * @param n
	 *            index of the substrate module to get
	 * @return substrate module at index n
	 * @throws TailoringSubstrateException
	 */
	public Module get(int n) throws TailoringSubstrateException {
		Module module = modules.get(n);
		if (module == null)
			throw new TailoringSubstrateException("Error: module " + n
					+ " in tailoring plan is null!");
		return module;
	}

	/**
	 * Get all substrate modules from this substrate package.
	 * 
	 * @return all substrate modules
	 */
	public List<Module> getAllModules() {
		return modules;
	}

}
