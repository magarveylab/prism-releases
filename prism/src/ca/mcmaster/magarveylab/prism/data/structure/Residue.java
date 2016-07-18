package ca.mcmaster.magarveylab.prism.data.structure;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.prism.data.Module;

/**
 * A single residue or monomer of a growing natural product scaffold.
 * 
 * @author skinnider
 *
 */
public class Residue {

	private IAtom ketone;
	private IAtom alphaCarbon;
	private IAtom nitrogen;
	private IAtomContainer structure;
	private Module module;

	/**
	 * Instantiate a new residue from a module.
	 * 
	 * @param module
	 *            module corresponding to this residue
	 */
	public Residue(Module module) {
		this.module = module;
	}

	/**
	 * Get the ketone carbon in this residue.
	 * 
	 * @return the ketone carbon
	 */
	public IAtom ketone() {
		return ketone;
	}

	/**
	 * Get the alpha carbon in this residue.
	 * 
	 * @return the alpha carbon
	 */
	public IAtom alphaCarbon() {
		return alphaCarbon;
	}

	/**
	 * Get the nitrogen in this residue.
	 * 
	 * @return the nitrogen
	 */
	public IAtom nitrogen() {
		return nitrogen;
	}

	/**
	 * Get the entire structure of this residue.
	 * 
	 * @return the residue's structure
	 */
	public IAtomContainer structure() {
		return structure;
	}

	/**
	 * Get the module associated with this residue.
	 * 
	 * @return the module corresponding to this residue
	 */
	public Module module() {
		return module;
	}

	/**
	 * Get the type of the module associated with this residue.
	 * 
	 * @return the residue type
	 */
	public ModuleTypes type() {
		return module.type();
	}

	/**
	 * Set the nitrogen atom within this residue.
	 * 
	 * @param nitrogen
	 *            the nitrogen atom
	 */
	public void setNitrogen(IAtom nitrogen) {
		this.nitrogen = nitrogen;
	}

	/**
	 * Set the ketone atom within this residue.
	 * 
	 * @param ketone
	 *            the ketone carbon
	 */
	public void setKetone(IAtom ketone) {
		this.ketone = ketone;
	}

	/**
	 * Set the alpha carbon within this residue.
	 * 
	 * @param alphaCarbon
	 *            the alpha carbon
	 */
	public void setAlphaCarbon(IAtom alphaCarbon) {
		this.alphaCarbon = alphaCarbon;
	}

	/**
	 * Set the structure of this residue.
	 * 
	 * @param structure
	 *            the residue structure
	 */
	public void setStructure(IAtomContainer structure) {
		this.structure = structure;
	}

}
