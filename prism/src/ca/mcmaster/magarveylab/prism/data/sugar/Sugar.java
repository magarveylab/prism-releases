package ca.mcmaster.magarveylab.prism.data.sugar;

import java.io.Serializable;

import ca.mcmaster.magarveylab.enums.SugarFamilies;
import ca.mcmaster.magarveylab.enums.interfaces.Structure;

/**
 * A generic sugar associated with natural product biosynthesis. 
 * @author skinnider
 *
 */
public abstract class Sugar implements Serializable {
	
	private static final long serialVersionUID = -1800418584985457039L;
	
	protected SugarFamilies family;
	
	public SugarFamilies family() {
		return family;
	}
	
	public boolean isDeoxygenated() {
		return family == SugarFamilies.DEOXY;
	}
	
	public boolean isHexose() {
		return family == SugarFamilies.HEXOSE;
	}
	
	public abstract String smiles();
	
	public abstract String name();

	public abstract Structure type();
	
}
