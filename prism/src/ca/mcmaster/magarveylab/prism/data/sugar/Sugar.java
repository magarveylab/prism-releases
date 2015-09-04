package ca.mcmaster.magarveylab.prism.data.sugar;

import ca.mcmaster.magarveylab.enums.Structure;
import ca.mcmaster.magarveylab.enums.SugarFamilies;

/**
 * A generic sugar associated with natural product biosynthesis. 
 * @author skinnider
 *
 */
public abstract class Sugar {
	
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
