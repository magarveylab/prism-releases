package ca.mcmaster.magarveylab.prism.data.sugar;

import ca.mcmaster.magarveylab.enums.DeoxySugars;
import ca.mcmaster.magarveylab.enums.SugarFamilies;
import ca.mcmaster.magarveylab.enums.domains.DeoxySugarDomains;

/**
 * A deoxygenated sugar associated with a natural product and biosynthesized by
 * a cassette of specialized deoxysugar biosynthesis genes.
 * 
 * @author skinnider
 *
 */
public class DeoxySugar extends Sugar {

	private static final long serialVersionUID = -7307662856702946915L;
	
	private DeoxySugars type;

	public DeoxySugar(DeoxySugars type) {
		this.type = type;
		this.family = SugarFamilies.DEOXY;
	}
	
	public DeoxySugars type() {
		return type;
	}
	
	public String smiles() {
		return type.smiles();
	}
	
	public String name() {
		return type.fullName();
	}
	
	public DeoxySugarDomains[] genes() {
		return type.genes();
	}
	
	@Override
	public String toString() {
		return type.toString();
	}
	
}
