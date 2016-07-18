package ca.mcmaster.magarveylab.prism.data.sugar;

import ca.mcmaster.magarveylab.enums.HexoseSugars;
import ca.mcmaster.magarveylab.enums.SugarFamilies;

/**
 * A hexose sugar, or microbial primary metabolite, associated with a natural
 * product.
 * 
 * @author skinnider
 *
 */
public class HexoseSugar extends Sugar {

	private static final long serialVersionUID = 7685888632881647197L;
	
	private HexoseSugars type;
	
	public HexoseSugar(HexoseSugars type) {
		this.type = type;
		this.family = SugarFamilies.HEXOSE;
	}
	
	public HexoseSugars type() {
		return type;
	}

	public String name() {
		return type.name();
	}
	
	public String smiles() {
		return type.smiles();
	}
	
	@Override
	public String toString() {
		return type.toString();
	}
	
}
