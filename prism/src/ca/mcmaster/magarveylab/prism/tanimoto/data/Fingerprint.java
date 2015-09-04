package ca.mcmaster.magarveylab.prism.tanimoto.data;

import java.util.BitSet;

/**
 * A chemical fingerprint.
 * @author skinnider
 *
 */
public class Fingerprint {
	
	private String type;
	private BitSet bitset;

	/**
	 * Instantiate a new chemical fingerprint.
	 * @param name			the name of this molecule
	 * @param type			the fingerprinter used to generate the BitSet
	 * @param bitset		the BitSet that represents this molecule's chemical fingerprint
	 */
	public Fingerprint(String type, BitSet bitset) {
		this.type = type;
		this.bitset = bitset;
	}
		
	/**
	 * Get the fingerprinter used to generate this fingerprint.
	 * @return	type
	 */
	public String type() {
		return type;
	}
	
	/**
	 * Get this molecule's fingerprint.
	 * @return	this molecule's fingerprint
	 */
	public BitSet bitset() {
		return bitset;
	}

}
