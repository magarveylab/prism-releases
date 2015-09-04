package ca.mcmaster.magarveylab.prism.database.data;

import java.util.BitSet;

public class SmallMolecule {

	private int id;
	private String name;
	private String smiles;
	private BitSet fcfp6Fingerprint;
	private BitSet ecfp6Fingerprint;
	
	public String name() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}
	
	public String smiles() {
		return smiles;
	}
	
	public void setSmiles(String smiles) {
		this.smiles = smiles;
	}
	
	public int id() {
		return id;
	}
	
	public void setID(int id) {
		this.id = id;
	}
	
	public BitSet fcfp6Fingerprint() {
		return fcfp6Fingerprint;
	}
	
	public void setFcfp6Fingerprint(BitSet fcfp6Fingerprint) {
		this.fcfp6Fingerprint = fcfp6Fingerprint;
	}
	
	public BitSet ecfp6Fingerprint() {
		return ecfp6Fingerprint;
	}
	
	public void setEcfp6Fingerprint(BitSet ecfp6) {
		this.ecfp6Fingerprint = ecfp6;
	}
	
}
