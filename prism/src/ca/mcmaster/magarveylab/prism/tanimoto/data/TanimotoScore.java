package ca.mcmaster.magarveylab.prism.tanimoto.data;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

import ca.mcmaster.magarveylab.prism.database.data.SmallMolecule;

/**
 * A Tanimoto coefficient comparing a cluster scaffold to a known natural
 * product.
 * 
 * @author skinnider
 *
 */
public class TanimotoScore implements Serializable {
	
	private static final long serialVersionUID = -7455099435535277149L;
	
	private SmallMolecule query;
	private SmallMolecule target;
	private Map<String, Float> scores = new HashMap<String, Float>();
	
	/**
	 * Instantiate a new Tanimoto score.
	 * @param query		cluster scaffold molecule
	 * @param target	database molecule being compared to the cluster scaffold
	 * @param type		name of the algorithm used to generate the two fingerprints 
	 * @param score		Tanimoto coefficient
	 */
	public TanimotoScore(SmallMolecule query, SmallMolecule target) {
		this.target = target;
		this.query = query;
	}
	
	/**
	 * Get the target (database) molecule this score represents a comparison with.
	 * @return	target molecule
	 */
	public SmallMolecule target() {
		return target;
	}
	
	/**
	 * Get query (predicted) molecule in this score.
	 * @return	query (scaffold) molecule
	 */
	public SmallMolecule query() {
		return query;
	}
	
	/**
	 * Get a Tanimoto coefficient associated with this query-target molecule pair.
	 * @param name	the name of the fingerprinting algorithm used to generate this score
	 * @return		the Tanimoto coefficient obtained from this fingerprinter 
	 */
	public float score(String name) {
		return scores.get(name);
	}
	
	/**
	 * Associate a new Tanimoto coefficient with this query-target molecule pair. 
	 * @param name		the name of the fingerprinting algorithm used to generate this score
	 * @param score		the Tanimoto coefficient obtained from these two fingerprints
	 */
	public void addScore(String name, Float score) {
		scores.put(name, score);
	}
	
}
