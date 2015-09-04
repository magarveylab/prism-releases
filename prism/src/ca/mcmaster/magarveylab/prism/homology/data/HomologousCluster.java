package ca.mcmaster.magarveylab.prism.homology.data;

/**
 * Data package representing a potentially homologous gene cluster.
 * 
 * @author skinnider
 *
 */
public class HomologousCluster {

	private String name;
	private double identityScore = 0.0d;
	private double domainScore = 0.0d;
	private double coverage = 0.0d;

	/**
	 * Instantiate a new homologous gene cluster.
	 * 
	 * @param name
	 *            name of the subject (database) cluster
	 */
	public HomologousCluster(String name) {
		this.name = name;
	}

	/**
	 * Get the name of this homologous gene cluster.
	 * 
	 * @return name of the subject (database) cluster
	 */
	public String name() {
		return name;
	}

	/**
	 * Get the weighted average identity score as determined by BLAST
	 * preferential voting.
	 * 
	 * @return the identity score
	 */
	public double identityScore() {
		return identityScore;
	}

	/**
	 * Set the weighted average identity score as determined by BLAST
	 * preferential voting.
	 * 
	 * @param identityScore
	 *            the identity score
	 */
	public void setIdentityScore(double identityScore) {
		this.identityScore = identityScore;
	}

	/**
	 * Get the domain score as determined by the PRISM cluster dereplication
	 * algorithm.
	 * 
	 * @return the domain score
	 */
	public double domainScore() {
		return domainScore;
	}

	/**
	 * Get the coefficient of the domain score divided by the total number of
	 * bases covered.
	 * 
	 * @return the average domain score per aligned residue
	 */
	public double averageDomainScore() {
		return coverage == 0.0d ? 0.0d : domainScore / coverage;
	}

	/**
	 * Set the domain score as determined by the PRISM cluster dereplication
	 * algorithm
	 * 
	 * @param domainScore
	 *            the domain score
	 */
	public void setDomainScore(double domainScore) {
		this.domainScore = domainScore;
	}

	/**
	 * Get the coverage of this homologous cluster, defined as the total length
	 * of the domain score alignment divided by the size of the query cluster.
	 * 
	 * @return domain score alignment coverage
	 */
	public double coverage() {
		return coverage;
	}

	/**
	 * Set the coverage of this homologous cluster, defined as the total length
	 * of the domain score alignment divided by the size of the query cluster.
	 * 
	 * @param coverage
	 *            domain score alignment coverage
	 */
	public void setCoverage(double coverage) {
		this.coverage = coverage;
	}

}
