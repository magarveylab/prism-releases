package ca.mcmaster.magarveylab.prism.data;

import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.enums.DomainFamilies;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.prism.blast.BlastSearchResult;
import ca.mcmaster.magarveylab.prism.data.sugar.Sugar;
import ca.mcmaster.magarveylab.prism.genome.data.HmmSearchResultAnnotation;

/**
 * A single domain within an orf. 
 * @author skinnider
 *
 */
public class Domain {
	
	private final int start;
	private final int end;
	private final double score;
	private String name;
	private String sequence;
	
	private DomainFamilies family;
	private DomainType type;

	private List<BlastSearchResult> blastResults = new ArrayList<BlastSearchResult>();
	private List<Substrate> substrates = new ArrayList<Substrate>();

	private Sugar sugar; // if GTr
	private BlastSearchResult reference; // for cluster dereplication 
	
	/**
	 * Instantiate a new domain from a HmmSearchResultAnnotation. 
	 * @param annotation	the annotation to convert 
	 */
	public Domain(final HmmSearchResultAnnotation annotation, DomainType type) {
		this.start = annotation.start();
		this.end = annotation.end();
		this.name = annotation.name();
		this.score = annotation.score();
		this.family = type.family();
		this.type = type;
	}
	
	/**
	 * Create a new deep copy of a domain. (Does not deep copy blastp search results!)
	 * @param domain	the domain to deep copy
	 */
	public Domain(final Domain domain) {
		this.start = domain.start();
		this.end = domain.end();
		this.score = domain.score();
		this.name = domain.name();
		this.family = domain.family();
		this.type = domain.type();
		this.substrates = domain.substrates();
	}
	
	/**
	 * Instantiate a new domain. 
	 * @param start		start of the domain
	 * @param end		end of the domain
	 * @param score		score associated with the domain 
	 * @param name		name of the domain
	 */
	public Domain(final int start, final int end, final double score, final String name) {
		this.start = start;
		this.end = end;
		this.score = score;
		this.name = name;
	}
		
	/**
	 * Set the sequence of this domain, given the start, end, and sequence of the parent open reading frame.
	 * @param sequence	the sequence of the parent orf
	 */
	public void setSequenceFromOrf(final String sequence) {
		this.sequence = sequence.substring(start, end - 1);
	}
	
	/**
	 * Print the relevant section of a domain's amino acid sequence, formatted
	 * for a multi-FASTA file (i.e. with trailing newline). Sequence is renamed
	 * to include domain start and end points.
	 * 
	 * @return the amino acid sequence in FASTA format
	 */
	public String printForFASTA() {
		return ">" + name + "\n" + sequence + "\n";
	}
	
	/**
	 * Get the point on the parent orf at which this domain starts.
	 * @return	domain start
	 */
	public int start() {
		return start;
	}
	
	/**
	 * Get the point on the parent orf at which this domain ends.
	 * @return	domain end 
	 */
	public int end() {
		return end;
	}
	
	/**
	 * Get the size of this domain in amino acids.
	 * @return	domain size 
	 */
	public int size() { 
		return end - start;
	}
	
	/**
	 * Get the name of this domain.
	 * @return	domain name
	 */
	public String name() {
		return name;
	}
	
	/**
	 * Get the score assigned to this domain by the process used to detect it.
	 * @return	domain score
	 */
	public double score() {
		return score;
	}
	
	/**
	 * Get the family of this domain, e.g. resistance or beta-lactam biosynthesis. 
	 * @return	domain family
	 */
	public DomainFamilies family() {
		return family;
	}
	
	/**
	 * Set the family of this domain.
	 * @param family	domain family 
	 */
	public void setFamily(DomainFamilies family) {
		this.family = family;
	}
	
	/**
	 * Get the family-specific type of this domain. 
	 * @return		the type of this domain
	 */
	public DomainType type() {
		return type;
	}
	
	/**
	 * Set the family-specific type of this domain.
	 * @param type	domain type
	 */
	public void setType(DomainType type) {
		this.type = type;
	}

	/**
	 * Add a new BLAST search result to this domain.
	 * @param result	the BLAST search result to add
	 */
	public void addBlastResult(final BlastSearchResult result) {
		blastResults.add(result);
	}
	
	/**
	 * Get all BLAST search results associated with this domain.
	 * @return	this domain's BLAST search results
	 */
	public List<BlastSearchResult> blastResults() {
		return blastResults;
	}

	/**
	 * Check whether another domain overlaps with this one.
	 * @param d2	the second domain
	 * @return		true if their start-end ranges overlap
	 */
	public boolean overlaps(final Domain d2) {
		final int x1 = start;
		final int x2 = end;
		final int y1 = d2.start();
		final int y2 = d2.end();
		return (x1 <= y2 && y1 <= x2);
	}
	
	/**
	 * Add a new substrate to this domain. 
	 * @param substrate	to add
	 */
	public void addSubstrate(Substrate substrate) {
		substrates.add(substrate);
	}
	
	/**
	 * Get all substrates associated with this domain.
	 * @return	all domain substrates
	 */
	public List<Substrate> substrates() {
		return substrates;
	}
	
	/**
	 * Set the entire list of substrates associated with this domain. 
	 * @param substrates	the new list of substrates
	 */
	public void setSubstrates(List<Substrate> substrates) {
		this.substrates = substrates;
	}
	
	/**
	 * Get the top-scoring substrate associated with this domain.
	 * @return	the top scoring substrate, or null if substrates is an empty list
	 */
	public Substrate topSubstrate() {
		return (substrates.size() > 0) ? substrates.get(0) : null;
	}

	/**
	 * If this is a glycosyltransferase domain, get its sugar substrate.
	 * @return	domain sugar substrate 
	 */
	public Sugar sugar() {
		return sugar;
	}
	
	/**
	 * If this is a glycosyltransferase domain, set its sugar substrate.
	 * @param sugar		domain sugar substrate 
	 */
	public void setSugar(Sugar sugar) {
		this.sugar = sugar;
	}

	/**
	 * Get the BLAST score of this domain to a reference domain of the same
	 * type, used to sort domains for cluster dereplication.
	 * 
	 * @return reference domain
	 */
	public BlastSearchResult referenceDomain() {
		return reference;
	}
	
	/**
	 * Get the score of the reference domain, or 0 if it is not set.
	 * 
	 * @return score or 0 if null
	 */
	public double referenceScore() {
		return reference == null ? 0.0d : reference.score();
	}

	/**
	 * Set the BLAST score of this domain to a reference domain of the same
	 * type, used to sort domains for cluster dereplication.
	 * 
	 * @param reference 	reference domain
	 */
	public void setReferenceDomain(BlastSearchResult reference) {
		this.reference = reference;
	}

}
