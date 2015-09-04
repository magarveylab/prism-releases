package ca.mcmaster.magarveylab.prism.data;

import java.util.ArrayList;
import java.util.List;

/**
 * Container for data about the organism associated with a query genome, if this is available in the FASTA header.
 * @author prees, skinnider
 *
 */
public class Organism {

	private String accession = new String();
	private String genus = new String();
	private String species = new String();
	private String strain = new String();
	private String rawfastaheader = new String();
	private String contig = new String();
	private List<String> classification = new ArrayList<String>();
	private String sequenceType = new String();
	private String source = new String();
	
	/**
	 * Instantiate a new, empty organism. 
	 */
	public Organism() {
	}
	
	/**
	 * Instantiate an organism from a FASTA header.
	 * 
	 * @param header
	 *            the entire, raw header
	 * @param accession
	 *            the accession number
	 * @param genus
	 *            the organism's genus
	 * @param species
	 *            the organism's species
	 * @param strain
	 *            the organism's strain
	 * @param sequencetype
	 *            the type of sequence
	 * @param contig
	 *            the name of the parent contig
	 */
	public Organism(String header, String accession, String genus,
			String species, String strain, String sequencetype, String contig) {
		this.accession = accession;
		this.genus = genus;
		this.species = species;
		this.strain = strain;
		this.rawfastaheader = header;
		this.sequenceType = sequencetype;
		this.contig = contig;
	}
	
	public Organism(String header, String accession, String genus,
			String species, String strain, List<String> classification,
			String source) {
		this.accession = accession;
		this.genus = genus;
		this.species = species;
		this.strain = strain;
		this.rawfastaheader = header;
		this.classification = classification;
		this.source = source;
	}
	
	public Organism(String rawfastaheader) {
		this.rawfastaheader = rawfastaheader;
	}
	
	public String contig() {
		return contig;
	}
	
	public void setContig(String contig) {
		this.contig = contig;
	}
	
	public String sequenceType() {
		return sequenceType;
	}
	
	public void setSequenceType(String sequenceType) {
		this.sequenceType = sequenceType;
	}
	
	public String rawfastaheader() {
		return rawfastaheader;
	}
	
	public void setRawFastaHeader(String rawFastaHeader) {
		this.rawfastaheader = rawFastaHeader;
	}
	
	public String accession() {
		return accession;
	}
	
	public void setAccession(String accession) {
		this.accession = accession;
	}
	
	public String genus() {
		return genus;
	}
	
	public void setGenus(String genus) {
		this.genus = genus;
	}
	
	public String species() {
		return species;
	}
	
	public void setSpecies(String species) {
		this.species = species;
	}
	
	public String strain() {
		return strain;
	}
	
	public void setStrain(String strain) {
		this.strain = strain;
	}
	
	public List<String> classification() {
		return classification;
	}
	
	public void setClassification(List<String> classification) {
		this.classification = classification;
	}
	
	public String source() {
		return source;
	}

	public void setSource(String source) {
		this.source = source; 
	}
	
}
