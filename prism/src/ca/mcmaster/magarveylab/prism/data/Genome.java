package ca.mcmaster.magarveylab.prism.data;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.prism.web.html.graph.CircularGenomeGraph;

/**
 * Data package for a user-submitted sequence file.
 * 
 * @author skinnider
 *
 */
public class Genome {

	private File file;
	private Organism organism;
	private CircularGenomeGraph graph;
	private List<Contig> contigs = new ArrayList<Contig>();
	private List<RnaSequence> ribosomalSequences = new ArrayList<RnaSequence>();

	/**
	 * Instantiate a new genome.
	 * 
	 * @param file
	 *            sequence file from which this genome is derived.
	 */
	public Genome(File file) {
		this.file = file;
	}

	/**
	 * Instantiate a new genome.
	 * 
	 * @param contigs
	 *            list of contigs that make up the genome.
	 */
	public Genome(List<Contig> contig) {
		this.contigs = contig;
	}

	/**
	 * Get the sequence file associated with this genome.
	 * 
	 * @return the genome sequence file
	 */
	public File file() {
		return file;
	}

	/**
	 * Set the sequence file associated with this genome.
	 * 
	 * @return the genome sequence file
	 */
	public void setFile(File f) {
		this.file = f;
	}

	/**
	 * Get the filepath of this genome's sequence file.
	 * 
	 * @return the filepath of the user-submitted sequence file
	 */
	public String filepath() {
		return file.getPath();
	}

	/**
	 * Get the filename of this genome's sequence file.
	 * 
	 * @return the filename of the user-submitted sequence file
	 */
	public String filename() {
		return file.getName();
	}

	/**
	 * Get the ribosomal sequences within this genome.
	 * 
	 * @return the genome's ribosomal sequences
	 */
	public List<RnaSequence> ribosomalSequences() {
		return this.ribosomalSequences;
	}

	/**
	 * Set the ribosomal sequences within this genome.
	 * 
	 * @param allRnaSeq
	 *            the genome's ribosomal sequences
	 */
	public void setRibosomalSequences(List<RnaSequence> allRnaSeq) {
		this.ribosomalSequences = allRnaSeq;
	}

	/**
	 * Add ribosomal sequences to this genome.
	 * 
	 * @param allRnaSeq
	 *            a set of ribosomal sequences
	 */
	public void addRibosomalSequences(List<RnaSequence> allRnaSeq) {
		this.ribosomalSequences.addAll(allRnaSeq);
	}

	/**
	 * Get all of the contigs associated with this genome.
	 * 
	 * @return all contigs
	 */
	public List<Contig> contigs() {
		return contigs;
	}

	/**
	 * Set the contigs associated with this genome.
	 * 
	 * @param contigs
	 *            all contigs
	 */
	public void setContigs(List<Contig> contigs) {
		this.contigs = contigs;
	}

	/**
	 * Get the circular graph of this genome.
	 * 
	 * @return the genome's circular graph
	 */
	public CircularGenomeGraph graph() {
		return graph;
	}

	/**
	 * Set the circular graph of this genome.
	 * 
	 * @param graph
	 *            the genome's circular graph
	 */
	public void setGraph(CircularGenomeGraph graph) {
		this.graph = graph;
	}

	/**
	 * Get the organism that has been associated with this genome.
	 * 
	 * @return the genome's organism
	 */
	public Organism organism() {
		return organism;
	}

	/**
	 * Set the organism to be associated with this genome.
	 * 
	 * @param organism
	 *            the genome's organism
	 */
	public void setOrganism(Organism organism) {
		this.organism = organism;
	}

	/**
	 * Get all clusters within this genome.
	 * 
	 * @return all secondary metabolism clusters in the genome
	 */
	public List<Cluster> clusters() {
		List<Cluster> clusters = new ArrayList<Cluster>();
		for (Contig f : contigs)
			clusters.addAll(f.clusters());
		return clusters;
	}

}
