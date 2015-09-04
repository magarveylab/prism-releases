package ca.mcmaster.magarveylab.prism.data;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import ca.mcmaster.magarveylab.prism.util.Sorter;

/**
 * Generic container for any contiguous sequence, consisting of a header and a sequence.
 * @author skinnider
 *
 */
public class Contig {

	private int index;
	private String header;
	private String sequence;
	private List<Orf> orfs = new ArrayList<Orf>();
	private List<Cluster> clusters = new ArrayList<Cluster>();
	private HashMap<String,String> files = new HashMap<String,String>();
	
	/**
	 * Instantiate a new contig.
	 * @param header	the FASTA or GenBank header for the contig sequence
	 * @param sequence	the contig sequence 
	 */
	public Contig(String header, String sequence) {
		this.header = header;
		this.sequence = sequence;
	}
	
	/**
	 * Get the header of this contig.
	 * @return	contig header
	 */
	public String header() {
		return header;
	}
	
	/**
	 * Get the nucleotide sequence of this contig.
	 * @return	nucleotide sequence 
	 */
	public String sequence() {
		return sequence;
	}
	
	/**
	 * Set the nucleotide sequence of this contig.
	 * @param sequence	nucleotide sequence
	 */
	public void setSequence(String sequence) {
		this.sequence = sequence;
	}
	
	/**
	 * Get all open reading frames on this contig. 
	 * @return	all contig orfs
	 */
	public List<Orf> orfs() {
		return orfs;
	}
	
	/**
	 * Set the open reading frames on this contig.
	 * @param orfs	contig orfs 
	 */
	public void setOrfs(List<Orf> orfs) {
		this.orfs = orfs;
	}
	
	/**
	 * Get all clusters on this contig.
	 * @return	contig clusters
	 */
	public List<Cluster> clusters() {
		return clusters;
	}
	
	/**
	 * Get the length of this contig's nucleotide sequence. 
	 * @return	nucleotide sequence length 
	 */
	public int length() {
		return sequence.length();
	}
	
	/**
	 * Get the index of this contig within the user-submitted sequence file.
	 * @return	contig index
	 */
	public int index() {
		return index;
	}

	/**
	 * Set the index of this contig within the user-submitted sequence file.
	 * @param index	contig index
	 */
	public void setIndex(int index) {
		this.index = index;
	}
	
	/**
	 * Associate a new filepath with this contig. 
	 * @param key	the key to access this filepath with
	 * @param value	the filepath itself
	 */
	public void setFile(String key, String value) {
		files.put(key, value);
	}
	
	/**
	 * Get a filepath associated with this contig.
	 * @param key	the key to access this filepath with
	 * @return
	 */
	public String getFile(String key) {
		return files.get(key);
	}

	/**
	 * Get all the orfs on this contig within the boundaries of a given
	 * cluster, regardless of whether they contain biosynthetic domains.
	 * 
	 * @param cluster 	the cluster in question
	 * @return 			all open reading frames within the boundaries of this cluster
	 */
	public List<Orf> getAllOrfs(Cluster cluster) {
		List<Orf> orfs = new ArrayList<Orf>();
	
		// get cluster start & end
		int start = cluster.start();
		int end = cluster.end();
	
		// if orf is within the cluster, add to return list
		for (Orf orf : this.orfs) {
			if (orf.start() > start && orf.end() < end)
				orfs.add(orf);
		}
	
		Sorter.sortOrfs(orfs);
		return orfs;
	}
	
}
