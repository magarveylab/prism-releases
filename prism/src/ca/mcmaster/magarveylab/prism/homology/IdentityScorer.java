package ca.mcmaster.magarveylab.prism.homology;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.prism.blast.BlastpSearch;
import ca.mcmaster.magarveylab.prism.blast.BlastSearchResult;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.homology.data.HomologousCluster;
import ca.mcmaster.magarveylab.prism.util.Sorter;
import ca.mcmaster.magarveylab.prism.util.Strings;
import ca.mcmaster.magarveylab.wasp.session.Session;

public class IdentityScorer {

	private String query;
	private String database;
	private Cluster cluster;
	private Session session;
	
	/**
	 * Instantiate a new cluster BLASTp search against a precompiled cluster BLAST database. 
	 * @param database	location of the database
	 * @param query		location of the query file
	 * @param cluster	the cluster this search is BLASTing against
	 * @param session	the current session
	 */
	public IdentityScorer(String database, String query, Cluster cluster, Session session) {
		this.database = database;
		this.query = query;
		this.cluster = cluster;
		this.session = session;
	}
	
	/**
	 * Execute a BLASTp search.
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public void run() throws IOException, InterruptedException {		
		BlastpSearch blast = new BlastpSearch(database, query, session);
		List<BlastSearchResult> results = blast.run(-1);

		List<BlastSearchResult> filtered = doPreferentialMatching(results);
		double score = score(filtered);

		String name = Strings.name(database);
		for (HomologousCluster hc : cluster.homologs())
			if (hc.name().equals(name))
				hc.setIdentityScore(score);
	}
	
	/**
	 * Restrict a list of BLASTP search results, in which all possible query-subject matches are read, to the highest-
	 * scoring matches of a unique query to a unique subject. 
	 * @param results
	 * @return
	 */
	private List<BlastSearchResult> doPreferentialMatching(List<BlastSearchResult> results) {
		List<BlastSearchResult> filtered = new ArrayList<BlastSearchResult>();
		List<String> queries = new ArrayList<String>();
		List<String> subjects = new ArrayList<String>();
		
		Sorter.sortBlastpResults(results);
		for (BlastSearchResult result : results) {
			String query = result.query();
			String subject = result.subject();
			if (!Strings.contains(query, queries) && !Strings.contains(subject, subjects)) {
				filtered.add(result);
				queries.add(query);
				subjects.add(subject);
			}
		}
	
		return filtered;
	}
	
	/**
	 * Compute the weighted identity index, I, of a homologous cluster.
	 * @param results	BlastpSearchResults of the homologous cluster 
	 * @return			identity index, I 
	 */
	public double score(List<BlastSearchResult> results) {
		double numerator = 0.0d;
		double denominator = 0.0d;
		for (BlastSearchResult result : results) {
			double identity = result.identity();
			double length = 1.0d * result.length();
			double coverage = 1.0d * result.coverage();
			
			numerator += identity * coverage;
			denominator += length;
		}
		double score = numerator / denominator;
		if (Double.isNaN(score)) 
			score = 0.0;
		
		return score;
	}

}
