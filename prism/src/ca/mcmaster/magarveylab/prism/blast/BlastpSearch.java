package ca.mcmaster.magarveylab.prism.blast;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.List;

import ca.mcmaster.magarveylab.prism.util.PrismProcessBuilder;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Execute a BLASTp search.
 * @author skinnider
 *
 */
public class BlastpSearch {
	
	private String query;
	private String database;
	private String eValue = "1E-05";
	private Session session;
	
	/**
	 * Instantiate a new BLASTp search against a precompiled BLAST database. 
	 * @param database	location of the database
	 * @param query		location of the query file
	 * @param session	the current session
	 */
	public BlastpSearch(String database, String query, Session session) {
		this.database = database;
		this.query = query;
		this.session = session;
	}

	/**
	 * Instantiate a new BLAST search against a precompiled BLAST database with a specific E-value cutoff. 
	 * @param database	location of the database
	 * @param query		location of the query file
	 * @param session	the current session
	 * @param evalue	the desired e-value 
	 */
	public BlastpSearch(String database, String query, Session session, String evalue) {
		this.database = database;
		this.query = query;
		this.session = session;
		this.eValue = evalue;
	}

	/**
	 * Execute a BLASTp search.
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public List<BlastSearchResult> run(int resultSize) throws IOException, InterruptedException {		
		// set executable csh script, chmod
		String executable = session.subDir("blast") + "blastp.sh";
		Runtime.getRuntime().exec("chmod +x " + executable);
		
		// use: blastn -db nt -query sequence.fsa -out results.out 
		String[] cmd = { executable, database, query, eValue };
		PrismProcessBuilder ppb = new PrismProcessBuilder(cmd);
		BufferedReader br = ppb.run();
		
		BlastSearchReader reader = new BlastSearchReader(br);
		List<BlastSearchResult> results = reader.read(resultSize);
		return results;
	}

}
