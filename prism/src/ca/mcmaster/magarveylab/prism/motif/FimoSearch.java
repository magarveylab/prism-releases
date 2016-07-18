package ca.mcmaster.magarveylab.prism.motif;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.List;

import ca.mcmaster.magarveylab.prism.util.PrismProcessBuilder;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Execute the FIMO executable to locate motifs within sequences.
 * 
 * @author skinnider
 *
 */
public class FimoSearch {
	
	/**
	 * The location of the FIMO executable file. 
	 */
	private String executable;
	
	/**
	 * The location of the motif file to search the domain sequence with. 
	 */
	private String motif;
	
	/**
	 * The location of the domain file to query against the motif.
	 */
	private String query;
	
	/**
	 * The current session. 
	 */
	private Session session;

	/**
	 * Instantiate a new FIMO search.
	 * 
	 * @param query
	 *            the query domain filepath
	 * @param motif
	 *            the ribosomal precursor motif filepath
	 * @param session
	 *            the current session
	 */
	public FimoSearch(String query, String motif, Session session) {
		this.query = query;
		this.motif = motif;
		this.session = session;
	}

	/**
	 * Execute a FIMO search, and associate all found motifs with their
	 * corresponding domains.
	 * 
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public List<FimoSearchResult> run() throws IOException, InterruptedException {
		// set the location of the executable file  
		setExecutable();
		
		// create command
		String[] cmd = { executable, motif, query };

		// execute command
		PrismProcessBuilder ppb = new PrismProcessBuilder(cmd);
		BufferedReader br = ppb.run();

		// read output
		FimoSearchReader reader = new FimoSearchReader(br);
		List<FimoSearchResult> results = reader.read();
		return results;
	}

	/**
	 * Set the location of the FIMO executable file.
	 * 
	 * @throws IOException
	 */
	public void setExecutable() throws IOException {
		// set executable
		executable = session.subDir("motifs") + "fimo.sh";
		// set executable permissions
		Runtime.getRuntime().exec("chmod +x " + executable);
	}
	
}