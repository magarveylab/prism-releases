package ca.mcmaster.magarveylab.prism.genome;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;

import ca.mcmaster.magarveylab.prism.util.PrismProcessBuilder;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Execute a hidden Markov model search against a protein FASTA database using hmmsearch (Finn et al. 2011).
 * @author skinnider
 *
 */
public class HmmSearch {

	protected String proteinFastaDatabase;
	protected String hmmModel;
	protected String hmmExecutable;
	protected Session session;
	private HmmSearchReader reader;
	
	/**
	 * Initiate a new hidden Markov model search against a protein FASTA database. 
	 * @param hmmModel	absolute path of the hidden Markov model file to search against; 
	 * 			output file will have the same name with the extension '.out'
	 * @param proteinFastaDatabase	the FASTA database to search
	 * @param session				the current PRISM session
	 */
	public HmmSearch(String hmmModel, String proteinFastaDatabase, Session session) {
		this.hmmModel = hmmModel;
		this.proteinFastaDatabase = proteinFastaDatabase;
		this.session = session;
	}

	/**
	 * Execute this hidden Markov model search. 
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public void run() throws IOException, InterruptedException {
		initialize();
		execute();
	}

	/**
	 * Set the locations of the hmmsearch executable and hmmsearch output. Set
	 * search reader.
	 * 
	 * @throws IOException
	 */
	protected void initialize() throws IOException {
		if (!new File(hmmModel).exists())
			throw new IOException(
					"Hidden Markov model file " + hmmModel + " does not exist");
		hmmExecutable = session.subDir("hmm") + "hmmsearch.sh";
		Runtime.getRuntime().exec("chmod +x " + hmmExecutable);
	}
	
	/**
	 * Create and start a process builder to execute the hmmsearch. 
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	protected void execute() throws IOException, InterruptedException {
		String[] cmd = { hmmExecutable, hmmModel, proteinFastaDatabase };
		PrismProcessBuilder ppb = new PrismProcessBuilder(cmd);
		BufferedReader br = ppb.run();
		reader = new HmmSearchReader(br);
	}
	
	/**
	 * Get the reader for this search.
	 * @return	the reader for this search
	 */
	public HmmSearchReader reader() {
		return reader;
	}
	
}
