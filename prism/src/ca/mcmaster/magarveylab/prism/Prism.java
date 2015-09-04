package ca.mcmaster.magarveylab.prism;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import org.openscience.cdk.exception.CDKException; 

import ca.mcmaster.magarveylab.prism.cluster.ClusterFinder;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.LibraryGenerator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Contig;
import ca.mcmaster.magarveylab.prism.data.Genome;
import ca.mcmaster.magarveylab.prism.data.Library;
import ca.mcmaster.magarveylab.prism.genome.GenomeSearch;
import ca.mcmaster.magarveylab.prism.orfs.OrfSearch;
import ca.mcmaster.magarveylab.prism.homology.HomologousClusterSearch;
import ca.mcmaster.magarveylab.prism.tanimoto.TanimotoSearch;
import ca.mcmaster.magarveylab.prism.tanimoto.data.TanimotoScore;
import ca.mcmaster.magarveylab.prism.util.PrismFileWriter;
import ca.mcmaster.magarveylab.prism.util.exception.BadSequenceException;
import ca.mcmaster.magarveylab.prism.util.exception.BadSmilesToFingerprinterException;
import ca.mcmaster.magarveylab.prism.util.exception.DatabaseConnectException;
import ca.mcmaster.magarveylab.prism.util.exception.NoSequenceException;
import ca.mcmaster.magarveylab.prism.util.exception.SugarGeneException;
import ca.mcmaster.magarveylab.prism.web.PrismConfig;
import ca.mcmaster.magarveylab.prism.web.html.PrismReport;
import ca.mcmaster.magarveylab.wasp.WebApplication;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * PRISM: PRediction Informatics for Secondary Metabolomes. <br>
 * Search a genome for novel genetically encoded secondary metabolites.
 * 
 * @author skinnider
 *
 */
public class Prism implements Runnable, WebApplication {

	private boolean terminated;
	private Session session;
	private PrismConfig config;
	private Genome genome;

	/**
	 * Instantiate a new PRISM search.
	 * 
	 * @param config
	 *            PRISM search configuration
	 * @param session
	 *            current session
	 */
	public Prism(PrismConfig config, Session session) {
		this.config = config;
		this.session = session;
	}
	
	/**
	 * Execute this PRISM search.
	 */
	public void run() {
		try {
			findOrfs();
			analyzeOrfs();
			findClusters();
			generateScaffolds();
			writeFiles();
			if (config.score)
				score();
		} catch (Exception e) {
			session.exceptionHandler().throwException(e);
		} finally {
			terminate();
		}
	}

	/**
	 * Find all orfs within a genome, and detect the organism associated with
	 * the genome file.
	 * 
	 * @throws Exception
	 * @throws InterruptedException
	 * @throws IOException
	 */
	public void findOrfs() throws IOException, InterruptedException, Exception {
		session.listener().addStage("Identifying genes",
				"Finding open reading frames...");
		
		File file = new File(config.input);
		System.out.println("config file initialized"); //TODO: remove after debugging
		if (file.length() > 50 * 1024 * 1024) // 50 MB check 
			throw new BadSequenceException("Could not run PRISM: "
					+ "your file exceeds the maximum file size of 50 MB!");
		
		genome = new Genome(file);
		System.out.println("genome variable initialized"); //TODO: remove after debugging
		OrfSearch os = new OrfSearch(genome, session);
		System.out.println("initialzed orf search"); //TODO: remove after debugging
		os.run();
		System.out.println("started orf search"); //TODO: remove after debugging
		System.out.println("Finished finding orfs.");
		if (genome == null || genome.contigs().size() == 0)
			throw new NoSequenceException("Could not detect any sequence "
					+ "information in user-submitted file!");
	}
	
	/**
	 * Use PRISM's library of BLAST databases and hidden Markov models to detect
	 * biosynthetic machinery within the genome.
	 * 
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws SugarGeneException
	 */
	public void analyzeOrfs() throws IOException, InterruptedException,
			SugarGeneException, Exception {
		System.out.println("Started analyzeOrfs()..");
		GenomeSearch gs = new GenomeSearch(genome, session);
		gs.run();
	}
	
	/**
	 * Group orfs into biosynthetic clusters.
	 * 
	 * @throws IOException
	 */
	public void findClusters() throws IOException {
		session.listener().addStage("Detecting clusters", "Clustering biosynthetic orfs...");

		ClusterFinder cf = new ClusterFinder();
		for (Contig contig : genome.contigs())
			cf.cluster(contig, config, session);

		System.out.println("[Prism] Detected " + genome.clusters().size() + " clusters");
	}

	/**
	 * Generate cluster scaffold libraries.
	 * 
	 * @throws Exception
	 */
	public void generateScaffolds() throws Exception {
		session.listener().addStage("Predicting products", 
				"Generating predicted scaffolds...");
			
		for (Contig contig : genome.contigs()) {
			for (Cluster cluster : contig.clusters()) { 
				try {
					LibraryGenerator lg = new LibraryGenerator(cluster, session);
					lg.generate();	
				} catch (Exception e) {
					session.listener().addStage("Structure prediction error<br>Cluster #" + cluster.index(), 
							"PRISM failed to generate a scaffold for cluster #" + cluster.index());
					e.printStackTrace();
				}
			}
		}
	}
	
	/**
	 * Write files required for HTML report.
	 * @throws Exception 
	 */
	public void writeFiles() throws Exception {
		session.listener().addStage("Writing files", "Generating output...");
		PrismFileWriter.writeAllFiles(genome, session);
	}

	/**
	 * Assess chemical similarity and genetic homology of identified clusters to known small molecules.
	 * @throws IOException
	 * @throws CDKException
	 * @throws SQLException 
	 * @throws BadSmilesToFingerprinterException 
	 * @throws InterruptedException 
	 * @throws DatabaseConnectException 
	 */
	public void score() throws IOException, CDKException, SQLException, BadSmilesToFingerprinterException, 
			InterruptedException, DatabaseConnectException {
		session.listener().addStage("Identifying known clusters", "Calculating homology to known clusters...");
		scoreHomology();
		scoreChemicalSimilarity();
	}
	
	/**
	 * Score clusters' genetic homology to known biosynthetic gene clusters. 
	 * @throws IOException
	 * @throws InterruptedException 
	 * @throws NumberFormatException 
	 */
	public void scoreHomology() throws IOException, NumberFormatException, InterruptedException {
		for (Contig contig : genome.contigs()) 
			for (Cluster cluster : contig.clusters()) {
				System.out.println("[Prism] Calculating cluster " + cluster.index() + " homology to known clusters");
				HomologousClusterSearch cs = new HomologousClusterSearch(cluster, session);
				cs.run();
			}
	}
	
	/**
	 * Score predicted cluster products' chemical similarity to known small molecules. 
	 * @throws IOException
	 * @throws CDKException
	 * @throws SQLException
	 * @throws InterruptedException
	 * @throws DatabaseConnectException 
	 */
	public void scoreChemicalSimilarity() throws IOException, CDKException, SQLException, InterruptedException, 
			DatabaseConnectException {
		for (Contig contig : genome.contigs()) {
			for (Cluster cluster : contig.clusters()) {
				session.listener().updateLastDetail("Calculating cluster " + cluster.index() + " Tanimoto coefficients...");
								
				Library library = cluster.library();
				if (library != null && library.scaffolds().size() > 0) {
					try {
						List<TanimotoScore> scores = TanimotoSearch.scoreFingerprints(cluster, session);
						cluster.scores().addAll(scores);
					} catch (BadSmilesToFingerprinterException e) {
						session.listener().updateLastDetail("Error generating fingerprints for cluster " + cluster.index());
						continue;
					}
				}
			}
		}
	}
	
	/**
	 * Check if this search is terminated.
	 * @return	true if the search has been terminated, false if it has not
	 */
	public boolean isTerminated() {
		return terminated;
	}
	
	/**
	 * Terminate this Prism genome search.
	 */
	public void terminate() {
		try {
			PrismReport report = (PrismReport) session.report();
			report.terminate();
	    	session.listener().addStage("Done", "Analysis complete!");
	    	Thread.sleep(1001);
			report.writeClusterPages();
			terminated = true;
			report.cancel();
		} catch (Exception e) {
			System.out.println("[Prism] Error writing report");
			session.exceptionHandler().throwException(e); 
		} finally {
			System.out.println("[Prism] Terminated.");
		}
	}
		
	/**
	 * Get the genome object created by this PRISM search. 
	 * @return	genome detected by this PRISM search
	 */
	public Genome genome() {
		return genome; 
	}
	
	/**
	 * Set the genome object associated with this PRISM search (for reading JSON). 
	 * @param genome	genome associated with this PRISM search
	 */
	public void setGenome(Genome genome) {
		this.genome = genome;
	}
	
	/**
	 * Get the configuration of this PRISM search.
	 * @return	PRISM search configuration
	 */
	public PrismConfig config() {
		return config;
	}
	
	/**
	 * Get the current session. 
	 * @return	current session 
	 */
	public Session session() {
		return session;
	}
	
}
