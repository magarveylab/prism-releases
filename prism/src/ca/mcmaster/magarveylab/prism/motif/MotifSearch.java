package ca.mcmaster.magarveylab.prism.motif;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.List;

import ca.mcmaster.magarveylab.enums.RibosomalPrecursorMotifs;
import ca.mcmaster.magarveylab.prism.cluster.analysis.OrfAnalyzer;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.util.exception.MotifException;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Search a domain for the occurrence of a motif.
 * 
 * @author skinnider
 *
 */
public class MotifSearch {

	private Domain domain;
	private RibosomalPrecursorMotifs type;
	private Cluster cluster;
	private Session session;

	/**
	 * Instantiate a new motif search given a query domain, a ribosomal
	 * precursor motif, and a cluster.
	 * 
	 * @param domain
	 *            the query domain
	 * @param type
	 *            the ribosomal precursor motif to search for in the query
	 *            domain
	 * @param cluster
	 *            the parent cluster
	 */
	public MotifSearch(Domain domain, RibosomalPrecursorMotifs type,
			Cluster cluster, Session session) {
		this.domain = domain;
		this.type = type;
		this.cluster = cluster;
		this.session = session;
	}

	/**
	 * Execute a FIMO search, and associate all found motifs with their
	 * corresponding domains.
	 * 
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void run() throws IOException, InterruptedException {
		String queryFilepath = session.dir() + "cluster_" + cluster.index()
				+ "_" + type.name() + ".fasta";
		String motifFilepath = session.subDir("motifs" + File.separator
				+ "leader_peptides")
				+ type.getFilename();

		writeQueryFile(queryFilepath);

		FimoSearch search = new FimoSearch(queryFilepath, motifFilepath,
				session);
		List<FimoSearchResult> results = search.run();

		matchResultsToDomains(results);
	}

	/**
	 * Write ribosomal precursor domains to a file to be read by FIMO.
	 * 
	 * @param queryFilepath
	 *            location to write the precursor domain file
	 * @throws IOException
	 */
	private void writeQueryFile(String queryFilepath) throws IOException {
		FileOutputStream is = new FileOutputStream(queryFilepath);
		OutputStreamWriter osw = new OutputStreamWriter(is);
		Writer w = new BufferedWriter(osw);
		w.write(">" + domain.name() + "\n"
				+ OrfAnalyzer.getParentOrf(domain, cluster).sequence());
		w.close();
	}

	/**
	 * Match FIMO search results to the query domain and instantiate new motif
	 * objects.
	 * 
	 * @param results
	 *            the results of the FIMO search
	 */
	private void matchResultsToDomains(List<FimoSearchResult> results) {
		// get domain sequence
		Orf orf = OrfAnalyzer.getParentOrf(domain, cluster);
		String query = orf.sequence();

		for (FimoSearchResult result : results) {
			try {
				Motif motif = new Motif(type, query, result);
				domain.addMotif(motif);
			} catch (MotifException e) {
				e.printStackTrace();
				continue;
			}
		}
	}

}
