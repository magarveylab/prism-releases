package ca.mcmaster.magarveylab.prism.homology;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ca.mcmaster.magarveylab.enums.domains.BetaLactamDomains;
import ca.mcmaster.magarveylab.enums.domains.DeoxySugarDomains;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.PrerequisiteDomains;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.enums.domains.TypeIIPolyketideDomains;
import ca.mcmaster.magarveylab.prism.blast.BlastpSearch;
import ca.mcmaster.magarveylab.prism.blast.BlastSearchResult;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.fasta.FastaUtil;
import ca.mcmaster.magarveylab.prism.fasta.FastaWriter;
import ca.mcmaster.magarveylab.prism.util.Sorter;
import ca.mcmaster.magarveylab.prism.util.Strings;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Use a MLST-style cluster dereplication algorithm, first presented by Udwary,
 * to calculate dereplication scores to known clusters.
 * 
 * @author skinnider
 *
 */
public class DomainScorer {

	private Cluster cluster;
	private Session session;

	/**
	 * Instantiate a new domain scorer for a putative biosynthetic gene cluster.
	 * @param cluster	query cluster identified by PRISM
	 * @param session	current session 
	 */
	public DomainScorer(Cluster cluster, Session session) {
		this.cluster = cluster;
		this.session = session;
	}
	
	/**
	 * Score a putative biosynthetic gene cluster against the MLST-type cluster database. 
	 * @return			list of BLAST results for each database cluster 
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public List<BlastSearchResult> score() throws IOException, InterruptedException {
		String file = session.dir() + "cluster_" + cluster.index() + "_domains.fasta";
		sortAndWriteDomains(cluster, file);
			
		// blast against database 
		session.listener().updateLastDetail("Calculating cluster " + cluster.index() + " homology scores...");
		String database = session.subDir("domainFasta") + "all";
		BlastpSearch blast = new BlastpSearch(database, file, session);
		List<BlastSearchResult> results = blast.run(-1);

		// parse results
		Sorter.sortBlastpResults(results);
		return results;
	}

	/**
	 * Write a concatenated sequence of domains to a file in order to
	 * dereplicate this cluster.
	 * 
	 * @param cluster
	 *            query cluster identified by PRISM
	 * @param path
	 *            filepath to write the file to
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public void sortAndWriteDomains(Cluster cluster, String path) throws IOException, InterruptedException {
		session.listener().updateLastDetail("Calculating reference domain scores for cluster " 
				+ cluster.index() + "...");
		List<Domain> sorted = new ArrayList<Domain>();
		
		// get domain types
		List<DomainType> types = new ArrayList<DomainType>();
		types.addAll(Arrays.asList(ThiotemplatedDomains.values()));
		types.addAll(Arrays.asList(TailoringDomains.values()));
		types.addAll(Arrays.asList(TypeIIPolyketideDomains.values()));
		types.addAll(Arrays.asList(DeoxySugarDomains.values()));
		types.addAll(Arrays.asList(PrerequisiteDomains.values()));
		types.addAll(Arrays.asList(BetaLactamDomains.values()));
		
		// sort and add domains
		for (DomainType type : types) {
			List<Domain> domains = cluster.domains(type);
			sortDomains(domains, type);
			sorted.addAll(domains);
		}
		
		// write sorted list to file 
		StringBuffer sb = new StringBuffer();
		for (Domain domain : sorted) {
			String sequence = FastaUtil.getDomainSequence(domain, cluster);
			if (sequence != null) {
				sb.append(sequence);
			} else {
				throw new NullPointerException("No sequence for domain " + domain.name() + "!");
			}
		}
		
		String name = Strings.name(path);
		FastaWriter.printToFasta(name, sb.toString(), path);
	}

	public void sortDomains(List<Domain> domains, DomainType type) throws IOException, InterruptedException {
		String database = session.subDir("reference") + type.family().toString().toLowerCase() + File.separator
				+ type.toString().toLowerCase();
		
		if (new File(database + ".phr").exists()
				&& domains.size() > 0) {
			String query = session.dir() + "cluster_" + cluster.index() + "_" + type.toString().toLowerCase() + ".fa";
			FastaWriter.printDomainsToFasta(domains, query);

			BlastpSearch blast = new BlastpSearch(database, query, session);
			List<BlastSearchResult> results = blast.run(-1);
			for (BlastSearchResult result : results) 
				for (Domain domain : domains) 
					if (result.query().equals(domain.name())) 
						domain.setReferenceDomain(result);
			Sorter.sortDomainsByReferenceScore(domains);
		} else {
			Sorter.sortDomainsByScore(domains);
		}
	}

}
