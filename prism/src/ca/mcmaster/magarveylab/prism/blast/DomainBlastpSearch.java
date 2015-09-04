package ca.mcmaster.magarveylab.prism.blast;

import java.io.IOException;
import java.util.List;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.prism.Prism;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Contig;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.web.PrismConfig;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Execute a BLASTp search on a thiotemplated domain-specific database.  
 * @author skinnider
 *
 */
public class DomainBlastpSearch {
	
	private String query;
	private String database;
	private Session session;
	private DomainType type;
	private Contig contig;
	private PrismConfig config;
	
	/**
	 * Instantiate a new BLASTP search against a precompiled BLAST database. 
	 * @param database	location of the database
	 * @param query		location of the query file
	 * @param type		the type of domain that this search is BLASTing against
	 * @param contig	the current contig
	 * @param session	the current session
	 */
	public DomainBlastpSearch(String database, String query, DomainType type, Contig contig, Session session) {
		this.database = database;
		this.query = query;
		this.type = type;
		this.contig = contig;
		this.session = session;
		Prism prism = (Prism) session.webapp();
		this.config = prism.config();
	}
	
	/**
	 * Run a thiotemplated domain-specific BLASTP search. 
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void run() throws IOException, InterruptedException {
		BlastpSearch blast = new BlastpSearch(database, query, session);
		List<BlastSearchResult> results = blast.run(config.display);
		matchResultsToDomains(results);
	}
	
	/**
	 * Match BLASTp results to their corresponding domains.
	 * @param results	list of results to match
	 */
	private void matchResultsToDomains(List<BlastSearchResult> results) {
		for (BlastSearchResult blastpResult : results) {
			for (Orf orf : contig.orfs()) 
				for (Domain domain : orf.domains(type))
					if (domain.name().equals(blastpResult.query()) && domain.score() > type.cutoff())
						domain.addBlastResult(blastpResult);
		}
	}

}
