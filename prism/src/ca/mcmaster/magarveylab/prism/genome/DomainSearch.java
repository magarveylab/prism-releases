package ca.mcmaster.magarveylab.prism.genome;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.prism.Prism;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Contig;
import ca.mcmaster.magarveylab.prism.data.Genome;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.fasta.FastaWriter;
import ca.mcmaster.magarveylab.prism.genome.data.HmmSearchResult;
import ca.mcmaster.magarveylab.prism.genome.data.HmmSearchResultAnnotation;
import ca.mcmaster.magarveylab.prism.util.Sorter;
import ca.mcmaster.magarveylab.prism.web.PrismConfig;
import ca.mcmaster.magarveylab.wasp.session.Session;

public class DomainSearch implements Runnable {
	
	protected Genome genome;
	protected Session session;
	protected PrismConfig config;
	protected HmmSearch search;
	protected String model;
	protected DomainType type;
	
	public DomainSearch(DomainType type, Genome genome, Session session) {
		this.genome = genome;
		this.session = session;
		Prism prism = (Prism) session.webapp();
		this.config = prism.config();
		this.type = type;
		this.model = getModelFile(type);
	}
	
	/**
	 * Get the location of the hidden Markov model file.
	 * @param type	domain type being searched for 
	 * @return		location of the hidden Markov model file
	 */
	public String getModelFile(DomainType type) {
		String family = type.family().toString().toLowerCase();
		return session.subDir("hmm") + family + File.separator + type.hmm();
	}
	
	public void run() {
		if (new File(model).exists()) {
			try {
				setStatusMessage();
				for (Contig contig : genome.contigs()) {
					runHmmSearch(contig);
					matchResultsToOrfs(contig);
					removeOverlap(contig);
					sortDomains(contig);
					printDomains(contig);
				}
			} catch (Exception e) {
				session.exceptionHandler().throwException(e);
			}
		}
	}
	
	/**
	 * Set the status message shown to the user during this domain search. 
	 */
	protected void setStatusMessage() {
		session.listener().updateLastDetail("Identifying " + type.fullName() + " domains...");
	}

	/**
	 * Execute the hmmsearch executable program given the arguments for this domain search.
	 * @param contig		contig to analyze
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	protected void runHmmSearch(Contig contig) throws IOException, InterruptedException {
		search = new HmmSearch(model, contig.getFile("orfs"), session);
		search.run();
	}
	
	/**
	 * Match results parsed from hmmsearch output text files to Orf objects within PRISM. 
	 * @param contig		contig to analyze
	 * @throws IOException
	 */
	protected void matchResultsToOrfs(Contig contig) throws IOException {
		List<HmmSearchResult> results = search.reader().read();
		for (HmmSearchResult result : results) {
			for (Orf orf : contig.orfs()) {
				if (orf.toString().equals(result.name) || orf.name().equals(result.name)) {
					parseResultForDomains(orf, result);
				}
			}
		}
	}
	
	/**
	 * Remove overlapping domains.
	 * @param contig	contig to analyze
	 */
	protected void removeOverlap(Contig contig) {
		for (Orf orf : contig.orfs())
			for (Domain d1 : orf.domains(type))
				for (Domain d2 : orf.domains(type)) 
					if (d1 != d2 && d1.overlaps(d2) && d1.score() < d2.score())
						orf.domains().remove(d1);
	}

	/**
	 * Sort all domains resulting from this domain search.
	 * @param contig		contig to analyze
	 */
	protected void sortDomains(Contig contig) {
		for (Orf orf : contig.orfs())
			Sorter.sortDomains(orf.domains());
	}
		
	/**
	 * Instantiate a new domain for each HmmSearchResultAnnotation in this HmmSearchResult. 
	 * @param type	the type of domain
	 * @param domainResult	the HmmSearchResult to parse
	 */
	protected void parseResultForDomains(Orf orf, HmmSearchResult domainResult) {
		for (HmmSearchResultAnnotation annotation : domainResult.annotations()) {
			Domain domain = new Domain(annotation, type);
			if (domain.score() > type.cutoff()) {
				domain.setSequenceFromOrf(orf.sequence());
				orf.domains().add(domain);				
			}
		}
	}

	/**
	 * Sort the substrates for this domain, if there are any. 
	 * @param contig	contig to analyze
	 */
	protected void sortSubstrates(Contig contig) {
		for (Orf orf : contig.orfs()) {
			List<Domain> domains = orf.substrateDomains();
			for (Domain domain : domains) 				
				Sorter.sortSubstrates(domain.substrates());
		}
	}
	
	/**
	 * Print all domains resulting from this domain search (on a given contig) to a file.
	 * @param contig		current contig being analyzed
	 * @throws IOException
	 */
	protected void printDomains(Contig contig) throws IOException {
		String path = session.dir() + "contig_" + contig.index() + "_" + type.toString().toLowerCase() + ".fa";
		List<Domain> domains = new ArrayList<Domain>();
		for (Orf orf : contig.orfs()) {
			List<Domain> orfDomains = orf.domains(type);
			domains.addAll(orfDomains);
		}

		if (domains.size() == 0)
			return;

		FastaWriter.printDomainsToFasta(domains, path);
	}

}
