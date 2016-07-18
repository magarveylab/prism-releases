package ca.mcmaster.magarveylab.prism.genome;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.interfaces.SubstrateType;
import ca.mcmaster.magarveylab.prism.Prism;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Contig;
import ca.mcmaster.magarveylab.prism.data.Genome;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.data.Substrate;
import ca.mcmaster.magarveylab.prism.enums.hmms.SubstrateHmm;
import ca.mcmaster.magarveylab.prism.fasta.FastaWriter;
import ca.mcmaster.magarveylab.prism.genome.data.HmmSearchResult;
import ca.mcmaster.magarveylab.prism.genome.data.HmmSearchResultAnnotation;
import ca.mcmaster.magarveylab.prism.util.Sorter;
import ca.mcmaster.magarveylab.prism.web.PrismConfig;
import ca.mcmaster.magarveylab.wasp.session.Session;

public class SubstrateDomainSearch {
	
	protected Genome genome;
	protected Session session;
	protected PrismConfig config;
	protected String model;
	protected DomainType type;
	protected SubstrateHmm[] substrates;
	protected String family;
	protected String typeString;
	
	public SubstrateDomainSearch(DomainType type, SubstrateHmm[] substrates, Genome genome, Session session) {
		this.genome = genome;
		this.session = session;
		Prism prism = (Prism) session.webapp();
		this.config = prism.config();
		this.type = type;
		this.substrates = substrates;
		this.model = getModelFile(type);
		this.family = type.family().toString().toLowerCase();
		this.typeString = type.toString().toLowerCase();
	}
	
	public void run() {
		try {
			findDomains();
			findSubstrates();	
			removeExcessSubstrates();
			removeOrphanDomains();
		} catch (Exception e) {
			session.exceptionHandler().throwException(e);
		}
	}

	
	/**
	 * Get the location of the hidden Markov model file.
	 * @param type	domain type being searched for 
	 * @return		location of the hidden Markov model file
	 */
	public String getModelFile(DomainType type) {
		return session.subDir("hmm") + family + File.separator + type.hmm();
	}
	
	/**
	 * Search biosynthetic orfs for putative domains. 
	 * @param contigs			FASTA item to analyze
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	private void findDomains() throws IOException, InterruptedException {
		session.listener().addStage("Analyzing " + type.fullName().toLowerCase() + " domains", 
				"Identifying " + type.fullName().toLowerCase() + " domains..."); 
		
		for (Contig contig : genome.contigs()) {
			HmmSearch search = new HmmSearch(session.subDir("hmm") + family + File.separator + type.hmm(), 
					contig.getFile("orfs"), session);
			search.run();

			List<HmmSearchResult> results = search.reader().read();
			matchDomainsToOrfs(contig, results);

			printDomainsToFasta(contig, contig.getFile(typeString));
		}
	}

	/**
	 * Locate putative substrate regions within substrate domains.  
	 * @param contigs			FASTA item to analyze
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	private void findSubstrates() throws IOException, InterruptedException {
		// Run ensemble of substrate HMMs on domain FASTA file
		for (SubstrateHmm substrate : substrates) {
			if (substrate.fullName().toString().length() > 0)
				session.listener().updateLastDetail("Finding " + type.fullName().toLowerCase()
						+ " domains with substrate " + substrate.fullName() + "...");

			for (Contig contig : genome.contigs()) {
				HmmSearch hmm = new HmmSearch(session.subDir("hmm") + family
						+ File.separator + typeString + File.separator + substrate.hmm(), 
						contig.getFile(typeString), session);
				hmm.run();

				List<HmmSearchResult> substrateResults = hmm.reader().read();
				matchSubstratesToDomains(contig, substrateResults, substrate);
				
				sortSubstrates(contig);
			}
		}
	}
	
	/**
	 * Matches a list of hidden Markov model search results representing domains to their corresponding orfs.
	 * @param contig	contig to analyze
	 * @param results	a list of domain HmmSearchResult objects 
	 */
	private void matchDomainsToOrfs(Contig contig, List<HmmSearchResult> results) {
		for (HmmSearchResult result : results) 
			for (Orf orf : contig.orfs()) 
				if (orf.toString().equals(result.name) || orf.name().equals(result.name)) 
					parseResultForDomains(orf, result);
	}
	
	/**
	 * Instantiate a new domain for each HmmSearchResultAnnotation} in this HmmSearchResult;
	 * set the sequence from the orf; and add it to the orf's scaffold. 
	 * @param orf				the orf in question
	 * @param result	the HmmSearchResult to parse
	 */
	public void parseResultForDomains(Orf orf, HmmSearchResult result) {
		for (HmmSearchResultAnnotation annotation : result.annotations()) {
			Domain d = new Domain(annotation, type); 
			if (d.score() > type.cutoff()) {
				String orfSequence = orf.sequence();
				String domainSequence = orfSequence.substring(
						d.start() - 1, d.end() - 1);
				d.setSequence(domainSequence);
				orf.domains().add(d);				
			}
		}
	}

	/**
	 * Print all domains to a FASTA file. 
	 * @param contig	contig to analyze
	 * @param path	the path of the FASTA file
	 * @throws IOException
	 */
	private void printDomainsToFasta(Contig contig, String path) throws IOException {
		List<Domain> domains = new ArrayList<Domain>();
		for (Orf orf : contig.orfs()) {
			List<Domain> typeDomains = orf.domains(type);
			domains.addAll(typeDomains);
		}

		if (domains.size() == 0)
			return;

		FastaWriter.printDomainsToFasta(domains, path);
	}

	/**
	 * Matches a list of hidden Markov model search results representing domain substrates to their corresponding 
	 * domains.
	 * @param contig	contig to analyze
	 * @param results	a list of HmmSearchResults representing substrates
	 * @param substrate			the substrate in question
	 */
	private void matchSubstratesToDomains(Contig contig, List<HmmSearchResult> results, 
			SubstrateType substrate) {
		for (HmmSearchResult result : results) 
			for (Orf orf : contig.orfs())
				// Match HmmSearchResult to the appropriate orf
				// (substrateResult.name() is orf.toString() w/domain start & end appended in this context)
				if (result.name.indexOf(orf.toString()) != -1 || result.name.indexOf(orf.name()) != -1) 
					// next match HmmSearchResult to appropriate domain
					for (Domain domain : orf.domains(type)) 
						if (result.name.indexOf(domain.name()) != -1) 
							for (HmmSearchResultAnnotation annotation : result.annotations()) {
								Substrate s = new Substrate(annotation, substrate);
								if (s.score() > type.cutoff()) 
									domain.addSubstrate(s);
							}
	}
	
	/**
	 * Sort identified substrates by score.
	 * @param contig		contig to analyze
	 */
	private void sortSubstrates(Contig contig) {
		for (Orf orf : contig.orfs()) {
			List<Domain> domains = orf.domains(type);
			for (Domain domain : domains)
				Sorter.sortSubstrates(domain.substrates());
		}
	}
	
	/**
	 * Remove domains with one or fewer substrates from the results. 
	 */
	private void removeOrphanDomains() {
		for (Contig contig : genome.contigs())
			for (Orf orf : contig.orfs()) 
				for (Domain domain : orf.domains(type))
					if (domain.substrates().size() < 1 || domain.topSubstrate().score() < type.cutoff())
						orf.domains().remove(domain);
	}

	/**
	 * Clear substrates which are beyond the display limit.  
	 */
	private void removeExcessSubstrates() {
		for (Contig contig : genome.contigs())
			for (Orf orf : contig.orfs())
				for (Domain domain : orf.domains(type)) {
					List<Substrate> substrates = domain.substrates();
					int start = config.display - 1;
					int end = substrates.size() - 1;
					if (end > start) 
						substrates.subList(start, end).clear();
				}
	}

}
