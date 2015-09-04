package ca.mcmaster.magarveylab.prism.genome;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import ca.mcmaster.magarveylab.enums.DomainFamilies;
import ca.mcmaster.magarveylab.enums.SubstrateDomainSearches;
import ca.mcmaster.magarveylab.enums.domains.AminoglycosideDomains;
import ca.mcmaster.magarveylab.enums.domains.BetaLactamDomains;
import ca.mcmaster.magarveylab.enums.domains.DeoxySugarDomains;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.NucleosideDomains;
import ca.mcmaster.magarveylab.enums.domains.PrerequisiteDomains;
import ca.mcmaster.magarveylab.enums.domains.RegulatorDomains;
import ca.mcmaster.magarveylab.enums.domains.ResistanceDomains;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.enums.domains.TypeIIPolyketideDomains;
import ca.mcmaster.magarveylab.enums.substrates.SubstrateType;
import ca.mcmaster.magarveylab.prism.Prism;
import ca.mcmaster.magarveylab.prism.blast.DomainBlastpSearch;
import ca.mcmaster.magarveylab.prism.data.Contig;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Genome;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.data.RnaSequence;
import ca.mcmaster.magarveylab.prism.fasta.FastaWriter;
import ca.mcmaster.magarveylab.prism.util.RibosomalSequence;
import ca.mcmaster.magarveylab.prism.util.Sorter;
import ca.mcmaster.magarveylab.prism.util.exception.RibosomalSequenceException;
import ca.mcmaster.magarveylab.prism.util.exception.SugarGeneException;
import ca.mcmaster.magarveylab.prism.web.PrismConfig;
import ca.mcmaster.magarveylab.wasp.session.Session;

public class GenomeSearch {

	private PrismConfig config;
	private Genome genome;
	private Session session;
	
	/**
	 * List of domains which can overlap with other domains. This includes
	 * domains which may be found inside other domains (e.g., the RiPP precursor
	 * recognition element), or domains for which the boundary is unclear (e.g.,
	 * fungal NR-PKS SAT domains).
	 */
	public final static DomainType[] canOverlap = new DomainType[] { 
		RibosomalDomains.RRE,
		ThiotemplatedDomains.REDUCTASE,
	};

	/**
	 * Instantiate a new genome analysis search.
	 * @param genome	the genome to analyze
	 * @param session	the current session
	 */
	public GenomeSearch(Genome genome, Session session) {
		this.genome = genome;
		this.session = session;
		Prism prism = (Prism) session.webapp();
		config = prism.config();
	}

	/**
	 * Detect adenylation domains and their substrates; acyltransferase domains and their substrates; and run
	 * other hidden Markov model domain searches within putatively biosynthetic orfs of a single contig. 
	 * @throws IOException 
	 * @throws InterruptedException 
	 * @throws SugarGeneException 
	 */
	public void run() throws InterruptedException, SugarGeneException, IOException, Exception {
		findRibosomalSequence();
		
		findSubstrateDomains();
		findPrerequisiteDomains();
		findOtherThiotemplatedDomains();
		findTypeIIPolyketideDomains();
		findTailoringDomains();
		executeBlastpAnalysis();
		findSugarGenes();
		findBetaLactamGenes();
		findAminoglycosideGenes();
		findResistanceDomains();
		findRibosomalDomains();
		findRegulatorDomains();
		findNucleosideGenes();

		removeOverlap();
	}
	
	/**
	 * Executed only if the find16s field in the config object is set to true. 
	 * Implements barrnap to search for 16s sequences in sequences where no 16s sequences are manually annotated.
	 * @throws IOException, InterruptedException, Exception
	 */
	public void findRibosomalSequence() throws IOException, InterruptedException, Exception {
		Prism prism = (Prism) session.webapp();
		PrismConfig config = prism.config();
		
		try {
			if (config.find16s) {
				RibosomalSequence rs = new RibosomalSequence(prism.genome().contigs(), config.input);
				rs.runProcess();
				List<RnaSequence> allRnaSeq = rs.parseRawOutput();
				//if there were not sequences added from the Genbank file and if there are some found from prediction
				if (prism.genome().ribosomalSequences().size() < 1 && !(allRnaSeq == null)) {
					genome.addRibosomalSequences(allRnaSeq);
				}
			}
		} catch (IOException e) {
			throw new RibosomalSequenceException("Could not run barrnap to identify 16S sequences!"); 
		}
	}

	/**
	 * Find substrate-containing domains within this genome.
	 */
	private void findSubstrateDomains() {
		for (SubstrateDomainSearches searchType : SubstrateDomainSearches.values()) {
			DomainType domain = searchType.type();
			SubstrateType[] substrates = searchType.substrates();

			SubstrateDomainSearch search = new SubstrateDomainSearch(domain, substrates, genome, session);
			search.run();
		}
	}
	
	/**
	 * Find type II polyketide domains within this genome.
	 */
	private void findTypeIIPolyketideDomains() {
		session.listener().addStage("Analyzing type II polyketide domains", 
				"Finding type II poleyktide domains...");
		for (TypeIIPolyketideDomains type : TypeIIPolyketideDomains.values()) {
			DomainSearch search = new DomainSearch(type, genome, session);
			search.run();
		}
	}

	/**
	 * Find tailoring domains within this genome.
	 */
	private void findTailoringDomains() {
		session.listener().addStage("Analyzing tailoring domains", 
				"Finding tailoring domains...");
		for (TailoringDomains type : TailoringDomains.values()) {
			DomainSearch search = new DomainSearch(type, genome, session);
			search.run();
		}
	}
	
	/**
	 * Find non-scaffold thiotemplated biosynthetic machinery within this genome.
	 * @throws IOException 
	 */
	private void findOtherThiotemplatedDomains() throws IOException {
		session.listener().addStage("Analyzing other thiotemplated domains", 
				"Finding other domains...");
		
		typeLoop:
		for (ThiotemplatedDomains type : ThiotemplatedDomains.values()) {
			// not substrate domains
			for (SubstrateDomainSearches searchType : SubstrateDomainSearches.values())
				if (type == searchType.type())
					continue typeLoop;	

			DomainSearch hmm = new DomainSearch(type, genome, session);
			hmm.run();
		}
	}

	/**
	 * Run BLAST analysis on domains identified within this genome. (Currently implemented: condensation, C-starter, and 
	 * ketosynthase).
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	private void executeBlastpAnalysis() throws IOException, InterruptedException {
		session.listener().updateLastDetail("Analyzing condensation and ketosynthase domains...");
		for (Contig contig : genome.contigs()) {
			// NAPDOS analysis
			DomainBlastpSearch condensationBlast = new DomainBlastpSearch(session.subDir("blast") + "condensationdb", 
					contig.getFile("condensation"), ThiotemplatedDomains.CONDENSATION, contig, session);
			condensationBlast.run();
			DomainBlastpSearch ketosynthaseBlast = new DomainBlastpSearch(session.subDir("blast") + "ketosynthasedb", 
					contig.getFile("ketosynthase"), ThiotemplatedDomains.KETOSYNTHASE, contig, session);
			ketosynthaseBlast.run();

			// type II polyketide analysis 
			DomainBlastpSearch ksbBlast = new DomainBlastpSearch(session.subDir("blast") + "CLF", 
					contig.getFile("CLF"), TypeIIPolyketideDomains.CLF, contig, session);
			ksbBlast.run();

			// c starter analysis
			FastaWriter.printCStarter(contig, contig.getFile("cstarter"));
			DomainBlastpSearch cstarterBlast = new DomainBlastpSearch(session.subDir("blast") + "cstarter", 
					contig.getFile("cstarter"), ThiotemplatedDomains.CONDENSATION, contig, session);
			cstarterBlast.run();
			
			// priming AT analysis 
			DomainBlastpSearch atBlast = new DomainBlastpSearch(session.subDir("blast") + "priming_AT", 
					contig.getFile("priming_AT"), TypeIIPolyketideDomains.PRIMING_AT, contig, session);
			atBlast.run();
			Sorter.sortBlastpResults(contig);
		}

		session.listener().updateLastDetail("Analyzing halogenation and glycosyltransferase domains...");
		for (Contig contig : genome.contigs()) {
			// PTM analysis
			DomainBlastpSearch chlorinaseBlast = new DomainBlastpSearch(session.subDir("blast") + "chlorination", 
					contig.getFile("chlorination"), TailoringDomains.CHLORINATION, contig, session);
			chlorinaseBlast.run();
			DomainBlastpSearch glycosyltransferaseBlast = new DomainBlastpSearch(session.subDir("blast") + "glycosyltransferase", 
					contig.getFile("glycosyltransferase"), TailoringDomains.GLYCOSYLTRANSFERASE, contig, session);
			glycosyltransferaseBlast.run();
		}
	}

	/**
	 * Find ribosomal domains within this genome.
	 */
	private void findRibosomalDomains() {
		if (config.ribosomal) {
			session.listener().addStage("Analyzing ribosomal domains", 
					"Finding ribosomal domains...");
			for (RibosomalDomains type : RibosomalDomains.values()) {
				DomainSearch hmm = new DomainSearch(type, genome, session);
				hmm.run();
			}
		}
	}

	/**
	 * Find biosynthetic domains which are considered prerequisites for rare substrate detection within this genome.
	 */
	private void findPrerequisiteDomains() {
		session.listener().addStage("Analyzing substrate biosynthesis domains", 
				"Finding rare substrate biosynthesis domains...");
		for (PrerequisiteDomains type : PrerequisiteDomains.values()) {
			DomainSearch hmm = new DomainSearch(type, genome, session);
			hmm.run();
		}
	}

	/**
	 * Find resistance genes within this genome.
	 */
	private void findResistanceDomains() {
		if (config.resistance) {
			session.listener().addStage("Analyzing resistance domains", 
					"Finding resistance domains...");
			for (ResistanceDomains type : ResistanceDomains.values()) {
				DomainSearch hmm = new DomainSearch(type, genome, session);
				hmm.run();
			}
		}
	}

	/**
	 * Find resistance genes within this genome.
	 */
	private void findRegulatorDomains() {
		if (config.regulation) {
			session.listener().addStage("Analyzing regulator domains", 
					"Finding regulator domains...");
			for (RegulatorDomains type : RegulatorDomains.values()) {
				DomainSearch hmm = new DomainSearch(type, genome, session);
				hmm.run();
			}
		}
	}

	/**
	 * Find sugar biosynthesis genes within this genome. 
	 * @throws IOException
	 * @throws SugarGeneException
	 */
	private void findSugarGenes() {
		session.listener().addStage("Analyzing sugar genes", 
				"Finding sugar biosynthesis genes...");
		for (DeoxySugarDomains gene : DeoxySugarDomains.values()) {
			DomainSearch sgs = new DomainSearch(gene, genome, session);
			sgs.run();
		}
	}
	
	/**
	 * Find beta-lactam subfamily biosynthesis genes within this genome.
	 */
	private void findBetaLactamGenes() {
		session.listener().addStage("Analyzing beta-lactam genes", 
				"Finding beta-lactam and monobactam biosynthesis genes...");
		for (BetaLactamDomains d : BetaLactamDomains.values()) {
			DomainSearch s = new DomainSearch(d, genome, session);
			s.run();
		}
	}
	
	/**
	 * Find aminoglycoside subfamily biosynthesis genes within this genome.
	 */
	private void findAminoglycosideGenes() {
		if (config.aminoglycoside) {
			session.listener().addStage("Analyzing aminoglycoside genes", 
					"Finding aminoglycoside biosynthesis genes...");
			for (AminoglycosideDomains d : AminoglycosideDomains.values()) {
				DomainSearch s = new DomainSearch(d, genome, session);
				s.run();
			}
		}
	}
	
	/**
	 * Find nucleoside subfamily biosynthesis genes within this genome.
	 */
	private void findNucleosideGenes() {
		if (config.nucleoside) {
			session.listener().addStage("Analyzing nucleoside genes", 
					"Finding nucleoside biosynthesis genes...");
			for (NucleosideDomains d : NucleosideDomains.values()) {
				DomainSearch s = new DomainSearch(d, genome, session);
				s.run();
			}
		}
	}

	/**
	 * Remove overlapping domains.
	 */
	private void removeOverlap() {
		for (Contig contig : genome.contigs()) {
			for (Orf orf : contig.orfs()) {
				List<Domain> domains = orf.domains();
				for (int i = 0; i < domains.size(); i++)
					for (int j = 0; j < domains.size(); j++) {
						Domain d1 = domains.get(i);
						Domain d2 = domains.get(j);
						if (d1.overlaps(d2) && d1 != d2
								&& Arrays.asList(canOverlap).indexOf(d1.type()) == -1
								&& d1.family() != DomainFamilies.REGULATOR
								&& Arrays.asList(canOverlap).indexOf(d2.type()) == -1
								&& d2.family() != DomainFamilies.REGULATOR)
							if (d2.score() > d1.score()) {
								domains.remove(d1);
								if (i > 0)
									i--;
								System.out.println("Removed " + d1.name());
							}
					}
			}
		}
	}

}
