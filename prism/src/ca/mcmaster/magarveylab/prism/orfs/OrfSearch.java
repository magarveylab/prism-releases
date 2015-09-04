package ca.mcmaster.magarveylab.prism.orfs;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava3.core.exceptions.CompoundNotFoundError;

import ca.mcmaster.magarveylab.prism.Prism;
import ca.mcmaster.magarveylab.prism.data.Contig;
import ca.mcmaster.magarveylab.prism.data.Genome;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.data.Organism;
import ca.mcmaster.magarveylab.prism.fasta.SequenceCleaner;
import ca.mcmaster.magarveylab.prism.fasta.FastaReader;
import ca.mcmaster.magarveylab.prism.fasta.FastaUtil;
import ca.mcmaster.magarveylab.prism.fasta.FastaWriter;
import ca.mcmaster.magarveylab.prism.util.GenbankParser;
import ca.mcmaster.magarveylab.prism.util.Sorter;
import ca.mcmaster.magarveylab.prism.util.exception.BadSequenceException;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Find open reading frames within a user-input sequence file, and return a new genome object.
 * @author skinnider
 *
 */
public class OrfSearch {
	
	private static final int minOrfSize = 75;
	
	private int counter = 1;
	private Genome genome;
	private Session session;

	private final Logger logger = Logger.getLogger(OrfSearch.class.getName());

	/**
	 * Instantiate a new open reading frame search.
	 * @param file		user-input sequence file
	 * @param session	current PRISM session
	 */
	public OrfSearch(Genome genome, Session session) {
		this.genome = genome;
		this.session = session;
	}

	/**
	 * Detect all open reading frames within a user-input sequence file. 
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws BadSequenceException 
	 * @throws Exception
	 */
	public void run() throws IOException, InterruptedException, BadSequenceException, Exception {
		String filename = genome.filepath();
		if (filename.matches(".+\\.gbk$") || filename.matches(".+\\.gb$") || filename.matches(".+\\.genbank$")) {
			logger.log(Level.INFO, "Set input file format to GenBank");
			GenbankParser.readGenbankFile(genome, filename, session);
			trimContigs(genome);
			findAllOrfs(genome);
		} else {
			logger.log(Level.INFO, "Set input file format to FASTA");
			FastaReader.readFastaFile(genome);
			trimContigs(genome);
			findAllOrfs(genome);
			configureFastaGenome(genome);
		}
	}
	
	/**
	 * Remove excess contigs within a user-submitted genome file.
	 * @param genome	genome to analyze
	 */
	public void trimContigs(Genome genome) {
		List<Contig> contigs = genome.contigs();
		if (contigs.size() >= 1_000) {
			session.listener().addStage("Trimming contigs", "Removing excess contigs from sequence file with " 
					+ contigs.size() + " contigs: maximum size 1,000 contigs!");
			Sorter.sortContigsBySize(contigs);
			contigs.subList(999, contigs.size() - 1).clear();
		}
	}
	
	/**
	 * Configure a genome input in FASTA format by detecting the organism, or simply setting
	 * the organism to be equal to the header of the first FASTA item. 
	 * @param genome
	 * @throws IOException
	 * @throws BadSequenceException
	 */
	private void configureFastaGenome(Genome genome) throws IOException, BadSequenceException {
		Prism prism = (Prism) session.webapp();
		
		List<Contig> contigs = genome.contigs();
		if (contigs.size() > 0) {
			Contig first = contigs.get(0);
			String header = first.header();
			Organism organism = null;
			if (prism.config().parseorganism) {
				organism = FastaUtil.parseOrganism(header);
			} else {
				organism = new Organism(header);
			}
			genome.setOrganism(organism);
		}
	}
	
	/**
	 * Find all open reading frames within the complete set of contigs in a user-submitted genome
	 * or other sequence file.
	 * @param genome	genome to analyze 
	 * @throws IOException
	 * @throws BadSequenceException
	 */
	private void findAllOrfs(Genome genome) throws IOException, BadSequenceException {
		List<Contig> contigs = genome.contigs();
		try {
			for (Contig contig : contigs) {
				List<Orf> orfs = findOrfs(contig);
				logger.log(Level.INFO, "Found " + orfs.size() + " orfs within contig " 
						+ (contigs.indexOf(contig) + 1) + " of " + contigs.size());
				contig.setOrfs(orfs);	
			}	
			
			// set files for each contig
			for (Contig contig : genome.contigs())
				FastaUtil.setFiles(contig, session);
			
			// print each contig's orfs to file
			for (Contig contig : genome.contigs())
				if (contig.orfs().size() > 0)
					FastaWriter.printOrfsToFasta(contig.orfs(), contig.getFile("orfs"));
		} catch (CompoundNotFoundError e) {
			throw new BadSequenceException();
		}
	}
	
	/**
	 * Manually find all orfs above a minimum size for a single contig. 
	 * @param contig	contig in question
	 * @return			All orfs in the contig
	 * @throws IOException 
	 */
	public List<Orf> findOrfs(Contig contig) throws CompoundNotFoundError {
		// clean sequence
		String clean = SequenceCleaner.clean(contig.sequence());
		contig.setSequence(clean);
		
		// find orfs
		List<Orf> orfs = scanForOrfs(contig);
		return orfs;
	}
		
	/**
	 * Manually find all orfs with length > 60 nt within a single contig. 
	 * @param contig		the contig to scan
	 * @param orfCounter	counter for total number of orfs identified
	 * @return				a list of all orfs in that contig
	 */
	public List<Orf> scanForOrfs(Contig contig) throws CompoundNotFoundError {
		List<Orf> orfs = new ArrayList<Orf>();
		String sequence = contig.sequence();
		int length = sequence.length();
		String reverseComplementarySequence = SequenceConverter.getReverseComplementaryDNASequence(sequence);
		List<Orf> forwardOrfs = scanFrameForOrfs(sequence);
		List<Orf> reverseOrfs = scanFrameForOrfs(reverseComplementarySequence);
		for (Orf orf : forwardOrfs)
			orf.setFrame("+");
		for (Orf orf : reverseOrfs) {
			orf.setFrame("-");
			int start = length - orf.start();
			int end = length - orf.end();
			orf.setStart(end);
			orf.setEnd(start);
		}
		orfs.addAll(forwardOrfs);
		orfs.addAll(reverseOrfs);
		Sorter.sortOrfs(orfs);
		return orfs;
	}
	
	/**
	 * Manually find all orfs with length > 60 nt within a single frame (forward or backward).
	 * @param sequence	the sequence to scan
	 * @param counter	counter for total number of orfs identified
	 * @return			a list of all orfs in that frame
	 */
	public List<Orf> scanFrameForOrfs(String sequence) throws CompoundNotFoundError {
		List<Orf> orfs = new ArrayList<Orf>();
		String[] startCodons = { "ATG", "GTG", "TTG" };
		String[] stopCodons = { "TAA", "TAG", "TGA" };
		int[] lastStop = { 0, 0, 0 };
		boolean doContinue = false;
		for (int i = 0; i < sequence.length() - 2; i++) {
			String codon = sequence.substring(i, i+3);
			if (Arrays.asList(startCodons).contains(codon) || i == 0) {
				// if index is less than last stop codon in same frame, a longer orf exists; ignore this substring
				if (i < lastStop[i % 3])
					continue;
				doContinue = false;
				for (int j = i; j < sequence.length() - 2; j += 3) {
					String secondCodon = sequence.substring(j, j+3);
					if (Arrays.asList(stopCodons).contains(secondCodon)) {
						if (j - i > minOrfSize) { // reject orfs shorter than 60 nt
							String orfSequence = sequence.substring(i,j+2);
							String orfAASequence = SequenceConverter.convertDNAToAA(orfSequence);
							String orfIndex = String.format("%05d", counter);
							String orfName = "orf" + orfIndex;
							Orf orf = new Orf(orfName, orfAASequence);
							orf.setStart(i);
							orf.setEnd(j+2);
							orfs.add(orf);
							counter++;
							lastStop[j%3] = j;
							doContinue = true;
							break;
						}
					}
				}
				if (doContinue)
					continue;
				// if a final start codon is found with no stop codon, save from start to end of sequence as an orf
				String orfSequence = sequence.substring(i);
				String orfAASequence = SequenceConverter.convertDNAToAA(orfSequence);
				String orfIndex = String.format("%05d", counter);
				String orfName = "orf" + orfIndex;
				Orf orf = new Orf(orfName, orfAASequence);
				orf.setStart(i);
				orf.setEnd(sequence.length());
				orfs.add(orf);
				counter++;
				lastStop[i%3] = sequence.length() + 1;
				break;
			}
		}
		return orfs;
	}
	
}
