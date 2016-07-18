package ca.mcmaster.magarveylab.prism.orfs;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava3.core.exceptions.CompoundNotFoundError;

import ca.mcmaster.magarveylab.enums.FileType;
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
import ca.mcmaster.magarveylab.prism.util.Files;
import ca.mcmaster.magarveylab.prism.util.exception.BadSequenceException;
import ca.mcmaster.magarveylab.prism.util.exception.ProdigalSearchException;
import ca.mcmaster.magarveylab.prism.web.PrismConfig;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Find open reading frames within a user-input sequence file, and return a new
 * genome object.
 * 
 * @author skinnider
 *
 */
public class OrfSearch {
	
	/**
	 * The minimum length of an open reading frame, in nucleotides. 
	 */
	private static final int minOrfSize = 60;
	
	private int counter = 1;
	private Genome genome;
	private PrismConfig config;
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
		Prism prism = (Prism) session.webapp();
		this.config = prism.config();
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

		// read contigs into genome object 
		FileType type = Files.getFileType(filename);
		if (type == FileType.GENBANK) {
			logger.log(Level.INFO, "Set input file format to GenBank");
			GenbankParser.readGenbankFile(genome, filename, session.root(), session.dir());
		} else if (type == FileType.FASTA) {
			logger.log(Level.INFO, "Set input file format to FASTA");
			FastaReader.readFastaFile(genome);
		}
		trimContigs(genome);

		// get orfs
		List<GenePredictionModes> modes = config.genePredictionModes;
		if (modes.contains(GenePredictionModes.PRODIGAL)) {
			try {
				ProdigalSearch prodigal = new ProdigalSearch(genome, session);
				prodigal.run();
			} catch (ProdigalSearchException e) {
				session.listener().addStage("Prodigal error",
								e.getMessage() + "<br>Finding all possible coding sequences...");
				if (!modes.contains(GenePredictionModes.ALL_ORFS))
					modes.add(GenePredictionModes.ALL_ORFS);
			}
		}
		
		if (modes.contains(GenePredictionModes.ALL_ORFS)) 
			// find  all potential coding sequences 
			findAllOrfs(genome);
		
		// remove overlap between different modes 
		removeOverlap();
		
		// set files for each contig
		for (Contig contig : genome.contigs())
			FastaUtil.setFiles(contig, session);

		// print each contig's orfs to file
		for (Contig contig : genome.contigs())
			if (contig.orfs().size() > 0)
				FastaWriter.printOrfsToFasta(contig.orfs(),
						contig.getFile("orfs"));
		
		// set organism information from FASTA genome  
		if (type == FileType.FASTA)
			configureFastaGenome(genome);
	}
	
	/**
	 * Remove excess contigs within a user-submitted genome file.
	 * @param genome	genome to analyze
	 */
	public void trimContigs(Genome genome) {
		List<Contig> contigs = genome.contigs();
		if (contigs.size() >= 1000) {
			session.listener().addStage("Trimming contigs", "Removing excess contigs from sequence file with " 
					+ contigs.size() + " contigs: maximum size 1,000 contigs!");
			Sorter.sortContigsBySize(contigs);
			contigs.subList(999, contigs.size() - 1).clear();
		}
	}
	
	/**
	 * Configure a genome input in FASTA format by detecting the organism, or
	 * simply setting the organism to be equal to the header of the first FASTA
	 * item.
	 * 
	 * @param genome
	 * @throws IOException
	 * @throws BadSequenceException
	 */
	private void configureFastaGenome(Genome genome) throws IOException,
			BadSequenceException {		
		List<Contig> contigs = genome.contigs();
		if (contigs.size() > 0) {
			Contig first = contigs.get(0);
			String header = first.header();
			Organism organism = null;
			if (config.parseorganism) {
				organism = FastaUtil.parseOrganism(header);
			} else {
				organism = new Organism(header);
			}
			genome.setOrganism(organism);
		}
	}
	
	/**
	 * Remove overlapping orfs detected by different gene prediction methods.
	 * Currently, orfs found by Prodigal supersede orfs detected by scanning an
	 * input file for all potential coding sequences. This can easily be
	 * extended in the future to include other open reading frame detection
	 * methods (e.g. for fungal genomes).
	 */
	private void removeOverlap() {
		// set priorities
		// Prodigal > all orfs
		Map<GenePredictionModes,Integer> priorities = new HashMap<GenePredictionModes,Integer>();
		priorities.put(GenePredictionModes.ALL_ORFS, 0);
		priorities.put(GenePredictionModes.PRODIGAL, 1);
		
		for (Contig contig : genome.contigs()) {
			List<Orf> orfs = contig.orfs();
			Iterator<Orf> itr1 = orfs.iterator();
			loop1: while (itr1.hasNext()) {
				Orf orf1 = itr1.next();
				Iterator<Orf> itr2 = orfs.iterator();
				while (itr2.hasNext()) {
					Orf orf2 = itr2.next();
					int start1 = orf1.frame().equals("+") ? orf1.start() : orf1.end();
					int start2 = orf2.frame().equals("+") ? orf2.start() : orf2.end();
					if (orf1 != orf2 && orf1.overlaps(orf2)
							&& orf1.frame().equals(orf2.frame())
							&& start1 % 3 == start2 % 3) {
						String partial1 = orf1.getPartial();
						String partial2 = orf2.getPartial();
						int priority1 = priorities.get(orf1.getMode());
						int priority2 = priorities.get(orf2.getMode());
						// remove partial orfs that overlap with complete orfs 
						if (partial1 != null && partial1.contains("1")
								&& partial2 != null && !partial2.contains("1")) {
//							System.out.println("Removed " + orf1.name() + " ("
//									+ orf1.start() + "-" + orf1.end()
//									+ "): orf is partial and overlaps "
//									+ orf2.name() + " (" + orf2.start() + "-"
//									+ orf2.end() + ")");
							itr1.remove();
							continue loop1;
						}
						// remove lower priority method orfs 
						if (priority2 > priority1) {
//							System.out.println("Removed " + orf1.name() + " ("
//									+ orf1.start() + "-" + orf1.end()
//									+ ") with frame " + orf1.frame() + " for overlap with " + orf2.name()
//									+ " (" + orf2.start() + "-" + orf2.end()
//									+ ") with frame " + orf1.frame());
							itr1.remove();
							continue loop1;
						}
					}
				}
			}
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
				logger.log(Level.INFO, "Found " + orfs.size() + " potential coding sequences within contig " 
						+ (contigs.indexOf(contig) + 1) + " of " + contigs.size());
				contig.addOrfs(orfs);	
			}	
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
							String orfSequence = sequence.substring(i,j+3);
							String orfAASequence = SequenceConverter.convertDNAToAA(orfSequence);
							String orfIndex = String.format("%05d", counter);
							String orfName = "orf" + orfIndex;
							Orf orf = new Orf(orfName, orfAASequence);
							orf.setStart(i);
							orf.setEnd(j+3);
							orf.setMode(GenePredictionModes.ALL_ORFS);
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
				orf.setMode(GenePredictionModes.ALL_ORFS);
				orfs.add(orf);
				counter++;
				lastStop[i%3] = sequence.length() + 1;
				break;
			}
		}
		return orfs;
	}
	
}
