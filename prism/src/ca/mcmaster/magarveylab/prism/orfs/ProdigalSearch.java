package ca.mcmaster.magarveylab.prism.orfs;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import ca.mcmaster.magarveylab.prism.data.Contig;
import ca.mcmaster.magarveylab.prism.data.Genome;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.util.PrismProcessBuilder;
import ca.mcmaster.magarveylab.prism.util.exception.ProdigalSearchException;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Use Prodigal to identify putative genes within a user-submitted sequence. <br>
 * <br>
 * See Prodigal website at http://prodigal.ornl.gov/ or Prodigal
 * documentation at https://github.com/hyattpd/prodigal/wiki/.<br>
 * <br>
 * Prodigal usage:<br>
 * -i INPUT<br>
 * -a PROTEIN_OUT<br>
 * -d NUCLEOTIDE_OUT<br>
 * -f OUTPUT_FORMAT
 * -o OUTPUT_FILE
 * @author skinnider
 *
 */
public class ProdigalSearch {
	
	private Genome genome;
	private Session session;
	private static final Logger logger = Logger.getLogger(ProdigalSearch.class.getName());

	/**
	 * Instantiate a new open reading frame search using Prodigal.<br>
	 * <br>
	 * See Prodigal website at http://prodigal.ornl.gov/ or Prodigal
	 * documentation at https://github.com/hyattpd/prodigal/wiki/
	 * 
	 * @param genome
	 *            current genome 
	 * @param session
	 *            current PRISM session
	 */
	public ProdigalSearch(Genome genome, Session session) {
		this.genome = genome;
		this.session = session;
	}

	public void run() throws IOException, InterruptedException,
			ProdigalSearchException {
		// get mode 
		String mode = getProdigalMode();
		
		String file = genome.filepath();
		String output = session.dir() + "prodigal_out.txt"; // for debugging
		String protein = session.dir() + "orfs_protein.fasta";
		String nucleotide = session.dir() + "orfs_nucleotide.fasta";
		
		String executable = session.subDir("orfs") + "prodigal.sh";
		Runtime.getRuntime().exec("chmod +x " + executable);
		
		String[] cmd = { executable, mode, output, file, protein, nucleotide };
		PrismProcessBuilder ppb = new PrismProcessBuilder(cmd);
		BufferedReader br = ppb.run();
		checkProdigalSearch(br);
		
		// read orfs 
		read(protein);
		
		// count orfs
		int count = 0;
		for (Contig c : genome.contigs())
			count += c.orfs().size();

		if (count == 0)
			// something probably went wrong
			throw new ProdigalSearchException("Error: "
					+ "Prodigal couldn't find any orfs in input sequence!");
		logger.log(Level.INFO, "Found " + count
				+ " orfs with Prodigal within "
				+ genome.contigs().size() + " contigs");
	}

	/**
	 * Determine the mode in which to run Prodigal. "Normal" mode (single) is
	 * appropriate when the length of the sequence is >500 kb and involves
	 * training Prodigal on the input sequence. "Anonymous" mode (meta) involves
	 * applying pre-calculated training files to the provided input sequence and
	 * is the default mode when the length of the longest contig is <500 kb.
	 * 
	 * @return the Prodigal mode to use (single or meta)
	 */
	public String getProdigalMode() {
		int largest = 0;
		for (Contig contig : genome.contigs())
			if (contig.length() > largest)
				largest = contig.length();
		if (largest > 500_000)
			return "single";
		return "meta";
	}
	
	/**
	 * Make sure Prodigal was able to execute successfully.
	 * 
	 * @param br
	 *            buffered input stream reader from the Process which ran
	 *            Prodigal
	 * @throws ProdigalSearchException
	 *             if Prodigal cannot read the input sequence, or the sequence
	 *             contains <20,000 characters
	 * @throws IOException
	 */
	public void checkProdigalSearch(BufferedReader br)
			throws ProdigalSearchException, IOException {
		String line;
		while ((line = br.readLine()) != null) {
			if (line.contains("Sequence read failed")) {
				throw new ProdigalSearchException("Error: "
						+ "Prodigal could not read input sequence!");
			} else if (line.contains("Sequence must be 20000 characters")) {
				throw new ProdigalSearchException("Error: "
						+ "sequence is less than 20,000 characters!");
			}
		}
	}

	public void read(String file) throws IOException, ProdigalSearchException {
		List<Contig> contigs = genome.contigs();
		BufferedReader br = new BufferedReader(new FileReader(file));
		Contig contig = null;
		Orf orf = null;
		StringBuffer sequence = new StringBuffer();
		String line = null;
		while ((line = br.readLine()) != null) {
			if (line.startsWith(">")) {
				// add the last orf, and reset 
				addOrf(orf, sequence, contig);
				contig = null; // important: reset contig here
				orf = null;

				// get name 
				String[] split = line.split("#");
				String name = split[0].trim();
				String orfContigName = name.substring(1, name.lastIndexOf("_"));
				for (Contig c : contigs) {
					String contigName = c.header().split("\\s+")[0];
					if (contigName.equals(orfContigName))
						contig = c;
				}
				
				String orfName = "orf" + name.substring(name.lastIndexOf("_"));
				orf = new Orf(orfName, "");
				
				// set start and end 
				String startRaw = split[1].trim();
				String endRaw = split[2].trim();
				String frameRaw = split[3].trim();
				int start = Integer.parseInt(startRaw) - 1;
				int end = Integer.parseInt(endRaw);
				String frame = frameRaw.contains("-") ? "-" : "+";
				orf.setStart(start);
				orf.setEnd(end);
				orf.setFrame(frame);
				
				// set partial
				String partialRaw = split[split.length - 1].trim();
				int partialIdx = partialRaw.indexOf("partial=");
				String partial = partialRaw.substring(partialIdx+8, partialIdx+10);
				orf.setPartial(partial);
				
				// set prediction mode
				orf.setMode(GenePredictionModes.PRODIGAL);
			} else {
				sequence.append(line.trim());
			}
		}
		// append last orf
		addOrf(orf, sequence, contig);
		
		br.close();
	}

	private void addOrf(Orf orf, StringBuffer sequence, Contig contig) {
		// can we add a new orf?
		if (orf != null && sequence.length() > 0 && contig != null) {
			// ignore the * (stop codon)
			String seq = sequence.toString();
			if (seq.endsWith("*"))
				seq = seq.substring(0, seq.length() - 1);

			// set amino acid sequence and add to contig
			orf.setAASequence(seq);
			contig.orfs().add(orf);

			// reset
			contig = null;
			orf = null;
			sequence.delete(0, sequence.length());
		}
	}

}
