package ca.mcmaster.magarveylab.prism.fasta;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.prism.data.Contig;
import ca.mcmaster.magarveylab.prism.data.Genome;

/**
 * Reads FASTA files.
 * 
 * @author skinnider
 *
 */
public class FastaReader {
	
	/**
	 * Read all contigs from a FASTA/multi-FASTA format sequence file. 
	 * @throws IOException
	 * @throws BadSequenceException
	 */
	public static void readFastaFile(Genome genome) throws IOException {
		File file = genome.file();
		List<Contig> contigs = FastaReader.readMultiFasta(file, 500);
		genome.setContigs(contigs);
	}

	/**
	 * Read a generic multi-FASTA file to a list of Fasta objects.
	 * 
	 * @param file
	 *            the location of the file to read
	 * @param minSize
	 *            the minimum sequence size to be read
	 * @return a list of FastaItem objects
	 * @throws IOException
	 */
	public static List<Contig> readMultiFasta(File file, int minSize)
			throws IOException {
		List<Contig> contigs = new ArrayList<Contig>();
		if (file.exists()) {
			int counter = 1;
			String line;
			StringBuffer name = new StringBuffer();
			StringBuffer sequence = new StringBuffer();
			InputStream in = new FileInputStream(file);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			while ((line = br.readLine()) != null) {
				if (line.startsWith(">")) {
					if (name.length() > 0) { // can't be first 
						if (sequence.length() > minSize) {
							Contig item = new Contig(name.toString(),
									sequence.toString());
							item.setIndex(counter);
							contigs.add(item);
							name.delete(0, name.length());
							sequence.delete(0, sequence.length());
							counter++;
						} else {
							name.delete(0, name.length());
							sequence.delete(0, sequence.length());
						}
					}
					String header = line.split(">")[1].trim();
			//		header = header.replace(" ", "_");
					name.append(header);
				} else {
					sequence.append(line.trim());
				}
			}
			// add the last orf (when there is no more > to mark the end)
			Contig contig = new Contig(name.toString(), sequence.toString());
			contig.setIndex(counter);
			contigs.add(contig);
			br.close();
		}
		return contigs;
	}

}
