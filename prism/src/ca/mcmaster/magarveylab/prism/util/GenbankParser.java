package ca.mcmaster.magarveylab.prism.util;

import ca.mcmaster.magarveylab.prism.data.Genome;
import ca.mcmaster.magarveylab.prism.data.Organism;
import ca.mcmaster.magarveylab.prism.data.Contig;
import ca.mcmaster.magarveylab.prism.data.RnaSequence;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava3.core.sequence.DNASequence;

/**
 * @author prees
 * Will parse a genbank file and generate a Genome object.
 */
public class GenbankParser {
		

	/**
	 * @param genome A Genome object to add information to
	 * @param file A string representing a path to the genbank file to parse
	 * @param session The session of the prism run.
	 * @return genome A genome object with the additional informated added.
	 */
	public static Genome readGenbankFile(Genome genome, String file, String root, String output) throws IOException, InterruptedException, Exception {	
		final Logger logger = Logger.getLogger(GenbankParser.class.getName());
		logger.log(Level.INFO, "Parsing Genbank file"); 

		Process processSeq = runBioPerl(file, "read_genbank_seq.pl", root, output);
		InputStream isSeq = processSeq.getInputStream();
		BufferedReader brSeq = new BufferedReader(new InputStreamReader(isSeq));

		StringBuilder seqSB = new StringBuilder();
		char[] chars = new char[8192];

		while (brSeq.read(chars) > 0) {
		    seqSB.append(String.valueOf(chars));
		//	chars = new char[chars.length];
		}

		String sequence = seqSB.toString();
		
		if (sequence.startsWith("NOSEQ")) {
			logger.log(Level.SEVERE, "No Sequence Found in Genbank file."); 
		}else {
			sequence = sequence.replaceAll("[^ACTGN]", "N");
			DNASequence seq = new DNASequence(sequence);

			Process process = runBioPerl(file, "read_genbank.pl", root, output);
			InputStream is = new BufferedInputStream(process.getInputStream());
			BufferedReader in = new BufferedReader(new InputStreamReader(is));
	
			String description = in.readLine();
			String accession = in.readLine();
			String classificationOne = in.readLine();
			List<String> classification = new ArrayList<String>(Arrays.asList(classificationOne.split("\\t")));
			
			String genus = in.readLine();
			String species = in.readLine();
			String strain = in.readLine();

			List<RnaSequence> allRnaSeq = new ArrayList<RnaSequence>();
			String line;
			boolean sixteen = false;

			while ((line = in.readLine()) != null) {

				if ("16S".equals(line)) {

					Integer start = Integer.valueOf(in.readLine());
					Integer end = Integer.valueOf(in.readLine());
					String strand = in.readLine();
					HashMap<String, String> data = new HashMap<String, String>();
					while ((line = in.readLine()) !=  "###") {
						if (line.equals("###")) {
							break;
						}
						String[] split = line.split("\\t");
						data.put(split[0], split[1]);	
					}
					
					String finalSequence;				
					//if strand is -1 means it needs the reverse complement from the original sequence
					if (strand.equals("-1")) {
						DNASequence sub = (DNASequence) seq.getSubSequence(start, end).getViewedSequence();
						finalSequence = sub.getReverseComplement().getSequenceAsString();
					}else {
						finalSequence = seq.getSequenceAsString();
					}
					

					String parsed = finalSequence.substring(start, end + 1);
					data.put("sequence", parsed);
					
					RnaSequence rna = new RnaSequence(start, end, data, strand);
					allRnaSeq.add(rna);
					sixteen = true;
				}
			}


			in.close();
		
			List<Contig> contigs = new ArrayList<Contig>();
			contigs.add(new Contig(description, sequence));


			Organism org = new Organism(description, accession, genus, species, strain, classification, "genbank");
			genome.setContigs(contigs);
			genome.setOrganism(org);
			genome.setFile(new File(file));
			
			if (sixteen) 
				genome.setRibosomalSequences(allRnaSeq);
			
		}

		return genome;
	}	
	
	/**
	 * @param file A string representing a path to the genbank file to parse
	 * @param script The name of the perl script to execute.
	 * @param session The session of the prism run.
	 * @return process The process that was executed.
	 */
	private static Process runBioPerl(String file, String script, String root, String output) throws IOException, InterruptedException {
		
		List<String> commandStringBuilder = new ArrayList<String>();
		commandStringBuilder.add("perl");
		commandStringBuilder.add(root + script);
		commandStringBuilder.add(file);
		System.out.println(commandStringBuilder.toString());

		ProcessBuilder pb = new ProcessBuilder(commandStringBuilder);
		pb.redirectError(new File(output + "errors.txt"));
		Process process = pb.start();
	//	process.waitFor();
		
		return process;
	}
}
