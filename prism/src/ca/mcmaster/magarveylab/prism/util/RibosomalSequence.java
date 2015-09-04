package ca.mcmaster.magarveylab.prism.util;

import ca.mcmaster.magarveylab.prism.data.Contig;
import ca.mcmaster.magarveylab.prism.data.RnaSequence;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.biojava3.core.sequence.DNASequence;

/**
 * @author prees
 * Implements the barrnap tool to identify 16s sequences and parse the results and return a objects representing those 16s sequences.
 */

public class RibosomalSequence {		
	File path;
	Process process;
	List<Contig> contigs;

	/**
	 * @param file A filepath as a String containing sequence information.
	 * @param List<Contig> A list of contigs that represent the same sequences contained in the file.
	 *
	 */
	public RibosomalSequence(List<Contig> contigs, String file) {
		this.path = new File(file);
		this.contigs = contigs;
	}
	

	public RibosomalSequence(File file) {
		this.path = file;
	}
	
	/**
	 * Execute the process to look for 16s sequences.
	 *
	 */
	public boolean runProcess() throws IOException, InterruptedException {
		ProcessBuilder pb = new ProcessBuilder(buildCommand());
		pb.redirectError(new File("errors.txt"));
		this.process = pb.start();
		this.process.waitFor();
		return true;
	}
	
	
	/**
	 * Execute the process to look for 16s sequences.
	 * @return A list of strings that are the command to be executed.
	 */
	private List<String> buildCommand() {
		List<String> commandStringBuilder = new ArrayList<String>();
		commandStringBuilder.add("barrnap");
		commandStringBuilder.add("--quiet");
		commandStringBuilder.add(this.path.toString());
		return commandStringBuilder;
	}
	
	
	/**
	 * Parse the results that were executed and return the results
	 * @return A list of RnaSequence objects containing any information the process found.
	 */
	public List<RnaSequence> parseRawOutput() throws IOException, Exception {
		BufferedReader in = new BufferedReader(new InputStreamReader(this.process.getInputStream()));
		String line = new String();
		ArrayList<RnaSequence> parsedSequences = new ArrayList<RnaSequence>();
		while((line = in.readLine()) != null) {
			if (!(line.startsWith("##"))) {
				String parts[] = line.split("\t");
				String[] attributes = parts[8].split(";");
				String name = new String();
				String product = new String();
				for (String attribute : attributes) {
					if (attribute.startsWith("Name=")) {
						name = attribute.substring(5);
					}else if (attribute.startsWith("product=")) {
						product = attribute.substring(8);
					}
				}
				
				if (name.equals("16S_rRNA")) {		
					StringBuilder full = new StringBuilder();
					for (Contig key :contigs) {
						full.append(key.sequence());
					}
					
					DNASequence seq = new DNASequence(full.toString());
					
					//if strand is -1 means it needs the reverse complement from the original sequence
					String finalSequence;
					if (parts[6].equals("-")) {
						DNASequence sub = (DNASequence) seq.getSubSequence(Integer.valueOf(parts[3]), Integer.valueOf(parts[4])).getViewedSequence();
						finalSequence = sub.getReverseComplement().getSequenceAsString();
					}else {
						finalSequence = seq.getSequenceAsString();
					}
					
	
					String parsed = finalSequence.substring(Integer.valueOf(parts[3]), Integer.valueOf(parts[4]) + 1);
					
				
					RnaSequence rnaSeq = new RnaSequence(Integer.valueOf(parts[3]), Integer.valueOf(parts[4]), parts[6], name, product);
					HashMap<String, String> data = new HashMap<String, String>();
					data.put("sequence", parsed);
					rnaSeq.setData(data);
					parsedSequences.add(rnaSeq);
				}
				
			}
		}
		return parsedSequences;
	}
	
}
