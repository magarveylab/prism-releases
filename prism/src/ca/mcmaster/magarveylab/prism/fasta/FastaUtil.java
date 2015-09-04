package ca.mcmaster.magarveylab.prism.fasta;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Contig;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.data.Organism;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Utilities class for FASTA file operations.
 * @author skinnider
 *
 */
public class FastaUtil {

	/**
	 * Generate the amino acid sequence of a given orf in FASTA format.
	 * @param orf	orf to generate
	 * @return		sequence in FASTA format, without trailing newline
	 */
	public static String formatOrfForFasta(Orf orf) {
		String sequence = formatSequenceForFasta(orf.sequence());
		return ">" + orf.name() + "\n" + sequence;
	}
	
	/**
	 * Get the amino acid sequence of a domain without prior knowledge of what orf it is on.
	 * @param domain	domain to get sequence for 
	 * @param cluster	parent cluster
	 * @return			the amino acid sequence of the query domain
	 */
	public static String getDomainSequence(Domain domain, Cluster cluster) {
		String sequence = null;
		for (Orf orf : cluster.orfs())
			if (orf.contains(domain))
				sequence = orf.sequence().substring(domain.start(), domain.end());
		return sequence;
	}
	
	/**
	 * Format any sequence for FASTA by inserting line breaks after every 50th character.
	 * @param sequence	sequence to format
	 * @return			sequence with line breaks
	 */
	public static String formatSequenceForFasta(String sequence) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < sequence.length(); i += 50) {
			String temp = (i + 50 > sequence.length()) ? sequence.substring(i) : sequence.substring(i, i + 50);
			sb.append(temp + "\n");
		}
		return sb.toString();
	}
	
	public static Organism parseOrganism(String header) {
		Organism o = new Organism(header);

		String[] data = new String[5];
		//String[] split = header.split("\\|");
		
		String split2 = header.substring(0, header.indexOf(" "));
		String other = header.substring(header.indexOf(" ") + 1);
		String[] split3 = split2.split("\\|");
		
		String accession = split3[3];	//split[3]		
		String pattern = "([^\\s]+)\\s?([^\\s]+)?([^,]+),?(.+)?";
		Pattern p = Pattern.compile(pattern);
		Matcher m = p.matcher(other); //split[4]
		
		if (m.find()) {
			data[0] = accession;
			if (m.group(1) != null)
				data[1] = m.group(1);
			if (m.group(1) != null)
				data[2] = m.group(2);
			if (m.group(1) != null)
				data[3] = m.group(3);
			if (m.group(1) != null)
				data[4] = m.group(4);
			
			String contig = new String();
			if (data[3].split(" ").length > 1  && data.length > 3) {
				String[] strain = data[3].split(" ");
				Pattern p2 = Pattern.compile("(?i)((contig|node|scaffold|ctg).+)$");
				Matcher m2 = p2.matcher(strain[strain.length -1]); //split[4]
				if (m2.find()) {
					contig = m2.group(1);
					int cut = strain[strain.length -1].indexOf(contig);
					if (cut > data[3].length()) {
						cut = data[3].length();
					}
					data[3] = data[3].substring(0, cut);
					StringBuilder newstrain = new StringBuilder();
					for (int i = 0; i < strain.length-1; i++) {
						newstrain.append(strain[i]);
						newstrain.append(" ");
					}
					data[3] = newstrain.toString();
				}
			}
			return new Organism(header, data[0], data[1], data[2], data[3], data[4], contig); 
		} else {
			return o;
		}
	}
	
	/**
	 * Initialize file locations for a contig.
	 * @param contig	contig in question
	 */
	public static void setFiles(Contig contig, Session session) {
		contig.setFile("orfs", session.dir() + "contig_" + contig.index() + "_orfs.fa");
		contig.setFile("adenylation", session.dir() + "contig_" + contig.index() + "_adenylation.fa");
		contig.setFile("acyl_adenylating", session.dir() + "contig_" + contig.index() + "_acyl_adenylating.fa");
		contig.setFile("acyltransferase", session.dir() + "contig_" + contig.index() + "_acyltransferase.fa");
		contig.setFile("priming_AT", session.dir() + "contig_" + contig.index() + "_priming_at.fa");
		contig.setFile("condensation", session.dir() + "contig_" + contig.index() + "_condensation.fa");
		contig.setFile("ketosynthase", session.dir() + "contig_" + contig.index() + "_ketosynthase.fa");
		contig.setFile("CLF", session.dir() + "contig_" + contig.index() + "_clf.fa");
		contig.setFile("cstarter", session.dir() + "contig_" + contig.index() + "_cstarter.fa");
		contig.setFile("chlorination", session.dir() + "contig_" + contig.index() + "_chlorination.fa");
		contig.setFile("glycosyltransferase", session.dir() + "contig_" + contig.index() + "_glycosyltransferase.fa");
	}
	
}