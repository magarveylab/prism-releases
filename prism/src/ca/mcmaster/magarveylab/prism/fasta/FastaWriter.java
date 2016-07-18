package ca.mcmaster.magarveylab.prism.fasta;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.prism.cluster.analysis.DomainAnalyzer;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Contig;
import ca.mcmaster.magarveylab.prism.data.Orf;

/**
 * Writes FASTA files.
 * @author skinnider
 *
 */
public class FastaWriter {

	/**
	 * Print a header and its associated sequence to a file.
	 * @param header		FASTA header to print
	 * @param sequence		sequence to print
	 * @param path			path of the FASTA file to print to
	 * @throws IOException
	 */
	public static void printToFasta(String header, String sequence, String path) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(path));
		bw.append(">" + header + "\n");
		bw.append(sequence + "\n");
		bw.close();
	}

	/**
	 * Print the amino acid sequences of a list of domains to a multi-FASTA file. 
	 * @param domains	the domains to print
	 * @param fasta		the location of the FASTA file
	 * @throws IOException
	 */
	public static void printDomainsToFasta(List<Domain> domains, String fasta) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(fasta));
		for (Domain d : domains) {
			bw.append(d.printForFASTA());
		}
		bw.close();
	}

	/**
	 * Print the amino acid sequence of a cluster's (scaffold) orfs to a multi-FASTA file.
	 * @param cluster	the cluster to print
	 * @param fasta		the location of the FASTA file
	 * @throws IOException
	 */
	public static void printOrfsToFasta(List<Orf> orfs, String fasta) throws IOException {
		File fastaFile = new File(fasta);
		if (!fastaFile.exists())
			fastaFile.createNewFile();
		BufferedWriter bw = new BufferedWriter(new FileWriter(fasta));
		for (Orf orf : orfs) {
			String orfFasta = FastaUtil.formatOrfForFasta(orf);
			bw.append(orfFasta);
		}
		bw.close(); 
	}

	/**
	 * Print all C starter domains to a file.
	 * @param genome	current genome
	 * @param file		location of file
	 * @throws IOException 
	 */
	public static void printCStarter(Contig contig, String file) throws IOException {
		List<Domain> cstarter = new ArrayList<Domain>();
		for (Orf orf : contig.orfs()) 
			for (Domain domain : orf.domains(ThiotemplatedDomains.CONDENSATION))
				if (DomainAnalyzer.isCStarter(domain))
					cstarter.add(domain);		
		printDomainsToFasta(cstarter, file);
	}
	
}
