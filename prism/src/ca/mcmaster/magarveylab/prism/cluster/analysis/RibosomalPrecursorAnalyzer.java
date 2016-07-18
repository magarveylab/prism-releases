package ca.mcmaster.magarveylab.prism.cluster.analysis;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import ca.mcmaster.magarveylab.enums.ClusterFamilies;
import ca.mcmaster.magarveylab.enums.DomainFamilies;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.analysis.type.RibosomalClusterTypeAnalyzer;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Contig;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.orfs.GenePredictionModes;
import ca.mcmaster.magarveylab.prism.util.Sorter;

/**
 * Analyzes putative precursor peptides for ribosomally synthesized and
 * post-translationally modified peptides (RiPPs).
 * 
 * @author skinnider
 *
 */
public class RibosomalPrecursorAnalyzer {

	/**
	 * Get all ribosomal peptide precursor domains within a cluster.
	 * 
	 * @param cluster
	 *            parent cluster
	 * @return all ribosomal peptide precursors
	 */
	public static List<Domain> getPrecursors(Cluster cluster) {
		List<Domain> precursors = new ArrayList<Domain>();
		List<DomainType> precursorTypes = getPrecursorTypes();
		for (Domain domain : cluster.domains()) {
			DomainType domainType = domain.type();
			for (DomainType precursorType : precursorTypes)
				if (domainType == precursorType) {
					precursors.add(domain);
					break;
				}
		}
		return precursors;
	}

	/**
	 * Get all lantipeptide precursor domains within a cluster.
	 * 
	 * @param cluster
	 *            parent cluster
	 * @return all lantipeptide precursors
	 */
	public static List<Domain> getLantipeptidePrecursors(Cluster cluster) {
		List<Domain> precursors = new ArrayList<Domain>();
		List<DomainType> lantipeptidePrecursors = getLantipeptidePrecursorTypes();
		for (Domain domain : cluster.domains()) {
			DomainType domainType = domain.type();
			for (DomainType precursorType : lantipeptidePrecursors)
				if (domainType == precursorType) {
					precursors.add(domain);
					break;
				}
		}
		return precursors;
	}

	/**
	 * Determine whether this cluster contains a lantipeptide precursor domain.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if the cluster contains a lantipeptide domain
	 */
	public static boolean containsLantipeptidePrecursor(Cluster cluster) {
		return getLantipeptidePrecursors(cluster).size() > 0;
	}

	/**
	 * Detect putative domains in a ribosomal natural product cluster. These
	 * domains are identified by heuristics rather than hidden Markov models
	 * when the presence of other ribosomal natural product biosynthesis domains
	 * suggests the presence of a cluster.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @param contig
	 *            the parent contig
	 */
	public static void detectPutativePrecursors(Cluster cluster, Contig contig) {
		// first, try prodigal orfs
		detectPutativePrecursors(cluster, contig, GenePredictionModes.PRODIGAL);
		if (getPrecursors(cluster).size() == 0)
			detectPutativePrecursors(cluster, contig,
					GenePredictionModes.ALL_ORFS);
	}

	/**
	 * Detect putative domains in a ribosomal natural product cluster using
	 * heuristics rather than hidden Markov models. This method prioritizes open
	 * reading frames identified with a specific method: e.g., with Prodigal or
	 * by scanning the sequence for all possible open reading frames.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @param contig
	 *            the parent contig
	 * @param mode
	 *            the method used to identify open reading frames
	 */
	public static void detectPutativePrecursors(Cluster cluster, Contig contig,
			GenePredictionModes mode) {
		for (Orf orf : contig.getAllOrfs(cluster, mode, 5_000)) {
			if (RibosomalClusterTypeAnalyzer.isLassoPeptideCluster(cluster)
					&& !cluster.contains(RibosomalDomains.Lasso_precursors)) {
				List<int[]> precursors = getPutativeLassoPeptidePrecursors(orf);
				for (int[] precursor : precursors) {
					System.out.println("Detected putative lasso peptide "
							+ "precursor in cluster " + cluster.index());
					Domain putative = new Domain(precursor[0] + 2,
							orf.length() / 3, -1.0d,
							"Putative_lasso_peptide_precursor_" + orf.name());
					putative.setFamily(DomainFamilies.RIBOSOMAL);
					putative.setType(RibosomalDomains.Putative_lasso_peptide);
					orf.add(putative);
					if (!cluster.orfs().contains(orf))
						cluster.orfs().add(orf);
				}
			} else if (RibosomalClusterTypeAnalyzer
					.isLinearAzoleCluster(cluster)
					&& getPrecursors(cluster).size() == 0
					&& isPutativeLAPPrecursor(orf)) {
				System.out.println("Detected putative LAP "
						+ "precursor in cluster " + cluster.index());
				Domain putative = new Domain(orf.length() / 6,
						orf.length() / 3 - 1, -1.0d, "Putative_LAP_precursor_"
								+ orf.name());
				putative.setFamily(DomainFamilies.RIBOSOMAL);
				putative.setType(RibosomalDomains.Putative_LAP);
				orf.add(putative);
				if (!cluster.orfs().contains(orf))
					cluster.orfs().add(orf);
			} else if (RibosomalClusterTypeAnalyzer
					.isLantipeptideCluster(cluster)
					&& getPrecursors(cluster).size() == 0
					&& isPutativeLantipeptidePrecursor(orf)) {
				System.out.println("Detected putative lantipeptide "
						+ "precursor in cluster " + cluster.index());
				Domain putative = new Domain(0, orf.length() / 3 - 1, -1.0d,
						"Putative_lantipeptide_precursor_" + orf.name());
				putative.setFamily(DomainFamilies.RIBOSOMAL);
				putative.setType(RibosomalDomains.Putative_lantipeptide);
				orf.add(putative);
				if (!cluster.orfs().contains(orf))
					cluster.orfs().add(orf);
			} else if (RibosomalClusterTypeAnalyzer
					.isThiopeptideCluster(cluster)
					&& getPrecursors(cluster).size() == 0
					&& isPutativeThiopeptidePrecursor(orf)) {
				System.out.println("Detected putative thiopeptide "
						+ "precursor in cluster " + cluster.index());
				Domain putative = new Domain(0, orf.length() / 3 - 1, -1.0d,
						"Putative_thiopeptide_precursor_" + orf.name());
				putative.setFamily(DomainFamilies.RIBOSOMAL);
				putative.setType(RibosomalDomains.Putative_thiopeptide);
				orf.add(putative);
				if (!cluster.orfs().contains(orf))
					cluster.orfs().add(orf);
			}
		}
		Sorter.sortOrfs(cluster.orfs());
	}

	/**
	 * Implement a simple heuristic strategy for identifying lantipeptide
	 * precursors when a combination of precursor-derived HMMs and TIGRFAM HMMs
	 * (for Nif11/NHLP-family leaders) fail to identify any precursor peptide.
	 * Open reading frames between 30-80 aa with at least 2 cysteines are
	 * considered potential lantipeptide precursors.
	 * 
	 * @param orf
	 *            open reading frame to analyze
	 * @return true if this orf meets heuristic criteria for lantipeptide
	 *         precursors
	 */
	public static boolean isPutativeLantipeptidePrecursor(Orf orf) {
		boolean flag = false;
		if (orf.length() > 240 || orf.length() < 90)
			return flag;

		int cys = 0;
		String seq = orf.sequence().toLowerCase();
		for (int i = 0; i < seq.length(); i++)
			if (seq.charAt(i) == 'c')
				cys++;
		return cys >= 2;
	}

	/**
	 * Implement a simple heuristic strategy for identifying thiopeptide
	 * precursors when a precursor-derived HMM fails to identify a precursor
	 * peptide. Open reading frames of length 30-90 aa are considered potential
	 * thiopeptide precursors when they have a stretch of 12 to 20 amino acids
	 * composed of at least 40% serines, threonines, and cysteines, and
	 * containing at least two serines.
	 * 
	 * @param orf
	 *            open reading frame to analyze
	 * @return true if this orf meets heuristic criteria for thiopeptide
	 *         precursors
	 */
	public static boolean isPutativeThiopeptidePrecursor(Orf orf) {
		boolean flag = false;
		if (orf.length() > 270 || orf.length() < 90)
			return flag;

		String sequence = orf.sequence();
		outerLoop: for (int length = 12; length <= 20; length++) {
			for (int i = 0; i < sequence.length() - length; i++) {
				String subsequence = sequence.substring(i, i + length);
				int s = 0, t = 0, c = 0;
				for (char ch : subsequence.toLowerCase().toCharArray()) {
					if (ch == 's')
						s++;
					if (ch == 't')
						t++;
					if (ch == 'c')
						c++;
				}
				int sum = s + t + c;
				double fraction = 1.0d * sum / length;
				if (fraction >= 0.4 && s >= 2) {
					flag = true;
					break outerLoop;
				}
			}
		}

		return flag;
	}
	
	/**
	 * Implement a simple heuristic strategy for identifying linear
	 * azol(in)e-containing peptide precursors. Open reading frames between
	 * 30-70 amino acids which contain a sequence of 7 amino acids of which all
	 * 7 are serine, threonine, or cysteine, or a sequence of 8 amino acids of
	 * which at least 7 are serine, threonine, or cysteine, are considered
	 * potential LAP precursors.
	 * 
	 * @param orf
	 *            open reading frame to analyze
	 * @return true if this orf meets heuristic criteria for linear
	 *         azole-containing peptide precursors
	 */
	
	public static boolean isPutativeLAPPrecursor(Orf orf) {
		boolean flag = false;
		if (orf.length() > 210 || orf.length() < 90)
			return flag;

		// 7 of 7 C/S/T
		String sequence = orf.sequence();
		for (int i = 0; i < sequence.length() - 7; i++) {
			String substring = sequence.substring(i, i + 7);
			int cst = 0;
			for (int j = 0; j < substring.length(); j++) {
				char c = substring.toLowerCase().charAt(j);
				if (c == 'c' || c == 's' || c == 't')
					cst++;
			}
			if (cst == 7)
				flag = true;
		}
		
		// 7+ of 8 C/S/T 
		for (int i = 0; i < sequence.length() - 8; i++) {
			String substring = sequence.substring(i, i + 8);
			int cst = 0;
			for (int j = 0; j < substring.length(); j++) {
				char c = substring.toLowerCase().charAt(j);
				if (c == 'c' || c == 's' || c == 't')
					cst++;
			}
			if (cst >= 7)
				flag = true;
		}

		return flag;
	}

	/**
	 * Implement a simple heuristic strategy for identifying lasso peptide
	 * precursors described by Maksimov, Pelczer, and Link (PNAS 2012). When a
	 * lasso peptide precursor cannot be identified by a hidden Markov model, a
	 * putative lasso peptide precursor may be identified in the cluster as an
	 * open reading frame < 80 amino acids in length (Maksimov et al. suggest
	 * 73), with a TX(S/G/C)X[6-8](D/E) motif 0-50 residues from the beginning
	 * and >5 residues from the end (Maksimov et al. suggest X[6-10] at a
	 * distance of 5-43 amino acids from the beginning and >5 residues from the
	 * end).
	 * 
	 * @param orf
	 *            open reading frame to analyze
	 * @return a list of integer arrays corresponding to the positions of the T,
	 *         G, and D/E residues which meets heuristic criteria for lasso
	 *         peptide precursors
	 */
	public static List<int[]> getPutativeLassoPeptidePrecursors(Orf orf) {
		List<int[]> precursors = new ArrayList<int[]>();
		String sequence = orf.sequence().toLowerCase();

		// PRISM often finds excessively long orfs, so get distance from
		// last *
		int asterisk = sequence.lastIndexOf('*');
		if (asterisk == -1)
			asterisk = 0;
		// cannot be >80 aa long
		String subsequence = sequence.substring(asterisk);
		if (subsequence.length() > 80)
			return precursors;

		// 2 ~ TXG, 10 ~ GX[6-8]D/E, 5 ~ end
		int max = sequence.length() - 2 - 10 - 5;
		for (int i = asterisk + 5; i < max && i <= (asterisk + 50); i++) {
			char c1 = sequence.charAt(i);
			if (c1 == 't') {
				// look for TxG
				int j = i + 2;
				char c2 = sequence.charAt(j);
				if (c2 == 'g' || c2 == 'c') {
					// look for Gx[6-8]D/E to form 7-9 membered rings
					for (int k = j + 6; k <= j + 8; k++) {
						char c3 = sequence.charAt(k);
						if (c3 == 'd' || c3 == 'e') {
							int[] precursor = new int[] { i, j, k };
							precursors.add(precursor);
							break;
						}
					}
				}
			}
		}

		return precursors;
	}

	/**
	 * Remove ribosomal precursors from non-ribosomal clusters.
	 * 
	 * @param cluster
	 *            cluster in question
	 */
	public static void removePrecursors(Cluster cluster) {
		List<DomainType> precursors = getPrecursorTypes();
		if (!cluster.families().contains(ClusterFamilies.RIBOSOMAL)) {
			for (Orf orf : cluster.orfs()) {
				Iterator<Domain> itr = orf.domains().iterator();
				while (itr.hasNext()) {
					Domain next = itr.next();
					if (Arrays.asList(precursors).indexOf(next.type()) != -1)
						itr.remove();
				}
			}
		}
	}

	/**
	 * Get all ribosomal peptide precursor domain types.
	 * 
	 * @return all ribosomal precursor types
	 */
	public static List<DomainType> getPrecursorTypes() {
		List<DomainType> types = new ArrayList<DomainType>();
		for (RibosomalDomains rd : RibosomalDomains.values())
			if (rd.isPrecursor())
				types.add(rd);
		return types;
	}

	/**
	 * Get all lantipeptide precursor domain types.
	 * 
	 * @return all lantipeptide precursor types
	 */
	public static List<DomainType> getLantipeptidePrecursorTypes() {
		List<DomainType> types = new ArrayList<DomainType>();
		for (RibosomalDomains rd : RibosomalDomains.values())
			if (rd.isLantipeptidePrecursor())
				types.add(rd);
		return types;
	}

}
