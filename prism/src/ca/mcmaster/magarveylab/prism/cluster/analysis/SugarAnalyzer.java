package ca.mcmaster.magarveylab.prism.cluster.analysis;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.paukov.combinatorics.Factory;
import org.paukov.combinatorics.Generator;
import org.paukov.combinatorics.ICombinatoricsVector;

import ca.mcmaster.magarveylab.enums.DeoxySugars;
import ca.mcmaster.magarveylab.enums.DomainFamilies;
import ca.mcmaster.magarveylab.enums.HexoseSugars;
import ca.mcmaster.magarveylab.enums.domains.DeoxySugarDomains;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.prism.blast.BlastSearchResult;
import ca.mcmaster.magarveylab.prism.combinatorialization.Combinations;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.CombinatorialData;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.sugar.DeoxySugar;
import ca.mcmaster.magarveylab.prism.data.sugar.HexoseSugar;
import ca.mcmaster.magarveylab.prism.data.sugar.Sugar;
import ca.mcmaster.magarveylab.prism.util.Sorter;

/**
 * Analyzes sugars for combinatorial library generation.
 * 
 * @author skinnider
 *
 */
public class SugarAnalyzer {

	/**
	 * Get all possible combinations of sugars based on analysis of a cluster's
	 * deoxysugar biosynthesis genes and glycosyltransferases.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return a list of list of sugars, where each list represents a possible
	 *         sugar combination
	 */
	public static List<List<Sugar>> getSugars(Cluster cluster) {
		List<List<Sugar>> sugars = new ArrayList<List<Sugar>>();

		// get deoxysugar size from GTrs 
		List<Domain> gtrs = cluster.domains(TailoringDomains.GLYCOSYLTRANSFERASE);
		int k = gtrs.size();
		for (Domain gtr : gtrs)
			if (isHexoseGlycosyltransferase(gtr))
				k--;
		k += cluster.domains(TailoringDomains.C_GLYCOSYLTRANSFERASE).size();

		// get all combinations of deoxysugars, size n
		if (k > 0) {
			System.out.println("Getting " + k + " combinations with repetitions");
			List<List<DeoxySugars>> combinations = getSugarCombinations(k);
			System.out.println("[SugarAnalyzer] Analyzing " + combinations.size() + " possible sugar combinations");

			// for each combination, calculate the # of sugar genes "left over" in the cluster
			int minimumUnusedFromCluster = 999;
			int minimumUnusedFromSugars = 999;
			for (List<DeoxySugars> combination : combinations) {
				List<Domain> domains = cluster.domains(DomainFamilies.SUGAR);

				// get all sugar genes 
				List<DeoxySugarDomains> genes = new ArrayList<DeoxySugarDomains>();
				for (DeoxySugars sugar : combination)
					for (DeoxySugarDomains gene : sugar.genes())
						if (genes.indexOf(gene) == -1)
							genes.add(gene);

				Iterator<Domain> itr1 = domains.iterator();
				while (itr1.hasNext()) {
					Domain domain = itr1.next();
					if (contains(domain, genes))
						itr1.remove();
				}

				Iterator<DeoxySugarDomains> itr2 = genes.iterator();
				while (itr2.hasNext()) {
					DeoxySugarDomains gene = itr2.next();
					if (contains(gene, cluster))
						itr2.remove();
				}

				int unusedFromCluster = domains.size();
				int unusedFromSugars = genes.size();

				if (unusedFromCluster <= minimumUnusedFromCluster && unusedFromSugars <= minimumUnusedFromSugars) {
					if (unusedFromCluster < minimumUnusedFromCluster || unusedFromSugars < minimumUnusedFromSugars)
						sugars.clear();
					minimumUnusedFromCluster = unusedFromCluster;
					minimumUnusedFromSugars = unusedFromSugars;
					Sorter.sortSugarsByName(combination);
					sugars.add(convertToSugars(combination));
				}
			}
			// print to console 
			System.out.println("[SugarAnalyzer] Minimized unused from cluster: " + minimumUnusedFromCluster);
			System.out.println("[SugarAnalyzer] Minimized unused from sugar: " + minimumUnusedFromSugars);
		}

		// add hexoses 
		List<Sugar> hexoses = getHexoseSugars(cluster);
		if (sugars.size() > 0) {
			for (Sugar hexose : hexoses)
				for (List<Sugar> list : sugars)
					list.add(hexose);
		} else if (hexoses.size() > 0) {
			sugars.add(new ArrayList<Sugar>());
			sugars.get(0).addAll(hexoses);
		}
		
		// set # sugar combinations
		CombinatorialData cd = cluster.combinatorialData();
		cd.setNumSugars(sugars.size());
		
		// if >100 combinations, only consider the first 100
		if (sugars.size() > 100)
			sugars = sugars.subList(0, 100);

		System.out.println("[SugarAnalyzer] Found " + sugars.size() + " combinations of " + k + " sugars");
		return sugars;
	}

	/**
	 * Convert a list of deoxy sugar types to a list of sugars. 
	 * @param combination	a combination of deoxy sugar types
	 * @return				a list of deoxy sugar objects 
	 */
	public static List<Sugar> convertToSugars(List<DeoxySugars> combination) {
		List<Sugar> sugars = new ArrayList<Sugar>();
		for (DeoxySugars deoxysugar : combination) {
			DeoxySugar sugar = new DeoxySugar(deoxysugar);
			sugars.add(sugar);
		}
		return sugars;
	}

	/**
	 * Get all combinations, with replacement/repetition, of all possible deoxy sugars. 
	 * @param k		size of each combination
	 * @return		all combinations with replacement 
	 */
	public static List<List<DeoxySugars>> getSugarCombinations(int k) {
		List<List<DeoxySugars>> combinations = new ArrayList<List<DeoxySugars>>();
		if (k < 4) {//evaluate exhaustively
			ICombinatoricsVector<DeoxySugars> initialVector = Factory.createVector(DeoxySugars.values());
			Generator<DeoxySugars> generator = Factory.createMultiCombinationGenerator(initialVector, k);
			for (ICombinatoricsVector<DeoxySugars> vector : generator) {
				List<DeoxySugars> combination = vector.getVector();
				combinations.add(combination);
			}
		} else {
			int n = DeoxySugars.values().length;
			List<int[]> sample = Combinations.sampleCombinations(n, k, 1_000);
			for (int[] combinatoric : sample) {
				List<DeoxySugars> combination = new ArrayList<DeoxySugars>();
				for (int index : combinatoric)
					combination.add(DeoxySugars.values()[index]);
				combinations.add(combination);
			}
		}
		return combinations;
	}
	
	/**
	 * Get all hexose sugars identified within this cluster by analyzing the cluster's glycosyltransferase domains. 
	 * @param cluster	query cluster
	 * @return			all hexose sugars in the cluster
	 */
	public static List<Sugar> getHexoseSugars(Cluster cluster) {
		List<Sugar> sugars = new ArrayList<Sugar>();
		List<Domain> gtrs = cluster.domains(TailoringDomains.GLYCOSYLTRANSFERASE);
		for (Domain gtr : gtrs) {
			if (isHexoseGlycosyltransferase(gtr)) {
				HexoseSugars type = getHexoseSugarType(gtr);
				HexoseSugar sugar = new HexoseSugar(type);
				sugars.add(sugar);
			}
		}
		return sugars;
	}

	/**
	 * Determine whether a cluster contains a deoxy sugar biosynthesis gene of a particular type
	 * @param gene		query deoxysugar gene
	 * @param cluster	cluster to analyze
	 * @return			true of the cluster contains a gene of the query type
	 */
	public static boolean contains(DeoxySugarDomains gene, Cluster cluster) {
		boolean flag = false;
		for (Domain domain : cluster.domains())
			if (domain.type() == gene)
				flag = true;
		return flag; 
	}

	/**
	 * Determine whether a list of deoxy sugar genes contains a gene of a particular type.
	 * @param domain	query sugar gene
	 * @param genes		list of deoxy sugar genes
	 * @return			true of the list contains a gene of the query type 
	 */
	public static boolean contains(Domain domain, List<DeoxySugarDomains> genes) {
		boolean flag = false;
		for (DeoxySugarDomains gene : genes)
			if (domain.type() == gene)
				flag = true;
		return flag;
	}

	/**
	 * Determine whether this glycosyltransferase domain represents a hexose-specific glycosyltransferase. 
	 * @param domain	glycosyltransferase domain to analyze
	 * @return			true if the top BLASTp hit is to a hexose GTr
	 */
	public static boolean isHexoseGlycosyltransferase(Domain domain) {
		boolean flag = false;
		if (domain.type() == TailoringDomains.GLYCOSYLTRANSFERASE && domain.blastResults().size() > 0) {
			BlastSearchResult top = domain.blastResults().get(0);
			String name = top.subject().toLowerCase();
			if (name.contains("glucose") || name.contains("mannose") || name.contains("gulose") 
					|| name.contains("glcnac"))
				flag = true;
		}
		return flag;
	}

	/**
	 * Convert a BLASTp search result subject string to a hexose sugar enum. 
	 * @param domain	glycosyltransferase domain to analyze
	 * @return			hexose sugar enum corresponding to this hexose sugar 
	 */
	public static HexoseSugars getHexoseSugarType(Domain domain) {
		HexoseSugars type = null;
		if (domain.type() == TailoringDomains.GLYCOSYLTRANSFERASE && domain.blastResults().size() > 0) {
			BlastSearchResult top = domain.blastResults().get(0);
			String name = top.subject().toLowerCase();
			if (name.contains("glucose"))
				type = HexoseSugars.GLUCOSE;
			if (name.contains("mannose")) 
				type = HexoseSugars.MANNOSE;
			if (name.contains("gulose"))
				type = HexoseSugars.GULOSE;
			if (name.contains("glcnac"))
				type = HexoseSugars.GLCNAC;
		}
		return type;
	}

}
