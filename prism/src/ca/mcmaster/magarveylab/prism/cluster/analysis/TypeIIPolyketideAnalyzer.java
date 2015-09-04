package ca.mcmaster.magarveylab.prism.cluster.analysis;

import java.util.List;

import ca.mcmaster.magarveylab.enums.clusters.ClusterType;
import ca.mcmaster.magarveylab.enums.clusters.TypeIIPolyketideClusterTypes;
import ca.mcmaster.magarveylab.enums.domains.TypeIIPolyketideDomains;
import ca.mcmaster.magarveylab.enums.substrates.TypeIIPolyketideStarters;
import ca.mcmaster.magarveylab.prism.blast.BlastSearchResult;
import ca.mcmaster.magarveylab.prism.cluster.analysis.type.TypeIIPolyketideClusterTypeAnalyzer;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;

/**
 * Analyze type II polyketide clusters.
 * @author skinnider
 *
 */
public class TypeIIPolyketideAnalyzer {

	/**
	 * Get the number of times a starter unit is extended by the type II
	 * polyketide synthase ketosynthase heterodimer. Only the highest-scoring
	 * chain length factor domain will be considered.
	 * 
	 * @param cluster
	 *            type II PKS cluster
	 * @return number of extension cycles, in ketide units
	 */
	public static int getChainLength(Cluster cluster) {
		int length = -1;
		double score = -1.0d;
		for (Domain td : cluster.domains(TypeIIPolyketideDomains.CLF))
			if (td.score() > score) {
				BlastSearchResult topBlastHit = td.blastResults().get(0);
				String[] split = topBlastHit.subject().split("_");
				String extensions = split[split.length-1];
				length = Integer.parseInt(extensions);
			}
		return length; 
	}

	/**
	 * Get the starter unit of a type II polyketide based on the domain content of the 
	 * module. Returns acetate by default, or a specific starter unit if the cluster contains
	 * an amidotransferase, phenylalanine ammonia-lyase, or priming acyltransferase.
	 * @param module	starter module
	 * @return			type II polyketide starter unit 
	 */
	public static TypeIIPolyketideStarters getStarter(Module module) {
		TypeIIPolyketideStarters start = TypeIIPolyketideStarters.ACETATE;
		if (module.contains(TypeIIPolyketideDomains.AMIDOTRANSFERASE)) {
			start = TypeIIPolyketideStarters.MALONAMATE;
		} else if (module.contains(TypeIIPolyketideDomains.PAL)) {
			start = TypeIIPolyketideStarters.BENZOATE;
		} else if (module.contains(TypeIIPolyketideDomains.PRIMING_AT)) {
			double score = -1;
			for (Domain domain : module.domains())
				if (domain.type() == TypeIIPolyketideDomains.PRIMING_AT 
						&& domain.score() > score
						&& domain.blastResults().size() > 0) {
					String name = domain.blastResults().get(0).subject();
					if (name.contains("propionate")) {
						start = TypeIIPolyketideStarters.PROPIONATE;
					} else if (name.contains("butyrate")) {
						start = TypeIIPolyketideStarters.BUTYRATE;
					} else if (name.contains("isobutyrate")) {
						start = TypeIIPolyketideStarters.ISOBUTYRATE;
					} else if (name.contains("2-methylbutyrate")) {
						start = TypeIIPolyketideStarters._2_METHYLBUTYRATE;
					} else if (name.contains("hexadienoate")) {
						start = TypeIIPolyketideStarters.HEXADIENOATE;
					}
				}
		} 
		return start;
	}

	/**
	 * Type II polyketide cyclases of clade 6B can catalyze two different types of cyclization reactions. Determine,
	 * based on the CLF, which they catalyze and set the type of the 6B cyclase. 
	 * @param cluster
	 */
	public static void checkCyclaseClade6B(Cluster cluster) {
		if (cluster.contains(TypeIIPolyketideDomains.CYCLASE_CLADE_6b)) {
			List<ClusterType> types = TypeIIPolyketideClusterTypeAnalyzer.getTypes(cluster);

			Domain cyclase = null;
			for (Domain domain : cluster.domains())
				if (domain.type() == TypeIIPolyketideDomains.CYCLASE_CLADE_6b)
					cyclase = domain;

			if (cyclase != null && types.size() > 0) {
				if (types.contains(TypeIIPolyketideClusterTypes.ANTHRACYCLINE)) {
					cyclase.setType(TypeIIPolyketideDomains.CYCLASE_CLADE_6b_SUBTYPE_1);
				} else if (types.contains(TypeIIPolyketideClusterTypes.TETRACYCLINE)
						|| types.contains(TypeIIPolyketideClusterTypes.AUREOLIC_ACID)) {
					cyclase.setType(TypeIIPolyketideDomains.CYCLASE_CLADE_6b_SUBTYPE_2);
				}
			}
		}
	}

	/**
	 * In pentangular polyphenol clusters, the oxytetracycline C6
	 * C-methyltransferase in fact acts at C8. Determine the correct substrate
	 * from the polyketide cyclization pattern.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 */
	public static void checkC6Methyltransferase(Cluster cluster) {
		if (cluster.contains(TypeIIPolyketideDomains.C6CMT))
			if (cluster.contains(TypeIIPolyketideDomains.CYCLASE_CLADE_5a))
				for (Domain domain : cluster.domains(TypeIIPolyketideDomains.C6CMT))
					domain.setType(TypeIIPolyketideDomains.C8CMT);
	}

	/**
	 * Determine whether this cluster represents a type II polyketide cluster.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if both type II KSa and KSb domains are found in the cluster
	 */
	public static boolean isTypeIIPolyketideCluster(Cluster cluster) {
		return cluster.contains(TypeIIPolyketideDomains.KSA) 
				&& cluster.contains(TypeIIPolyketideDomains.CLF);
	}

}
