package ca.mcmaster.magarveylab.prism.cluster.analysis.type;

import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.enums.OrfTypes;
import ca.mcmaster.magarveylab.enums.clusters.ClusterType;
import ca.mcmaster.magarveylab.enums.clusters.ThiotemplatedClusterTypes;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Orf;

/**
 * Analyze thiotemplated clusters.
 * @author skinnider
 *
 */
public class ThiotemplatedClusterTypeAnalyzer {

	public static List<ClusterType> getPolyketideTypes(Cluster cluster) {
		List<ClusterType> types = new ArrayList<ClusterType>();

		// fungal type I?
		if (isFungalTypeICluster(cluster)) {
			types.add(ThiotemplatedClusterTypes.FUNGAL_TYPE_I);
			return types;
		}
		
		// enediyne?
		List<ClusterType> enediyne = getEnediyneTypes(cluster);
		if (enediyne.size() > 0) {
			types.addAll(enediyne);
			return types;
		}

		if (!isThiotemplatedCluster(cluster))
			return types; 

		boolean pks = false;
		for (Orf orf : cluster.orfs()) 
			if (orf.type() == OrfTypes.PKS)
				pks = true;
		if (pks == true)
			types.add(ThiotemplatedClusterTypes.PKS);
		
		return types;
	}
	
	public static List<ClusterType> getNonribosomalPeptideTypes(Cluster cluster) {
		List<ClusterType> types = new ArrayList<ClusterType>();

		if (!isThiotemplatedCluster(cluster))
			return types; 
		
		boolean nrps = false;
		List<Orf> orfs = cluster.orfs();
		for (Orf orf : orfs) 
			if (orf.type() == OrfTypes.NRPS)
				nrps = true;
		if (nrps == true) 
			types.add(ThiotemplatedClusterTypes.NRPS);
		
		return types;
	}
	
	public static List<ClusterType> getOtherTypes(Cluster cluster) {
		List<ClusterType> pks = getPolyketideTypes(cluster);
		List<ClusterType> nrps = getNonribosomalPeptideTypes(cluster);
		List<ClusterType> types = new ArrayList<ClusterType>();
		if (nrps.size() == 0 && pks.size() == 0 && isThiotemplatedCluster(cluster))
			types.add(ThiotemplatedClusterTypes.NULL);
		return types;
	}

	/**
	 * Determine whether this is a thiotemplated (NRPS, PKS, or hybrid) cluster.
	 * The method returns true if the cluster contains at least one scaffold
	 * (adenylation, acyltransferase, adenylation-ketoreductase) domain and one
	 * bond-forming (condensation or ketosynthase) domain.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if thiotemplated
	 */
	public static boolean isThiotemplatedCluster(Cluster cluster) {
		boolean scaffold = false;
		boolean cOrKS = false;
		for (Orf orf : cluster.orfs()) {
			if (orf.substrateDomains().size() > 0) {
				scaffold = true;
			}
			if (orf.domains(ThiotemplatedDomains.CONDENSATION).size() > 0
					|| orf.domains(ThiotemplatedDomains.KETOSYNTHASE).size() > 0) {
				cOrKS = true;
			}
		}
		return (cOrKS && scaffold);
	}
	
	/**
	 * Determine whether this is a fungal iterative type I polyketide synthase
	 * cluster. This method returns true if the cluster contains at least one
	 * ketosynthase domain, and at least one product template domain.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if fungal type I
	 */
	public static boolean isFungalTypeICluster(Cluster cluster) {
		return cluster.contains(ThiotemplatedDomains.KETOSYNTHASE)
				&& (cluster.contains(ThiotemplatedDomains.PRODUCT_TEMPLATE_I)
				|| cluster.contains(ThiotemplatedDomains.PRODUCT_TEMPLATE_II)
				|| cluster.contains(ThiotemplatedDomains.PRODUCT_TEMPLATE_III)
				|| cluster.contains(ThiotemplatedDomains.PRODUCT_TEMPLATE_IV)
				|| cluster.contains(ThiotemplatedDomains.PRODUCT_TEMPLATE_V));
	}

	/**
	 * Determine the type (9- and/or 10-membered ring) of a putative enediyne
	 * cluster.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return the enediyne cluster type(s), or an empty list if this cluster is
	 *         not an enediyne cluster
	 */
	public static List<ClusterType> getEnediyneTypes(Cluster cluster) {
		List<ClusterType> types = new ArrayList<ClusterType>();
		if (cluster.contains(ThiotemplatedDomains.ENEDIYNE_PPTASE)) 
			for (Domain domain : cluster.domains(ThiotemplatedDomains.KETOSYNTHASE))
				if (domain.blastResults().size() > 0) 
					if (domain.blastResults().get(0).subject().contains("_ene9")
							|| domain.blastResults().get(0).subject().contains("SP01_2830_1")
							|| domain.blastResults().get(0).subject().contains("SP02_0350_1")	
							|| domain.blastResults().get(0).subject().contains("Stro2697_1")) {
						types.add(ThiotemplatedClusterTypes.ENEDIYNE_9_MEMBERED);
					} else if (domain.blastResults().get(0).subject().contains("_ene10")) {
						types.add(ThiotemplatedClusterTypes.ENEDIYNE_10_MEMBERED);
					}
		return types;
	}

}
