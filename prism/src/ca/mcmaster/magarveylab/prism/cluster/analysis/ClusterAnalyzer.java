package ca.mcmaster.magarveylab.prism.cluster.analysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import ca.mcmaster.magarveylab.enums.ClusterFamilies;
import ca.mcmaster.magarveylab.enums.Frames;
import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.SubstratePrerequisites;
import ca.mcmaster.magarveylab.enums.clusters.ClusterType;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.enums.substrates.AdenylationSubstrates;
import ca.mcmaster.magarveylab.enums.substrates.SubstrateType;
import ca.mcmaster.magarveylab.prism.cluster.analysis.type.ThiotemplatedClusterTypeAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.analysis.type.TypeIIPolyketideClusterTypeAnalyzer;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.data.Substrate;
import ca.mcmaster.magarveylab.prism.homology.data.HomologousCluster;
import ca.mcmaster.magarveylab.prism.tanimoto.data.TanimotoScore;
import ca.mcmaster.magarveylab.prism.util.Sorter;
import ca.mcmaster.magarveylab.prism.util.Strings;
import ca.mcmaster.magarveylab.prism.web.PrismConfig;

/**
 * Analyzes clusters for combinatorial library generation.
 * @author skinnider
 *
 */
public class ClusterAnalyzer {

	/**
	 * Get the orf that starts this cluster. Will return FAAL > C* > non-nitrogenous A > starter A > KS1.
	 * @param orfs		list of orfs
	 * @return			start orf, or null if no discernible start orf
	 */
	public static Orf start(List<Orf> orfs) {
		Orf start = null;
		boolean cstarter = false;
		boolean acylAdenylate = false;

		for (Orf orf : orfs) {
			if (OrfAnalyzer.isStarterKS(orf) && !cstarter && !acylAdenylate)
				start = orf;
			if (OrfAnalyzer.isCStarter(orf) && !acylAdenylate) {
				start = orf;
				cstarter = true;
			}
			if (OrfAnalyzer.isStarterAcylAdenylate(orf)) {
				if (acylAdenylate) { // if there are two starter units, there is no start
					start = null;
					break;
				}
				start = orf;
				acylAdenylate = true;
			}
		}

		if (start != null)
			System.out.println("[ClusterAnalyzer] Located start orf: " + start.name());
		return start;
	}
	
	/**
	 * Get the orf that ends this cluster. 
	 * @param orfs		list of orfs
	 * @return			end orf, or null if no discernible end orf
	 */
	public static Orf end(List<Orf> orfs) {
		Orf end = null;
		List<Orf> ends = new ArrayList<Orf>();
		for (Orf orf : orfs) 
			if (OrfAnalyzer.isEnd(orf))
				ends.add(orf);
		if (ends.size() == 1)
			end = ends.get(0);
		if (end != null)
			System.out.println("[ClusterAnalyzer] Located end orf: " + end.name());
		return end;
	}

	/**
	 * Detect the frame of the majority of this cluster's scaffold orfs for use in scaffold construction.
	 * @param cluster
	 */
	public static void detectFrame(Cluster cluster) {
		int positive = 0, negative = 0;
		List<Orf> scaffoldOrfs = cluster.moduleOrfs();
		for (Orf orf : scaffoldOrfs) {
			if (orf.frame().indexOf("+") != -1) {
				positive++;
			} else {
				negative++;
			}
		}
		Frames frame = Frames.NONE;
		if (positive > 0 && negative == 0)
			frame = Frames.POSITIVE;
		if (positive == 0 && negative > 0)
			frame = Frames.NEGATIVE;

		// check for start/end enzymes
		Orf start = ClusterAnalyzer.start(cluster.orfs());
		Orf end = ClusterAnalyzer.end(cluster.orfs());
		List<Orf> moduleOrfs = cluster.moduleOrfs();
		if (frame == Frames.NEGATIVE)
			Collections.reverse(moduleOrfs);
		if (start != null && moduleOrfs.indexOf(start) != 0) {
			System.out.println("[ClusterAnalyzer] Error: biosynthetic start orf " + start.name() + " is not at the start of "
					+ "cluster " + cluster.index() + ": index " + moduleOrfs.indexOf(start) 
					+ ". Changed frame from " + frame + " to NONE.");
			frame = Frames.NONE;
		} else if (end != null && moduleOrfs.indexOf(end) != moduleOrfs.size() - 1) {
			System.out.println("[ClusterAnalyzer] Error: biosynthetic end orf " + end.name() + " is not at the end of "
					+ "cluster " + cluster.index() + ": index " + moduleOrfs.indexOf(end) + ", size " + moduleOrfs.size()
					+ ". Changed frame from " + frame + " to NONE.");
			frame = Frames.NONE;
		}

		cluster.setFrame(frame);
		System.out.println("[ClusterAnalyzer] Set cluster " + cluster.index() + " frame to " + cluster.frame().toString());
	}
	
	/**
	 * Analyze cluster to determine whether modules can be extended within the growing natural product scaffold.
	 * Fatty acids, for instance, can only start biosynthesis. 
	 * @param cluster	cluster to analyze 
	 */
	public static void setExtendability(Cluster cluster) {
		for (Module module : cluster.modules()) {
			Domain scaffold = module.scaffold();
			if (scaffold == null || scaffold.substrates().size() == 0)
				continue;
			Substrate top = scaffold.topSubstrate();
			String smiles = top.smiles();
			if (smiles.indexOf("F") == -1)
				module.setCanExtend(false);
		}
	}

	/**
	 * Remove prolyl-AMP ligases if the cluster does not contain a proline dehydrogenase domain. 
	 * @param cluster	cluster to analyze 
	 */
	public static void checkPyrroleModules(Cluster cluster) {
		if (!cluster.contains(TailoringDomains.PROLINE_DEHYDROGENASE)) {
			for (Orf orf : cluster.orfs()) {
				Iterator<Module> itr = orf.modules().iterator();
				while (itr.hasNext()) {
					Module next = itr.next();
					if (next.type() == ModuleTypes.ACYL_ADENYLATE 
							&& next.scaffold() != null 
							&& next.scaffold().topSubstrate().type() == AdenylationSubstrates.PROLINE_3)
						itr.remove();
				}
			}
		}
	}

	/**
	 * Analyze cluster to determine whether all A-T, or trans-adenylation, modules actually work in trans, or merely
	 * represent an adenylation module with condensation domain absent. (Ensure there are both insertion and trans-A
	 * modules, and that trans-A modules are on their own orfs). 
	 * @param cluster	cluster to analyze 
	 */
	public static void checkTransAdenylationModules(Cluster cluster) {
		List<Module> insertion = cluster.modules(ModuleTypes.TRANS_ADENYLATION_INSERTION);
		System.out.println("[ClusterAnalyzer] Found " + insertion.size() + " adenylation insertion modules");
		if (insertion.size() > 1) 
			for (int i = 1; i < insertion.size(); i++)
				insertion.get(i).inactivate();

		for (Orf orf : cluster.orfs()) {
			List<Module> orfTransA = orf.modules(ModuleTypes.TRANS_ADENYLATION);
			List<Module> orfInsertion = orf.modules(ModuleTypes.TRANS_ADENYLATION_INSERTION);

			if (orf.modules().size() == 1 && insertion.size() == 0) {
				// if a trans A domain is on its own orf AND there are no insertion modules, it is a starter A domain
				for (Module module : orfInsertion) {
					module.inactivate();
					System.out.println("[ClusterAnalyzer] Changed trans-adenylation insertion module on orf " + orf.name() 
							+ " to " + module.type());
				}
				for (Module module : orfTransA) {
					module.setType(ModuleTypes.ADENYLATION);
					System.out.println("[ClusterAnalyzer] Changed trans-adenylation module on orf " + orf.name() 
							+ " to " + module.type());
				}
			} else if (orf.modules().size() > 1) {
				// if a trans A domain isn't on its own orf, it can't be a trans A domain
				for (int i = 0; i < orfTransA.size(); i++) {
					Module module = orfTransA.get(i);
					module.setType(ModuleTypes.ADENYLATION);
					System.out.println("[ClusterAnalyzer] Changed trans-adenylation module on orf " + orf.name() 
							+ " to " + module.type());
				}
			}
		}
	}

	/**
	 * Analyze cluster to inactivate trans-AT insertion modules in the absence of a trans-AT domain. 
	 * @param cluster	cluster to analyze
	 */
	public static void checkTransAcyltransferaseModules(Cluster cluster) {
		List<Module> transATInsertion = cluster.modules(ModuleTypes.TRANS_AT_INSERTION);
		List<Module> transAT = cluster.modules(ModuleTypes.TRANS_AT);
		if (transAT.size() == 0 && transATInsertion.size() > 0)
			for (Module module : transATInsertion)
				module.inactivate();
	}

	/**
	 * Ensure that a cluster does not have overlapping start modules (i.e., C-starter, salicylic acid/
	 * phenylacetate adenylation, and FAAL modules) as this will cause an error. Priority: FAAL > C*. 
	 * @param cluster	cluster to check
	 */
	public static void checkStarterModules(Cluster cluster) {
		List<Module> cStarter = cluster.modules(ModuleTypes.C_STARTER);
		List<Module> faal = cluster.modules(ModuleTypes.ACYL_ADENYLATE);
		
		// only start!!
		Iterator<Module> itr = faal.iterator();
		while (itr.hasNext()) {
			Module next = itr.next();
			if (next.canExtend())
				itr.remove();
		}

		// fix C*/FAAL overlap
		if (faal.size() > 0 && cStarter.size() > 0)
			for (Module module : cStarter) {
				System.out.println("[ClusterAnalyzer] Inactivated C-starter module with " + module.scaffold().name());
				module.inactivate();		
			}
	}

	/**
	 * Inactivate modules with a scaffold domain that hits to a hidden Markov model associated with non-proteinogenic 
	 * amino acid or starter unit biosynthesis. 
	 * @param cluster	cluster to analyze 
	 */
	public static void checkDockingModules(Cluster cluster) {
		for (Orf orf : cluster.orfs()) 
			for (Module module : orf.modules()) 
				if (module.scaffold() != null && module.scaffold().type() == ThiotemplatedDomains.ADENYLATION
						&& (module.scaffold().topSubstrate().type() 
								== AdenylationSubstrates.GLYCOPEPTIDE_STARTER_UNIT_BIOSYNTHESIS
						|| module.scaffold().topSubstrate().type() 
								== AdenylationSubstrates.QUINOMYCIN_STARTER_UNIT_BIOSYNTHESIS)) {
					module.inactivate();
					System.out.println("[ClusterAnalyzer] Inactivated putative docking domain " 
							+ module.scaffold().topSubstrate().type() + " in " + orf.name());
				}
	}

	/**
	 * Check substrates which require prerequisite domains to be present to be definitively assigned. 
	 * For instance, detection of enduracididine requires a PLP-dependent aminotransferase. 
	 * @param cluster	cluster to analyze 
	 */
	public static void checkSubstratePrerequisites(Cluster cluster) {
		for (SubstratePrerequisites prerequisites : SubstratePrerequisites.values()) {
			SubstrateType substrate = prerequisites.substrate();
			DomainType[] domains = prerequisites.domains();
			
			List<Domain> substrateDomains = new ArrayList<Domain>();
			for (Domain domain : cluster.domains())
				if (domain.topSubstrate() != null && domain.topSubstrate().type() == substrate)
					substrateDomains.add(domain);

			boolean containsPrerequisites = true;
			for (DomainType domainType : domains)
				if (cluster.domains(domainType).size() == 0)
					containsPrerequisites = false;

			if (!containsPrerequisites)
				for (Domain domain : substrateDomains)
					domain.substrates().remove(0);
		}
	}

	/**
	 * Determine whether this cluster is a known biosynthetic gene cluster and get the name of its top match. Homology 
	 * takes priority over Tanimoto scoring. 
	 * @param cluster	cluster to check
	 * @param config	current PRISM configuration
	 * @return	the string "Unknown" if the cluster is unknown, or the name of the top-scoring cluster if known
	 */
	public static String isKnown(Cluster cluster, PrismConfig config) {
		List<TanimotoScore> tanimoto = cluster.scores();
		List<HomologousCluster> homologs = cluster.homologs();

		if (homologs.size() > 0) {
			HomologousCluster h = homologs.get(0);
			Sorter.sortHomologousClustersByIdentityScore(homologs);
			if (h.identityScore() > config.homologyCutoff)
				return Strings.removeUnderscores(Strings.capitalizeFirstLetter(h.name()));
			Sorter.sortHomologousClustersByDomainScore(homologs);
			h = homologs.get(0);
			if (h.averageDomainScore() > config.homologyCutoff
					&& h.coverage() > 0.25)
				return Strings.capitalize(h.name());
		}

		if (tanimoto.size() > 0) {
			TanimotoScore t = tanimoto.get(0);
			if (t.score("ecfp6") > config.tanimotoCutoff
					|| t.score("fcfp6") > config.tanimotoCutoff) 
				return Strings.capitalize(t.target().name());
		}

		return "Unknown";
	}

	/**
	 * Determine whether a cluster represents a trans-AT cluster, to execute
	 * special trans-AT-specific analysis.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if both trans AT domain and trans AT insertion modules found
	 */
	public static boolean isTransAT(Cluster cluster) {
		List<Module> transAT = cluster.modules(ModuleTypes.TRANS_AT);
		List<Module> transATInsertion = cluster.modules(ModuleTypes.TRANS_AT_INSERTION);
		return (transAT.size() > 0 && transATInsertion.size() > 0);
	}

	/**
	 * Assign a biosynthetic family (or list of biosynthetic families) and a
	 * family-specific subtype (or list of subtypes) to a natural product
	 * biosynthetic gene cluster.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 */
	public static void setClusterType(Cluster cluster) {
		List<ClusterType> typeII = TypeIIPolyketideClusterTypeAnalyzer.getTypes(cluster);
		List<ClusterType> nrps = ThiotemplatedClusterTypeAnalyzer.getNonribosomalPeptideTypes(cluster);
		List<ClusterType> pks = ThiotemplatedClusterTypeAnalyzer.getPolyketideTypes(cluster);
		List<ClusterType> otherThiotemplated = ThiotemplatedClusterTypeAnalyzer.getOtherTypes(cluster);
		
		if (nrps.size() > 0) 
			cluster.addFamily(ClusterFamilies.NONRIBOSOMAL_PEPTIDE);
		if (pks.size() > 0) 
			cluster.addFamily(ClusterFamilies.TYPE_I_POLYKETIDE);
		if (typeII.size() > 0) 
			cluster.addFamily(ClusterFamilies.TYPE_II_POLYKETIDE);
		if (otherThiotemplated.size() > 0) 
			cluster.addFamily(ClusterFamilies.NULL);
		
		cluster.addTypes(typeII);
		cluster.addTypes(nrps);
		cluster.addTypes(pks);
		cluster.addTypes(otherThiotemplated);
	}

}
