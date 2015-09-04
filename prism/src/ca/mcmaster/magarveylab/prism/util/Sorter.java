package ca.mcmaster.magarveylab.prism.util;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import ca.mcmaster.magarveylab.enums.DeoxySugars;
import ca.mcmaster.magarveylab.prism.blast.BlastSearchResult;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Contig;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.data.Substrate;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.sugar.Sugar;
import ca.mcmaster.magarveylab.prism.homology.data.HomologousCluster;
import ca.mcmaster.magarveylab.prism.tanimoto.data.TanimotoScore;

/**
 * Custom sorting functions.
 * 
 * @author skinnider
 *
 */
public class Sorter {
	
	public static void sortRibosomalModulesByIndex(List<Module> modules) {
		Collections.sort(modules, new Comparator<Module>() {
			@Override
			public int compare(Module m1, Module m2) {
				return Integer.compare(Integer.parseInt(m1.scaffold().name().split("_")[1]), 
						Integer.parseInt(m2.scaffold().name().split("_")[1]));
			}
		});
	}

	/**
	 * Sort a list of contigs by the size of their nucleotide sequences.
	 * 
	 * @param contigs
	 *            contigs to sort
	 */
	public static void sortContigsBySize(List<Contig> contigs) {
		Collections.sort(contigs, new Comparator<Contig>() {
			@Override
			public int compare(Contig c1, Contig c2) {
				return Integer.compare(c1.length(), c2.length()) * -1;
			}
		});
	}

	/**
	 * Sort a list of substrates by their scores.
	 * 
	 * @param substrates
	 *            substrates to sort
	 * @return sorted list of substrates
	 */
	public static void sortSubstrates(List<Substrate> substrates) {
		Collections.sort(substrates, new Comparator<Substrate>() {
			@Override
			public int compare(Substrate s1, Substrate s2) {
				return Double.compare(s1.score(), s2.score()) * -1;
			}
		});
	}

	/**
	 * Sort a list of domains by start point.
	 * 
	 * @param domains
	 *            domains to sort
	 */
	public static void sortDomains(List<Domain> domains) {
		Collections.sort(domains, new Comparator<Domain>() {
			@Override
			public int compare(Domain d1, Domain d2) {
				return Integer.compare(d1.start(), d2.start());
			}
		});
	}

	/**
	 * Sort a list of domains by score (highest first).
	 * 
	 * @param domains
	 *            domains to sort
	 */
	public static void sortDomainsByScore(List<Domain> domains) {
		Collections.sort(domains, new Comparator<Domain>() {
			@Override
			public int compare(Domain d1, Domain d2) {
				return Double.compare(d1.score(), d2.score()) * -1;
			}
		});
	}

	/**
	 * Sort a list of domains by score to a reference domain (highest first).
	 * 
	 * @param domains
	 *            domains to sort
	 */
	public static void sortDomainsByReferenceScore(List<Domain> domains) {
		Collections.sort(domains, new Comparator<Domain>() {
			@Override
			public int compare(Domain d1, Domain d2) {
				return Double.compare(d1.referenceScore(), d2.referenceScore()) * -1;
			}
		});
	}

	/**
	 * Sort a list of orfs by start point.
	 * 
	 * @param orfs
	 *            orfs to sort
	 */
	public static void sortOrfs(List<Orf> orfs) {
		Collections.sort(orfs, new Comparator<Orf>() {
			@Override
			public int compare(Orf o1, Orf o2) {
				return Integer.compare(o1.start(), o2.start());
			}
		});
	}

	/**
	 * Sort a list of modules by scaffold domain start (or, in the absence of a
	 * scaffold domain, the first module).
	 * 
	 * @param modules
	 *            modules to sort
	 */
	public static void sortModules(List<Module> modules) {
		Collections.sort(modules, new Comparator<Module>() {
			@Override
			public int compare(Module m1, Module m2) {
				Domain s1 = m1.scaffold();
				if (s1 == null)
					s1 = m1.first();
				Domain s2 = m2.scaffold();
				if (s2 == null)
					s2 = m2.first();
				return Integer.compare(s1.start(), s2.start());
			}
		});
	}

	/**
	 * Sort a list of Tanimoto scores by the ECFP6 coefficient (highest first).
	 * 
	 * @param scores
	 *            scores to sort
	 */
	public static void sortTanimotoScoresByEcfp6Coefficient(List<TanimotoScore> scores) {
		Collections.sort(scores, new Comparator<TanimotoScore>() {
			@Override
			public int compare(TanimotoScore s1, TanimotoScore s2) {
				return Float.compare(s1.score("ecfp6"), s2.score("ecfp6")) * -1;
			}
		});
	}

	/**
	 * Sort a list of Blastp results by bit score (highest first).
	 * 
	 * @param results
	 *            results to sort
	 */
	public static void sortBlastpResults(List<BlastSearchResult> results) {
		Collections.sort(results, new Comparator<BlastSearchResult>() {
			@Override
			public int compare(BlastSearchResult r1, BlastSearchResult r2) {
				return Double.compare(r1.score(), r2.score()) * -1;
			}
		});
	}

	/**
	 * Sort all Blastp results. (Currently limited to sorting condensation
	 * results, to remove NaPDoS/C-starter overlap.)
	 * 
	 * @param contig
	 *            contig to sort
	 */
	public static void sortBlastpResults(Contig contig) {
		for (Orf orf : contig.orfs())
			for (Domain domain : orf.domains())
				if (domain.blastResults().size() > 0)
					sortBlastpResults(domain.blastResults());
	}

	/**
	 * Sort a list of homologous clusters by identity score (highest first).
	 * 
	 * @param homologs
	 *            homologous clusters to sort
	 */
	public static void sortHomologousClustersByIdentityScore(List<HomologousCluster> homologs) {
		Collections.sort(homologs, new Comparator<HomologousCluster>() {
			@Override
			public int compare(HomologousCluster h1, HomologousCluster h2) {
				return Double.compare(h1.identityScore(), h2.identityScore())
						* -1;
			}
		});
	}

	/**
	 * Sort a list of homologous clusters by domain score (highest first).
	 * 
	 * @param homologs
	 *            homologous clusters to sort
	 */
	public static void sortHomologousClustersByDomainScore(List<HomologousCluster> homologs) {
		Collections.sort(homologs, new Comparator<HomologousCluster>() {
			@Override
			public int compare(HomologousCluster h1, HomologousCluster h2) {
				return Double.compare(h1.domainScore(), h2.domainScore()) * -1;
			}
		});
	}

	/**
	 * Sort a list of homologous clusters by average domain score (highest first).
	 * 
	 * @param homologs
	 *            homologous clusters to sort
	 */
	public static void sortHomologousClustersByAverageDomainScore(List<HomologousCluster> homologs) {
		Collections.sort(homologs, new Comparator<HomologousCluster>() {
			@Override
			public int compare(HomologousCluster h1, HomologousCluster h2) {
				return Double.compare(h1.averageDomainScore(),
						h2.averageDomainScore()) * -1;
			}
		});
	}

	/**
	 * Sort a list of sugars by their name (alphabetical).
	 * 
	 * @param sugars
	 *            sugars to sort
	 */
	public static void sortSugarsByName(List<DeoxySugars> sugars) {
		Collections.sort(sugars, new Comparator<DeoxySugars>() {
			@Override
			public int compare(DeoxySugars s1, DeoxySugars s2) {
				return s1.toString().compareTo(s2.toString());
			}
		});
	}

	/**
	 * Sort a list of sugar combinations by the name of the first sugar
	 * (alphabetical).
	 * 
	 * @param combinations
	 *            sugar combinations to sort
	 */
	public static void sortSugarCombinationsByName(
			List<List<Sugar>> combinations) {
		Collections.sort(combinations, new Comparator<List<Sugar>>() {
			@Override
			public int compare(List<Sugar> s1, List<Sugar> s2) {
				if (s1 == null)
					return -1;
				if (s2 == null)
					return 1;
				return s1.get(0).toString().compareTo(s2.get(0).toString());
			}
		});
	}

	/**
	 * Sort a list of reaction plans by their priority, i.e. the order in which
	 * they should be executed.
	 * 
	 * @param plans
	 *            reaction plans to sort
	 */
	public static void sortPlansByPriority(List<ReactionPlan> plans) {
		Collections.sort(plans, new Comparator<ReactionPlan>() {
			@Override
			public int compare(ReactionPlan p1, ReactionPlan p2) {
				return Integer.compare(p1.priority(), p2.priority());
			}
		});
	}

}
