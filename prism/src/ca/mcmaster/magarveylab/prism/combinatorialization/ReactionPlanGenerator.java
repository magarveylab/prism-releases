package ca.mcmaster.magarveylab.prism.combinatorialization;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.openscience.cdk.exception.CDKException;

import ca.mcmaster.magarveylab.enums.ReactionPriorities;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.annotation.AnnotatorUtil;
import ca.mcmaster.magarveylab.prism.cluster.reactions.ReactionUtil;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.util.exception.ClassInstantiationException;

/**
 * Generates reaction plans for combinatorialization. 
 * @author skinnider
 *
 */
public class ReactionPlanGenerator {

	/**
	 * Get all combinations of reaction plans for a given permutation of
	 * biosynthetic modules within a biosynthetic gene cluster.
	 * 
	 * @param permutation
	 *            a given biosynthetic module permutation
	 * @param cluster
	 *            the cluster corresponding to this permutation
	 * @return a list of lists, where each inner list represents a reaction plan
	 *         combination
	 * @throws IOException
	 * @throws ClassInstantiationException
	 * @throws CDKException
	 */
	public static List<List<ReactionPlan>> generateReactionPlans(
			List<Module> permutation, Cluster cluster) throws IOException,
			ClassInstantiationException, CDKException {
		List<ReactionPlan> plans = getReactionPlans(permutation, cluster);
		List<List<ReactionPlan>> combinations = combinatorializeReactionPlans(
				plans, permutation);
		return combinations;
	}
	
	/**
	 * Get all possible reaction plans for a given permutation of biosynthetic
	 * modules within a biosynthetic gene clusters, for combinatorialization.
	 * 
	 * @param permutation
	 *            a given biosynthetic module permutation
	 * @param cluster
	 *            the cluster corresponding to this permutation
	 * @return a list of reaction plans
	 * @throws IOException
	 * @throws ClassInstantiationException
	 * @throws CDKException
	 */
	public static List<ReactionPlan> getReactionPlans(List<Module> permutation, Cluster cluster) 
			throws IOException, ClassInstantiationException, CDKException {
		List<ReactionPlan> plans = new ArrayList<ReactionPlan>();

		DomainType[] unique = new DomainType[] {
				// domains that can only act once per cluster
				RibosomalDomains.LanB, RibosomalDomains.LanC,
				RibosomalDomains.LanM, RibosomalDomains.LanKC, RibosomalDomains.PatD,
				RibosomalDomains.McbB, RibosomalDomains.McbC,
				RibosomalDomains.McbD, RibosomalDomains.LazB,
				RibosomalDomains.LazC, RibosomalDomains.LazC_b,
				RibosomalDomains.LazE, RibosomalDomains.LazF,
				RibosomalDomains.YmBC_a, RibosomalDomains.YmBC_b,
				RibosomalDomains.BotCD, RibosomalDomains.GodF,
				RibosomalDomains.GodG, RibosomalDomains.CypH,
				RibosomalDomains.CypL, RibosomalDomains.PoyBC,
				RibosomalDomains.PoyF, RibosomalDomains.PoyE,
				RibosomalDomains.PoyI, RibosomalDomains.AlbA,
				RibosomalDomains.TfxB, RibosomalDomains.SunA,
				RibosomalDomains.SunS, RibosomalDomains.PatF,
				RibosomalDomains.ProcA };
		boolean[] used = new boolean[unique.length];
		Arrays.fill(used, false);

		for (Domain domain : cluster.domains()) {
			int idx = Arrays.asList(unique).indexOf(domain.type());
			if (idx != -1) {
				if (used[idx] == true) {
					continue;
				} else {
					used[idx] = true;
				}
			}
			
			List<SubstrateSet> substrates = AnnotatorUtil.findSubstrates(
					domain, permutation, cluster);
			for (SubstrateSet substrate : substrates) {
				ReactionPriorities priority = ReactionUtil.getPriority(domain);
				if (priority != null) {
					ReactionPlan reactionPlan = new ReactionPlan(domain,
							substrate, priority);
					plans.add(reactionPlan);
				}
			}

			System.out.println("Identified " + substrates.size()
					+ " substrates for " + domain.type());
		}
		checkDuplicateReactionPlans(plans);
		return plans;
	}
	
	/**
	 * Remove duplicated reaction plans which correspond to reactions that can
	 * only occur once, e.g. isopenicillin N synthase, isopenicillin N
	 * acyltransferase, deacetoxycephalosporin C synthase, P450s, or polyketide
	 * cyclases.
	 * 
	 * @param plans
	 *            a list of reaction plans
	 */
	public static void checkDuplicateReactionPlans(List<ReactionPlan> plans) {
		ReactionPriorities[] redundant = new ReactionPriorities[] {
				ReactionPriorities.P450A, ReactionPriorities.P450B,
				ReactionPriorities.P450C, ReactionPriorities.P450D,
				ReactionPriorities.IPNS, ReactionPriorities.IAT,
				ReactionPriorities.DOACS, ReactionPriorities.ABM,
				ReactionPriorities.ANGULAR_CYCLASE_1,
				ReactionPriorities.ANGULAR_CYCLASE_2,
				ReactionPriorities.C9_KETOREDUCTASE,
				ReactionPriorities.C11_KETOREDUCTASE,
				ReactionPriorities.C15_KETOREDUCTASE,
				ReactionPriorities.C17_KETOREDUCTASE,
				ReactionPriorities.C19_KETOREDUCTASE,
				ReactionPriorities.FIRST_RING_CYCLASE_1,
				ReactionPriorities.FIRST_RING_CYCLASE_2,
				ReactionPriorities.FIRST_RING_CYCLASE_3,
				ReactionPriorities.FOURTH_RING_CYCLASE_1,
				ReactionPriorities.FOURTH_RING_CYCLASE_2,
				ReactionPriorities.FOURTH_RING_CYCLASE_3,
				ReactionPriorities.FOURTH_RING_CYCLASE_4,
				ReactionPriorities.THIRD_RING_CYCLASE_1,
				ReactionPriorities.THIRD_RING_CYCLASE_2,
				ReactionPriorities.THIRD_RING_CYCLASE_3,
		};
		boolean[] used = new boolean[redundant.length];
		Arrays.fill(used, false);
		
		Iterator<ReactionPlan> itr = plans.iterator();
		while (itr.hasNext()) {
			ReactionPlan next = itr.next();	
			ReactionPriorities type = next.reaction();
			for (int i = 0; i < redundant.length; i++) 
				if (type == redundant[i])
					if (used[i] == false) {
						used[i] = true;
					} else {
						System.out.println("Removing " + type + " reaction plan");
						itr.remove();
					}
		}
	}

	/**
	 * Combinatorialize a raw list of reaction plans.
	 * 
	 * @param plans
	 *            all possible reaction plans for a given biosynthetic module
	 *            permutation
	 * @return a combinatorialized list of reaction plan combinations
	 */
	public static List<List<ReactionPlan>> combinatorializeReactionPlans(
			List<ReactionPlan> plans, List<Module> permutation) {
		System.out.println("[ReactionPlanGenerator] Combinatorializing "
				+ "a list of " + plans.size() + " reaction plans");

		List<List<ReactionPlan>> input = new ArrayList<List<ReactionPlan>>();
		Map<Domain, List<ReactionPlan>> map = convertReactionPlansToMap(plans);
		for (Map.Entry<Domain, List<ReactionPlan>> entry : map.entrySet()) {
			input.add(entry.getValue());
			System.out.println("Generated map key for domain "
					+ entry.getKey().type() + " with "
					+ entry.getValue().size() + " plans");
		}

		long total = 1;
		for (List<ReactionPlan> p : map.values()) {
			total *= p.size();
		}
		
		List<List<ReactionPlan>> output = new ArrayList<List<ReactionPlan>>();
		if (Math.abs(total) > 1000) {
			System.out.println("[ReactionPlanGenerator] Sampling from a space of " + total);
			output = sampleReactionPlanCombinations(input, permutation, 1_000);
		} else {
			System.out.println("[ReactionPlanGenerator] Generating all combinations from a space of " + total);
			Combinations.combinations(input, output, 0, new ArrayList<ReactionPlan>());
		}
		
		System.out.println("[ReactionPlanGenerator] Generated " + output.size()
				+ " combinations of reaction plans");
		
		return output;
	}
	
	public static List<List<ReactionPlan>> sampleReactionPlanCombinations(
			List<List<ReactionPlan>> lists, List<Module> permutation, int limit) {
		Random random = new Random(0);
		int i = 0, used = 0, failed = 0;
		List<List<ReactionPlan>> output = new ArrayList<List<ReactionPlan>>();
		List<int[]> combinations = new ArrayList<int[]>();
		while (i < limit
				&& used < 150_000
				&& failed < 100_000_000) {
			if (failed % 5_000_000 == 0 && failed > 0)
				System.out.println("[ReactionPlanGenerator] Failed plan " + failed + " (" + i + " successes)");
			
			int[] combination = new int[lists.size()];
			for (int j = 0; j < lists.size(); j++) {
				List<ReactionPlan> list = lists.get(j);
				int r = random.nextInt(list.size());
				combination[j] = r;
			}
			if (!Combinations.usedCombination(combination, combinations)
					|| combinations.size() == 0) {
				// get reaction plans
				List<ReactionPlan> list = new ArrayList<ReactionPlan>();
				for (int a = 0; a < combination.length; a++) {
					int c = combination[a];
					list.add(lists.get(a).get(c));
				}		

				// check overlap
				ReactionOverlapChecker roc = new ReactionOverlapChecker(list,
						permutation);
				if (roc.check()) {
					combinations.add(combination);
					output.add(list);
					i++;
				} else {
					failed++;
				}
			} else {
				used++;
			}
			
			if (used >= 150_000)
				System.out.println("[ReactionPlanGenerator] Exiting reaction plan generation: too many plans used");
			if (failed >= 100_000_000)
				System.out.println("[ReactionPlanGenerator] Exiting reaction plan generation: too many plans failed");
		}
		
		System.out.println("[ReactionPlanGenerator] Generated " + output.size() + " sampled reaction plan combinations");
		
		return output;
	}

	/**
	 * Convert a raw list of all reaction domains, unsorted by reaction domain,
	 * into a map of reaction domains to plans matching them to a potential
	 * substrate.
	 * 
	 * @param plans
	 *            all possible reaction plans for a given biosynthetic module
	 *            permutation
	 * @return a map of reaction domains to their reaction plans
	 */
	public static Map<Domain, List<ReactionPlan>> convertReactionPlansToMap(
			List<ReactionPlan> plans) {
		Map<Domain, List<ReactionPlan>> map = new HashMap<Domain,List<ReactionPlan>>();
		for (ReactionPlan plan : plans) {
			if (map.containsKey(plan.domain())) {
				map.get(plan.domain()).add(plan);
			} else {
				List<ReactionPlan> list = new ArrayList<ReactionPlan>();
				list.add(plan);
				map.put(plan.domain(), list);
			}
		}
		return map;
	}

}
