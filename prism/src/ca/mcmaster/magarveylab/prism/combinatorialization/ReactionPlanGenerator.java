package ca.mcmaster.magarveylab.prism.combinatorialization;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.exception.CDKException;

import ca.mcmaster.magarveylab.enums.ReactionPriorities;
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
	 * Get all combinations of reaction plans for a given permutation of biosynthetic modules within a biosynthetic 
	 * gene cluster.
	 * @param permutation	a given biosynthetic module permutation
	 * @param cluster		the cluster corresponding to this permutation
	 * @return				a list of lists, where each inner list represents a reaction plan combination 
	 * @throws IOException
	 * @throws ClassInstantiationException 
	 * @throws CDKException 
	 */
	public static List<List<ReactionPlan>> generateReactionPlans(List<Module> permutation, Cluster cluster) 
			throws IOException, ClassInstantiationException, CDKException {
		List<ReactionPlan> plans = getReactionPlans(permutation, cluster);
		List<List<ReactionPlan>> combinations = combinatorializeReactionPlans(plans);
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
		for (Domain domain : cluster.domains()) {
			List<SubstrateSet> substrates = AnnotatorUtil.findSubstrates(domain, permutation, cluster);
			for (SubstrateSet substrate : substrates) {
				ReactionPriorities priority = ReactionUtil.getPriority(domain);
				if (priority != null) {
					ReactionPlan reactionPlan = new ReactionPlan(domain, substrate, priority);
					plans.add(reactionPlan);
				}
			}
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
				ReactionPriorities.P450A,
				ReactionPriorities.P450B,
				ReactionPriorities.P450C,
				ReactionPriorities.P450D,
				ReactionPriorities.IPNS,
				ReactionPriorities.IAT,
				ReactionPriorities.DOACS,
				ReactionPriorities.ABM,
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
				ReactionPriorities.THIRD_RING_CYCLASE_3
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
						itr.remove();
					}
		}
	}
	
	/**
	 * Combinatorialize a raw list of reaction plans.
	 * @param plans		all possible reaction plans for a given biosynthetic module permutation 
	 * @return			a combinatorialized list of reaction plan combinations 
	 */
	public static List<List<ReactionPlan>> combinatorializeReactionPlans(List<ReactionPlan> plans) {
		List<List<ReactionPlan>> input = new ArrayList<List<ReactionPlan>>();
		List<List<ReactionPlan>> output = new ArrayList<List<ReactionPlan>>();

		Map<Domain,List<ReactionPlan>> map = convertReactionPlansToMap(plans);
		for (Map.Entry<Domain,List<ReactionPlan>> entry : map.entrySet()) 
			input.add(entry.getValue());
		
		Combinations.combinations(input, output, 0, new ArrayList<ReactionPlan>());
		return output; 
	}
	
	/**
	 * Convert a raw list of all reaction domains, unsorted by reaction domain, into a map of reaction domains
	 * to plans matching them to a potential substrate. 
	 * @param plans		all possible reaction plans for a given biosynthetic module permutation
	 * @return			a map of reaction domains to their reaction plans 
	 */
	public static Map<Domain,List<ReactionPlan>> convertReactionPlansToMap(List<ReactionPlan> plans) {
		Map<Domain,List<ReactionPlan>> map = new HashMap<Domain,List<ReactionPlan>>();
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
