package ca.mcmaster.magarveylab.prism.combinatorialization;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.CDKException;

import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.prism.cluster.analysis.CyclizationAnalyzer;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.CombinatorialData;
import ca.mcmaster.magarveylab.prism.data.Cyclization;
import ca.mcmaster.magarveylab.prism.data.Module; 
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.CombinatorialPlan;
import ca.mcmaster.magarveylab.prism.data.sugar.Sugar;
import ca.mcmaster.magarveylab.prism.util.exception.ClassInstantiationException;

/**
 * Generates combinatorial plans.<p>
 * A combinatorial plan is the set of instructions for the construction of a
 * single scaffold within a combinatorial scaffold library. It represents the
 * intersection of a biosynthetic module permutation, a sugar combination, a
 * cyclization pattern, and a set of tailoring reaction plans.
 * 
 * @author skinnider
 *
 */
public class CombinatorialPlanGenerator {

	/**
	 * Generate all combinatorial plans for the biosynthesis of a natural
	 * product scaffold, based on all generated module permutations within a
	 * cluster. A single plan consists of a single module permutation, a single
	 * sugar combination, a single cyclization pattern, and a single unique set
	 * of tailoring reaction plans.
	 * 
	 * @param cluster
	 *            current biosynthetic cluster
	 * @param permutations
	 *            list of module permutations
	 * @return all possible combinatorialization schemes for all module
	 *         permutations
	 * @throws IOException
	 * @throws ClassInstantiationException
	 * @throws CDKException 
	 */
	public static List<CombinatorialPlan> generateAllPlans(Cluster cluster, List<List<Module>> permutations) 
			throws IOException, ClassInstantiationException, CDKException {
		List<CombinatorialPlan> plans = new ArrayList<CombinatorialPlan>();
		CombinatorialData cd = cluster.combinatorialData();
		
		// generate plans 
		int i = 0;
		for (List<Module> permutation : permutations) {
			List<CombinatorialPlan> permutationPlans = generatePlans(cluster, permutation);
			for (CombinatorialPlan permutationPlan : permutationPlans) {
				if (i > 1000) {
					cd.setNumCombinatorialPlans(1001);
					break;					
				}
				plans.add(permutationPlan);
				i++;
			}
		}
		
		if (cd.getNumCombinatorialPlans() != 1001)
			cd.setNumCombinatorialPlans(plans.size());
		System.out.println("[CombinatorialPlanGenerator] Generated "
				+ plans.size() + " combinatorial plans for cluster "
				+ cluster.index());
		return plans;
	}

	/**
	 * Generate a list of possible combinatorialization schemes from a cluster's
	 * module permutation.
	 * 
	 * @param cluster
	 *            current biosynthetic cluster
	 * @param permutation
	 *            generated module permutation
	 * @return all possible combinatorialization schemes for this module
	 *         permutation
	 * @throws IOException
	 * @throws ClassInstantiationException
	 * @throws CDKException 
	 */
	public static List<CombinatorialPlan> generatePlans(Cluster cluster, List<Module> permutation) 
			throws IOException, ClassInstantiationException, CDKException {
		List<CombinatorialPlan> plans = new ArrayList<CombinatorialPlan>();

		List<List<ReactionPlan>> reactionPlans = ReactionPlanGenerator.generateReactionPlans(permutation, cluster);
		List<Cyclization> cyclizations = CyclizationAnalyzer.getAllCyclizations(permutation, cluster);
		List<List<Sugar>> sugars = cluster.sugars();

		// increment number of possible reaction plans 
		CombinatorialData cd = cluster.combinatorialData();
		int reactions = cd.getNumReactions();
		cd.setNumReactions(reactions + reactionPlans.size());

		int i = reactionPlans.size() == 0 ? 1 : reactionPlans.size();
		int j = sugars.size() == 0 ? 1 : sugars.size();
		int k = cyclizations.size() == 0 ? 1 : cyclizations.size();
		
		int numPlans = 0;
		for (int a = 0; a < i; a++) {
			for (int b = 0; b < j; b++) {
				for (int c = 0; c < k; c++) {
					if (numPlans > 1_000)
						break;
					
					List<ReactionPlan> reactionPlan = reactionPlans.size() == 0 ? new ArrayList<ReactionPlan>()
							: reactionPlans.get(a);
					List<Sugar> sugarCombination = sugars.size() == 0 ? new ArrayList<Sugar>()
							: sugars.get(b);

					if (reactionPlan.size() > 0) {
						// copy reaction plans
						List<ReactionPlan> reactionCopy = new ArrayList<ReactionPlan>();
						for (ReactionPlan reaction : reactionPlan)
							reactionCopy.add(new ReactionPlan(reaction));
						reactionPlan = reactionCopy;
					}
					
					// set sugars
					if (sugars.size() > 0 && reactionPlan != null) {
						int x = 0;
						for (ReactionPlan reaction : reactionPlan) { 
							if (reaction.type() == TailoringDomains.GLYCOSYLTRANSFERASE
									|| reaction.type() == TailoringDomains.C_GLYCOSYLTRANSFERASE) {
								Sugar sugar = sugarCombination.get(x);
								reaction.setSmiles(sugar.smiles());
								x++;
							}
						}
					}
					
					Cyclization cyclization = cyclizations.size() == 0 ? null : cyclizations.get(c);
					CombinatorialPlan scheme = new CombinatorialPlan(reactionPlan, cyclization, permutation);
					plans.add(scheme);
					numPlans++;
				}	
			}
		}

		return plans;
	}

}
