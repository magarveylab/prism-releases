package ca.mcmaster.magarveylab.prism.cluster.reactions;

import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ca.mcmaster.magarveylab.enums.ReactionPriorities;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.enums.Reactions;
import ca.mcmaster.magarveylab.prism.util.exception.ClassInstantiationException;

/**
 * Utility functions involving reactions.
 * 
 * @author skinnider
 *
 */
public class ReactionUtil {

	/**
	 * Get the reaction for a reaction plan.
	 * 
	 * @param plan
	 *            reaction plan to analyze
	 * @param scaffold
	 *            scaffold on which to execute the reaction
	 * @return
	 * @throws ClassInstantiationException
	 */
	public static Reaction getReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster)
			throws ClassInstantiationException {
		Domain domain = plan.domain();
		DomainType type = domain.type();
		Reaction reaction = null;

		for (Reactions reactionType : Reactions.values()) {
			Reaction r;
			try {
				Constructor<? extends GenericReaction> c = reactionType.reactionClass()
						.getConstructor(ReactionPlan.class, Scaffold.class, Cluster.class);
				r = c.newInstance(plan, scaffold, cluster);
				if (Arrays.asList(r.domains()).contains(type))
					reaction = r;
			} catch (Exception e) {
				throw new ClassInstantiationException("Could not instantiate reaction with type " + reactionType);
			}
		}

		if (reaction == null)
			throw new ClassInstantiationException("Could not get reaction for tailoring plan with domain type "
					+ type.fullName());

		return reaction;
	}
	
	/**
	 * Get the priority of the reaction associated with a domain. 
	 * @param domain	domain to analyze
	 * @return			the priority 
	 */
	public static ReactionPriorities getPriority(Domain domain) {
		ReactionPriorities priority = null;
		DomainType type = domain.type();
		for (ReactionPriorities rp : ReactionPriorities.values())
			if (rp.domain() == type)
				priority = rp;
		return priority;
	}

	/**
	 * Get all unique annotators for all reactions (i.e. redundant annotators
	 * will be removed).
	 * 
	 * @return all annotators
	 * @throws ClassInstantiationException
	 */
	public static List<Annotator> getAllAnnotators()
			throws ClassInstantiationException {
		List<Annotator> all = new ArrayList<Annotator>();
		for (Reactions reactionType : Reactions.values()) {
			Annotator[] annotators = reactionType.annotators();
			annotatorLoop: for (Annotator annotator : annotators) {
				for (Annotator a : all)
					if (a.getClass().equals(annotator.getClass()))
						continue annotatorLoop;
				all.add(annotator);
			}
		}
		return all;
	}

}
