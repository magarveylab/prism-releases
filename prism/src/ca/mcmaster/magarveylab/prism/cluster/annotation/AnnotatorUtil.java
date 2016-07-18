package ca.mcmaster.magarveylab.prism.cluster.annotation;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.openscience.cdk.exception.CDKException;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.prism.cluster.reactions.ReactionUtil; 
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.prism.util.exception.ClassInstantiationException;

public class AnnotatorUtil {
	
	/**
	 * Find all potential substrates for a domain associated with a biosynthetic
	 * reaction within a given module permutation.
	 * 
	 * @param domain
	 *            domain to analyze
	 * @param permutation
	 *            module permutation to analyze
	 * @return all potential substrate sets 
	 * @throws IOException
	 * @throws ClassInstantiationException
	 * @throws CDKException
	 */
	public static List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster) 
			throws IOException, ClassInstantiationException, CDKException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		List<Annotator> annotators = ReactionUtil.getAllAnnotators();
		DomainType type = domain.type();
		for (Annotator annotator : annotators)
			if (Arrays.asList(annotator.domains()).contains(type))  {
				substrates.addAll(annotator.findSubstrates(domain, permutation, cluster));
				System.out.println("Got annotator " + annotator.getClass().getName() + " for domain " + type);
			}
		return substrates;
	}

	/**
	 * Convert a list of modules to a list of substrate sets. 	
	 * @param modules	modules to convert
	 * @return			converted modules, as a list of substrate sets (where each set has size = 1)
	 */
	public static List<SubstrateSet> convertModulesToSubstrateSets(List<Module> modules) {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		for (Module module : modules) {
			SubstrateSet substrate = new SubstrateSet(module);
			substrates.add(substrate);
		}
		return substrates;
	}

}
