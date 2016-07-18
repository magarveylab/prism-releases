package ca.mcmaster.magarveylab.prism.cluster.analysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.CDKException;

import ca.mcmaster.magarveylab.enums.CyclizationPatterns;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.CombinatorialData;
import ca.mcmaster.magarveylab.prism.data.Cyclization;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.enums.hmms.AdenylationHmms;

/**
 * Performs analysis of natural product macrocyclization patterns.
 * 
 * @author skinnider
 *
 */
public class CyclizationAnalyzer {

	/**
	 * Get all potential cyclization patterns for this putative natural product
	 * scaffold.
	 * 
	 * @param permutation
	 *            module pattern to analyze
	 * @param cluster
	 *            parent cluster
	 * @return all potential cyclization patterns
	 * @throws IOException
	 * @throws CDKException 
	 */
	public static List<Cyclization> getAllCyclizations(List<Module> permutation, Cluster cluster) 
			throws IOException, CDKException {
		List<Cyclization> cyclizations = new ArrayList<Cyclization>();
		
		// identify potential cyclizations 
		cyclizations.addAll(getMacrolactamCyclizations(permutation, cluster));
		cyclizations.addAll(getMacrolactoneCyclizations(permutation, cluster));
		cyclizations.addAll(getLinearCyclizations(permutation, cluster));
		cyclizations.addAll(getImineCyclizations(permutation, cluster));
		
		// set # of cyclizations 
		CombinatorialData cd = cluster.combinatorialData();
		cd.setNumCyclizations(cyclizations.size());
		
		return cyclizations;
	}

	/**
	 * Get a list of potential macrolactam cyclizations for a natural product
	 * scaffold, defined here as the first module in a permutation when it is an
	 * adenylation (or acyl-adenylating) module.
	 * 
	 * @param permutation
	 *            permutation to analyze
	 * @return potential macrolactam cyclization(s)
	 */
	public static List<Cyclization> getMacrolactamCyclizations(List<Module> permutation, Cluster cluster) {
		List<Cyclization> cyclizations = new ArrayList<Cyclization>();
		
		if (cluster.contains(ThiotemplatedDomains.REDUCTASE)
				|| (cluster.contains(TailoringDomains.P450A)
						&& cluster.contains(TailoringDomains.P450B) 
						&& cluster.contains(TailoringDomains.P450C)))
			return cyclizations;

		// get first module (if A domain)
		if (permutation.size() > 3) {
			Module module = permutation.get(0);
			if (module.isAdenylationModule()) {
				Cyclization c = new Cyclization(module, CyclizationPatterns.LACTAM);
				cyclizations.add(c);
			}
		}

		return cyclizations;
	}
	
	/**
	 * Get a list of potential macrolactone cyclizations within this putative
	 * natural product, defined as serines, threonines, fatty acids with free
	 * hydroxyls, and beta-hydroxylated amino acids.
	 * 
	 * @param permutation
	 *            permutation to analyze
	 * @return potential macrolactone cyclization(s)
	 * @throws IOException
	 * @throws CDKException 
	 */
	public static List<Cyclization> getMacrolactoneCyclizations(List<Module> permutation, Cluster cluster) 
			throws IOException, CDKException {
		List<Cyclization> cyclizations = new ArrayList<Cyclization>();
		
		if (cluster.contains(ThiotemplatedDomains.REDUCTASE)
				|| (cluster.contains(TailoringDomains.P450A)
						&& cluster.contains(TailoringDomains.P450B) 
						&& cluster.contains(TailoringDomains.P450C)))
			return cyclizations;
		
		List<Module> modules = new ArrayList<Module>();
		
		modules.addAll(ModuleAnalyzer.cyclizationHydroxyls(permutation));
		// get beta-hydroxylated AAs
		modules.addAll(ModuleAnalyzer.modules(permutation, AdenylationHmms.BETA_HYDROXY_PHENYLALANINE));
		modules.addAll(ModuleAnalyzer.modules(permutation, AdenylationHmms.BETA_HYDROXY_ASPARAGINE));
		modules.addAll(ModuleAnalyzer.modules(permutation, AdenylationHmms.BETA_HYDROXY_ASPARTATE));
		modules.addAll(ModuleAnalyzer.modules(permutation, AdenylationHmms.BETA_HYDROXY_LEUCINE));

		for (Module module : modules) {
			Cyclization c = new Cyclization(module, CyclizationPatterns.LACTONE);
			cyclizations.add(c);
		}
		
		return cyclizations;
	}
	
	/**
	 * Get a list of potential cyclization patterns which represent the absence
	 * of a cyclization, i.e. a linear natural product.
	 * 
	 * @param permutation
	 *            permutation to analyze
	 * @param cluster
	 *            potential linear cyclization(s)
	 * @return
	 */
	public static List<Cyclization> getLinearCyclizations(List<Module> permutation, Cluster cluster) {
		List<Cyclization> cyclizations = new ArrayList<Cyclization>();
		
		Cyclization cyclization = null;
		if (cluster.contains(ThiotemplatedDomains.REDUCTASE)) {
			cyclization = new Cyclization(null, CyclizationPatterns.LINEAR_ALDEHYDE);
		} else {
			cyclization = new Cyclization(null, CyclizationPatterns.LINEAR);
		}
		cyclizations.add(cyclization);
		
		return cyclizations;
	}

	/**
	 * Get a list of potential imine cyclization patterns within this putative
	 * natural product--only executed if the cluster contains a reductase
	 * domain.
	 * 
	 * @param permutation
	 *            permutation to analyze
	 * @param cluster
	 *            potential cyclic imine cyclization patterns
	 * @return
	 */
	public static List<Cyclization> getImineCyclizations(List<Module> permutation,
			Cluster cluster) {
		List<Cyclization> cyclizations = new ArrayList<Cyclization>();

		if (cluster.contains(ThiotemplatedDomains.REDUCTASE)
				&& permutation.size() > 0) {
			Module first = permutation.get(0);
			if (first.isAdenylationModule()) {
				Cyclization cyclization = new Cyclization(first, CyclizationPatterns.IMINE);
				cyclizations.add(cyclization);
			}
		}
		
		return cyclizations;
	}

}
