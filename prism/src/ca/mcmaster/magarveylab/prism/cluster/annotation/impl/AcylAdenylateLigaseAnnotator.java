package ca.mcmaster.magarveylab.prism.cluster.annotation.impl;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.CDKException;

import ca.mcmaster.magarveylab.enums.ClusterFamilies;
import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.clusters.ThiotemplatedClusterTypes;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.enums.substrates.TypeIIPolyketideStarters;
import ca.mcmaster.magarveylab.prism.cluster.analysis.TypeIIPolyketideAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;

/**
 * Find potential sites at which an acyl adenylating ligase could attach a
 * substrate, in cluster families (enediynes and type II polyketides) where
 * these are not integrated in scaffold synthesis.
 * 
 * @author skinnider
 *
 */
public class AcylAdenylateLigaseAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { ThiotemplatedDomains.ACYL_ADENYLATING };
	}

	public List<SubstrateSet> findSubstrates(Domain domain,
			List<Module> permutation, Cluster cluster)
			throws IOException, CDKException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		
		// in type II clusters with a starter unit: acyl-CoA ligase is used to activate starter 
		for (Module module : permutation)
			if (module.type() == ModuleTypes.TYPE_II_PKS_STARTER 
					&& TypeIIPolyketideAnalyzer.getStarter(module) != TypeIIPolyketideStarters.ACETATE)
				return substrates;
		
		// add starter-only ligases which are not in the permutation
		Module module = null;
		for (Module m : permutation)
			if (m.contains(domain)) 
				module = m;
		if (module == null) {
			Annotator annotator = new HydroxylAnnotator();
			substrates.addAll(annotator.findSubstrates(domain, permutation, cluster));
		} else {
			System.out.println("No substrates for domain " + domain.name() + ": in module " + module.scaffold().name());
		}
		
		// get all modules for enediynes and type II PKSs
		if (!cluster.types().contains(ThiotemplatedClusterTypes.ENEDIYNE_9_MEMBERED)
				&& !cluster.types().contains(ThiotemplatedClusterTypes.ENEDIYNE_10_MEMBERED)
				&& !cluster.families().contains(ClusterFamilies.TYPE_II_POLYKETIDE))
			return substrates;
		if (domain.topSubstrate() == null)
			return substrates;

		for (Module m : permutation)
			substrates.add(new SubstrateSet(m));
		return substrates;
	}

}
