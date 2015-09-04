package ca.mcmaster.magarveylab.prism.cluster.annotation.impl;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;

/**
 * Find substrates for product template domains in clade I (C2/C7 monocyclic
 * polyketide synthases).
 * 
 * @author skinnider
 *
 */
public class ProductTemplateCladeIAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { ThiotemplatedDomains.PRODUCT_TEMPLATE_I };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		
		if (permutation.size() < 4) {
			System.out.println("Could not find substrates for clade I product template domain: "
					+ "fewer than 4 modules in permutation");
			return substrates;
		}
		
		List<Module> modules = permutation.subList(0, 4);
		SubstrateSet substrate = new SubstrateSet(modules);
		substrates.add(substrate);

		return substrates;
	}

}
