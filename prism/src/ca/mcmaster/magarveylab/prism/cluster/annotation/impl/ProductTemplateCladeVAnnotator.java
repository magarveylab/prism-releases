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
 * Find substrates for product template domains in clade V (C6/C11 and C4/C13
 * cyclization, as in the emodin and asperthecin clusters).
 * 
 * @author skinnider
 *
 */
public class ProductTemplateCladeVAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { ThiotemplatedDomains.PRODUCT_TEMPLATE_V };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		
		if (permutation.size() < 7) {
			System.out.println("Could not find substrates for clade V product template domain: "
					+ "fewer than 7 modules in permutation");
			return substrates;
		}
		
		List<Module> modules = permutation.subList(1, 7);
		SubstrateSet substrate = new SubstrateSet(modules);
		substrates.add(substrate);

		return substrates;
	}

}
