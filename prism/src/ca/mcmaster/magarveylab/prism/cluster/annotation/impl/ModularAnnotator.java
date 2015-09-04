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
 * Find substrates for enzymatic domains for which the substrate is located in the same module as the domain itself. 
 * @author skinnider
 *
 */
public class ModularAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { ThiotemplatedDomains.N_METHYLTRANSFERASE,
				ThiotemplatedDomains.C_METHYLTRANSFERASE,
				ThiotemplatedDomains.O_METHYLTRANSFERASE,
				ThiotemplatedDomains.NITROREDUCTASE };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		Module module = null;
		for (Module m : permutation)
			if (m.contains(domain))
				module = m;

		if (module != null) {
			SubstrateSet substrate = new SubstrateSet(module);
			substrates.add(substrate);
		} 

		return substrates;
	}

}
