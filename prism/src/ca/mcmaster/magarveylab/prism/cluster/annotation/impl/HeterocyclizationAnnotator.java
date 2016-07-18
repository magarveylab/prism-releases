package ca.mcmaster.magarveylab.prism.cluster.annotation.impl;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.prism.cluster.analysis.DomainAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.prism.enums.hmms.AdenylationHmms;

public class HeterocyclizationAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { ThiotemplatedDomains.CONDENSATION };
	}

	public List<SubstrateSet> findSubstrates(Domain domain,
			List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		if (domain.type() == ThiotemplatedDomains.CONDENSATION) {
			// condensation domains which are not cyclization domains can be ignored
			if (!DomainAnalyzer.isCyclization(domain))
				return substrates;
			// cyclization domains which do not surround Ser or Thr can be ignored
			for (Module module : permutation)
				if (module.contains(domain))
					if (module.scaffold() != null 
							&& module.scaffold().topSubstrate().type() != AdenylationHmms.THREONINE_1
							&& module.scaffold().topSubstrate().type() != AdenylationHmms.THREONINE_2
							&& module.scaffold().topSubstrate().type() != AdenylationHmms.CYSTEINE_1
							&& module.scaffold().topSubstrate().type() != AdenylationHmms.CYSTEINE_2)
						return substrates;
		}

		Module module = null;
		for (Module m : permutation)
			if (m.contains(domain))
				module = m;
		
		if (module != null) {
			int idx = permutation.indexOf(module);
			if (idx > 0) {
				Module last = permutation.get(idx-1);
				SubstrateSet substrate = new SubstrateSet(last, module);
				substrates.add(substrate);
			}
		}

		return substrates;
	}

}
