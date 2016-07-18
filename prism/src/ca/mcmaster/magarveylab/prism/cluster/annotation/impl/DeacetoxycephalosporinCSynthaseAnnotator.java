package ca.mcmaster.magarveylab.prism.cluster.annotation.impl;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.domains.BetaLactamDomains;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.prism.enums.hmms.AdenylationHmms;

public class DeacetoxycephalosporinCSynthaseAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { BetaLactamDomains.DOACS };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster) 
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		for (int i = 0; i < permutation.size() - 1; i++) {
			Module module = permutation.get(i);
			Module next = permutation.get(i+1);
			if (module.isAdenylationModule() && next.isAdenylationModule()) {
				Domain d1 = module.scaffold();
				Domain d2 = next.scaffold();
				if ((d1.topSubstrate().type() == AdenylationHmms.CYSTEINE_1 
						|| d1.topSubstrate().type() == AdenylationHmms.CYSTEINE_2) 
						&& (d2.topSubstrate().type() == AdenylationHmms.VALINE_1
						|| d2.topSubstrate().type() == AdenylationHmms.VALINE_2
						|| d2.topSubstrate().type() == AdenylationHmms.VALINE_3)) {
					Module cysteine = module;
					Module valine = next;
					SubstrateSet substrate = new SubstrateSet(cysteine, valine);
					substrates.add(substrate);
				}
			}
		}
		return substrates;
	}
	
}
