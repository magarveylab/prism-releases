package ca.mcmaster.magarveylab.prism.cluster.annotation.impl;

import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.Substrate;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.prism.enums.hmms.AdenylationHmms;

public class P450CAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { TailoringDomains.P450C };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster) {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		if (permutation.size() < 7)
			return substrates;
		
		int x = 0;
		Module first = null, second = null;
		for (Module module : permutation) 
			if (module.type() == ModuleTypes.ADENYLATION) {
				if (x == 4)
					first = module;
				if (x == 6)
					second = module;
				x++;
			}
		if (first == null || second == null)
			return substrates;

		Domain firstScaffold = first.scaffold();
		Domain secondScaffold = second.scaffold();
		if (firstScaffold == null || secondScaffold == null)
			return substrates;
		Substrate firstSubstrate = firstScaffold.topSubstrate();
		if (firstSubstrate == null || !(firstSubstrate.type() == AdenylationHmms.TYROSINE_1 
				|| firstSubstrate.type() == AdenylationHmms.TYROSINE_2 
				|| firstSubstrate.type() == AdenylationHmms.BETA_HYDROXYTYROSINE
				|| firstSubstrate.type() == AdenylationHmms._4_HYDROXY_PHENYLGLYCINE))
			return substrates;
		Substrate secondSubstrate = secondScaffold.topSubstrate();
		if (secondSubstrate == null || secondSubstrate.type() != AdenylationHmms._3_5_DIHYDROXYPHENYLGLYCINE)
			return substrates;

		SubstrateSet substrate = new SubstrateSet(first, second);
		substrates.add(substrate); 
		return substrates;
	}
	
}
