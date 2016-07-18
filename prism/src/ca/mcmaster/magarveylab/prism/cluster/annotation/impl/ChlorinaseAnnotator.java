package ca.mcmaster.magarveylab.prism.cluster.annotation.impl;

import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.prism.blast.BlastSearchResult;
import ca.mcmaster.magarveylab.prism.cluster.analysis.ModuleAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.AnnotatorUtil;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.prism.enums.hmms.AdenylationHmms;

public class ChlorinaseAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { TailoringDomains.CHLORINATION };
	}

	/**
	 * Find potential substrates for a chlorinase domain given its BLAST search results.
	 * @param domain		chlorinase domain to analyze
	 * @param permutation	permutation of moduels being analyzed
	 * @return				all potential substrates for this chlorination reaction
	 */
	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster) {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		List<Module> modules = new ArrayList<Module>();
		
		for (BlastSearchResult hit : domain.blastResults()) {
			String name = hit.subject().toLowerCase();
			System.out.println("[ChlorinaseAnnotator] Finding substrates for " + name);
			
			if (name.contains("tyrosine") || name.contains("phenylglycine")) {
				modules = ModuleAnalyzer.sixMemberedAromatic(permutation);
			} else if (name.contains("tryptophan")) {
				modules = ModuleAnalyzer.modules(permutation, AdenylationHmms.TRYPTOPHAN);
			} else if (name.contains("leucine")) {
				modules = ModuleAnalyzer.bcaa(permutation);
			} else if (name.contains("threonine")) {
				modules.addAll(ModuleAnalyzer.modules(permutation, AdenylationHmms.THREONINE_1));
				modules.addAll(ModuleAnalyzer.modules(permutation, AdenylationHmms.THREONINE_2));
				modules.addAll(ModuleAnalyzer.modules(permutation, AdenylationHmms.ALLO_THREONINE));
			} else if (name.contains("histidine")) {
				modules = ModuleAnalyzer.modules(permutation, AdenylationHmms.HISTIDINE_1);
				modules.addAll(ModuleAnalyzer.modules(permutation, AdenylationHmms.HISTIDINE_2));
			} else if (name.contains("c12") || name.contains("c14")) {
				for (Module module : permutation)
					if (module.type() == ModuleTypes.ACYLTRANSFERASE 
							|| module.type() == ModuleTypes.TYPE_II_PKS)
						modules.add(module);
			} else if (name.contains("starter")) {
				if (permutation.size() > 0) {
					Module substrate = permutation.get(0);
					modules.add(substrate);
				}
			} else if (name.contains("proline")) { 
				if (permutation.size() > 0) {
					modules.addAll(ModuleAnalyzer.modules(permutation, AdenylationHmms.PROLINE_1));
					modules.addAll(ModuleAnalyzer.modules(permutation, AdenylationHmms.PROLINE_2));
					modules.addAll(ModuleAnalyzer.modules(permutation, AdenylationHmms.PROLINE_3));
					modules.addAll(ModuleAnalyzer.modules(permutation, AdenylationHmms.METHYL_PROLINE));
					Module substrate = permutation.get(0);
					modules.add(substrate);
				}
			} else {
				continue;
			}
			
			substrates = AnnotatorUtil.convertModulesToSubstrateSets(modules);
			if (substrates.size() > 0)
				break;
		}
		
		System.out.println("[ChlorinaseAnnotator] Found " + substrates.size() + " substrates for " + domain.type() 
				+ " domain");
		return substrates;
	}
	
}
