package ca.mcmaster.magarveylab.prism.cluster.annotation.impl;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.prism.cluster.analysis.ModuleAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.AnnotatorUtil;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.prism.enums.hmms.AdenylationHmms;

public class ProlineDehydrogenaseAnnotator implements Annotator {
	
	public DomainType[] domains() {
		return new DomainType[] { TailoringDomains.PROLINE_DEHYDROGENASE };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster) 
			throws InvalidSmilesException, IOException {
		List<Module> modules = new ArrayList<Module>();
		modules.addAll(ModuleAnalyzer.modules(permutation, AdenylationHmms.PROLINE_1));
		modules.addAll(ModuleAnalyzer.modules(permutation, AdenylationHmms.PROLINE_2));
		modules.addAll(ModuleAnalyzer.modules(permutation, AdenylationHmms.PROLINE_3));
		modules.addAll(ModuleAnalyzer.modules(permutation, AdenylationHmms.METHYL_PROLINE));
		List<SubstrateSet> substrates = AnnotatorUtil.convertModulesToSubstrateSets(modules);
		return substrates;
	}
	
}
