package ca.mcmaster.magarveylab.prism.cluster.annotation.impl;

import java.io.IOException;
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

public class TryptophanDioxygenaseAnnotator implements Annotator {
	
	public DomainType[] domains() {
		return new DomainType[] { TailoringDomains.TRYPTOPHAN_DIOXYGENASE };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster) 
			throws InvalidSmilesException, IOException {
		List<Module> modules = ModuleAnalyzer.modules(permutation, AdenylationHmms.TRYPTOPHAN);
		List<SubstrateSet> substrates = AnnotatorUtil.convertModulesToSubstrateSets(modules);
		return substrates;
	}
	
}
