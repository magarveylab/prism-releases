package ca.mcmaster.magarveylab.prism.cluster.annotation.impl;

import java.io.IOException;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.domains.BetaLactamDomains;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.substrates.AdenylationSubstrates;
import ca.mcmaster.magarveylab.prism.cluster.analysis.ModuleAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.AnnotatorUtil;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;

/**
 * Find sites at which the isopenicillin N acyltransferase can react (i.e.,
 * 2-amino-adipic acid residues).
 * 
 * @author skinnider
 *
 */
public class IsopenicillinNAcyltransferaseAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { BetaLactamDomains.IAT };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster) 
			throws InvalidSmilesException, IOException {
		List<Module> modules = ModuleAnalyzer.modules(permutation, AdenylationSubstrates._2_AMINO_ADIPIC_ACID);
		List<SubstrateSet> substrates = AnnotatorUtil.convertModulesToSubstrateSets(modules);
		return substrates;
	}

}
