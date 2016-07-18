package ca.mcmaster.magarveylab.prism.cluster.annotation.ribosomal;

import java.io.IOException;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.analysis.RibosomalClusterAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.AnnotatorUtil;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;

/**
 * Get all cysteine substrates. 
 * 
 * @author skinnider
 *
 */
public class CysteineAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { 
				RibosomalDomains.TpdI,
				RibosomalDomains.TpdL,
				RibosomalDomains.TpdM };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<Module> cys = RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.CYSTEINE, permutation);
		List<SubstrateSet> substrates = AnnotatorUtil.convertModulesToSubstrateSets(cys);
		return substrates;
	}

}
