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
 * Get potential phenylalanine substrates for the YM-216391 phenylalanine
 * beta-hydroxylase YmE, thiopeptide phenylalanine beta-hydroxylases, and the
 * bottromycin beta-C-methyltransferase BotRMT1.
 * 
 * @author skinnider
 *
 */
public class PhenylalanineAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.YmE,
				RibosomalDomains.TpdQ,
				RibosomalDomains.BotRMT1 };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<Module> phe = RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.PHENYLALANINE, permutation);
		List<SubstrateSet> substrates = AnnotatorUtil.convertModulesToSubstrateSets(phe);
		return substrates;
	}

}
