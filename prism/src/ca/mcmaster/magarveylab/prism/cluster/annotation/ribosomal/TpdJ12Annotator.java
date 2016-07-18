package ca.mcmaster.magarveylab.prism.cluster.annotation.ribosomal;

import java.io.IOException;
import java.util.ArrayList;
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
 * Get potential substrates for the thiomuracin biosynthesis enzymes TpdJ1 and
 * TpdJ2, whose function is unclear but which perform phenylalanine
 * beta-hydroxylation and isoleucine dihydroxylation
 * 
 * @author skinnider
 *
 */
public class TpdJ12Annotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.TpdJ1, 
				RibosomalDomains.TpdJ2 };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<Module> ile = RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.ISOLEUCINE, permutation);
		List<Module> phe = RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.PHENYLALANINE, permutation);
		
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>(); 
		
		substrates.addAll(AnnotatorUtil.convertModulesToSubstrateSets(ile));
		substrates.addAll(AnnotatorUtil.convertModulesToSubstrateSets(phe));
		return substrates;
	}

}
