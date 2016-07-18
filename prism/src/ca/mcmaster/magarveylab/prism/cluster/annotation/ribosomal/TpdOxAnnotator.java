package ca.mcmaster.magarveylab.prism.cluster.annotation.ribosomal;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.analysis.RibosomalClusterAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.cluster.reactions.RibosomalUtil;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;

/**
 * Find potential substrates for thiopeptide cytochrome P450 oxygenases. Because
 * there is no clear relationship between sequence and function, and because the
 * functions of many of these oxygenases has not been experimentally elucidated,
 * multiple possible reactions are considered, including valine
 * beta-hydroxylation, isoleucine epoxidation to hydroxymethylproline,
 * phenylalanine beta-hydroxylation,
 * 
 * @author skinnider
 *
 */
public class TpdOxAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.LanJ };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		// get all subsets 
		List<Module> alanine = RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.ALANINE, permutation);
		List<List<Module>> subsets = RibosomalUtil.getSubsets(alanine);

		// convert to modules
		for (List<Module> subset : subsets) {
			SubstrateSet substrate = new SubstrateSet(subset);
			substrates.add(substrate);
		}

		return substrates;
	}

}
