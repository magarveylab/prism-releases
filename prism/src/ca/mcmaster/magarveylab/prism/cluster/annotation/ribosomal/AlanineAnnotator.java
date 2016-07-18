package ca.mcmaster.magarveylab.prism.cluster.annotation.ribosomal;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;
import ca.mcmaster.magarveylab.prism.cluster.analysis.RibosomalClusterAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.cluster.reactions.RibosomalUtil;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;

/**
 * Get potential alanine substrates for the lantipeptide oxidative enzyme which
 * converts dehydroalanine to D-alanine.
 * 
 * @author skinnider
 *
 */
public class AlanineAnnotator implements Annotator {

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
