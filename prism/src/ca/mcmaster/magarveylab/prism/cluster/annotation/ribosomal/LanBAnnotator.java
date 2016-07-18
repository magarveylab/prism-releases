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
 * Get potential points of reaction for the LanB lantibiotic dehydratase, which
 * converts serine or threonine to dehydroalanine or dehydroaminobutyric acid,
 * respectively. Also used to find substrates for the thiopeptide dehydratase
 * LazB.
 * 
 * @author skinnider
 *
 */
public class LanBAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.LanB, RibosomalDomains.LazB,
				RibosomalDomains.GodF };
	}

	public List<SubstrateSet> findSubstrates(Domain domain,
			List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		// if LanKC or LanM, reject
		if (cluster.contains(RibosomalDomains.LanKC)
				|| cluster.contains(RibosomalDomains.LanM))
			return substrates;
		// goadsporin needs both GodF/GodG
		if (domain.type() == RibosomalDomains.GodF
				&& !cluster.contains(RibosomalDomains.GodG))
			return substrates;

		List<Module> input = new ArrayList<Module>();
		List<Module> ser = RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.SERINE, permutation);
		List<Module> thr = RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.THREONINE, permutation);
		input.addAll(ser);
		input.addAll(thr);

		// get all subsets
		List<List<Module>> subsets = null;
		if (input.size() > 2) {
			subsets = RibosomalUtil.getSubsetsWithMinSize(input, 2);
		} else {
			subsets = RibosomalUtil.getSubsets(input);
		}

		// convert to modules
		for (List<Module> subset : subsets) {
			SubstrateSet substrate = new SubstrateSet(subset);
			substrates.add(substrate);
		}

		return substrates;
	}

}
