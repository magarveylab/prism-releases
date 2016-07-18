package ca.mcmaster.magarveylab.prism.cluster.annotation.ribosomal;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.analysis.RibosomalClusterAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.combinatorialization.Combinations;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;

/**
 * Get potential substrates for lysine-tryptophan crosslinking in streptide
 * clusters.
 * 
 * @author skinnider
 *
 */
public class StrBAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.StrB };
	}

	public List<SubstrateSet> findSubstrates(Domain domain,
			List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		List<Module> k = RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.LYSINE, permutation);
		List<Module> w = RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.TRYPTOPHAN, permutation);

		if (k.size() == 0 || w.size() == 0)
			return substrates;

		List<List<Module>> input = new ArrayList<List<Module>>();
		input.add(k);
		input.add(w);
		List<List<Module>> output = Combinations.getCombinations(input);
		
		for (List<Module> combination : output) {
			int idx1 = permutation.indexOf(combination.get(0));
			int idx2 = permutation.indexOf(combination.get(1));
			if (Math.abs(idx1 - idx2) > 1) {
				SubstrateSet substrate = new SubstrateSet(combination);
				substrates.add(substrate);
			}
		}

		return substrates;
	}

}
