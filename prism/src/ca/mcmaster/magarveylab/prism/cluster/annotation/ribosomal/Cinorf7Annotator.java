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
 * Find potential sites at which homologs of the cinnamycin orf 7 can catalyze
 * lysinoalanine formation.
 * 
 * @author skinnider
 *
 */
public class Cinorf7Annotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.Cinorf7 };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		
		// get all serine
		List<Module> serine = RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.SERINE, permutation);

		// get all lysine
		List<Module> lysine = RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.LYSINE, permutation);

		// get all combinations
		List<List<Module>> input = new ArrayList<List<Module>>();
		input.add(serine);
		input.add(lysine);
		List<List<Module>> output = Combinations.getCombinations(input);
		
		for (List<Module> modules : output) {
			SubstrateSet substrate = new SubstrateSet();
			substrate.add(modules.get(0)); // serine
			substrate.add(modules.get(1)); // lysine
			substrates.add(substrate);
		}

		return substrates;
	}

}
