package ca.mcmaster.magarveylab.prism.cluster.annotation.ribosomal;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;
import org.paukov.combinatorics.Factory;
import org.paukov.combinatorics.Generator;
import org.paukov.combinatorics.ICombinatoricsVector;

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
 * Get all combinations of two serine/threonine residues with two
 * glutamate/aspartate residues, for microviridin and marinostatin omega-ester
 * bond formation as catalyzed by the ATP grasp enzyme MdnC.
 * 
 * @author skinnider
 *
 */
public class MdnCAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.MdnC };
	}

	public List<SubstrateSet> findSubstrates(Domain domain,
			List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		List<Module> st = new ArrayList<Module>();
		st.addAll(RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.SERINE, permutation));
		st.addAll(RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.THREONINE, permutation));

		List<Module> de = new ArrayList<Module>();
		de.addAll(RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.ASPARTATE, permutation));
		de.addAll(RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.GLUTAMATE, permutation));

		if (st.size() == 0)
			return substrates;
		if (de.size() == 0)
			return substrates;

		List<List<Module>> input = new ArrayList<List<Module>>();
		input.add(st);
		input.add(de);
		List<List<Module>> output = Combinations.getCombinations(input);

		ICombinatoricsVector<List<Module>> vector = Factory
				.createVector(output);
		Generator<List<Module>> generator = Factory
				.createSimpleCombinationGenerator(vector, 2);
		for (ICombinatoricsVector<List<Module>> subset : generator) {
			List<List<Module>> combination = subset.getVector();
			SubstrateSet substrate = new SubstrateSet();
			for (List<Module> pair : combination)
				substrate.addAll(pair);
			substrates.add(substrate);
		}

		return substrates;
	}

}
