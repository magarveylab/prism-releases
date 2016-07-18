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
import ca.mcmaster.magarveylab.prism.combinatorialization.Permutations;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;

/**
 * Get all potential sites of disulfide bond formation in lassopeptides (where
 * 1, 2, or 3 disulfide bonds are considered).
 * 
 * @author skinnider
 *
 */
public class LassoPeptideDisulfideAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.Lasso_precursors };
	}

	public List<SubstrateSet> findSubstrates(Domain domain,
			List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		List<Module> cys = RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.CYSTEINE, permutation);

		List<List<Module>> subsets = new ArrayList<List<Module>>();

		ICombinatoricsVector<Module> vector = Factory.createVector(cys);
		Generator<Module> gen0 = Factory.createSimpleCombinationGenerator(
				vector, 2);
		Generator<Module> gen1 = Factory.createSimpleCombinationGenerator(
				vector, 4);
		Generator<Module> gen2 = Factory.createSimpleCombinationGenerator(
				vector, 6);
		for (ICombinatoricsVector<Module> subset : gen0)
			subsets.add(subset.getVector());
		for (ICombinatoricsVector<Module> subset : gen1) {
			// get permutations
			List<Module> combination = subset.getVector();
			List<SubstrateSet> permutations = getPermutations(combination);
			for (SubstrateSet s : permutations)
				substrates.add(s);
		}
		for (ICombinatoricsVector<Module> subset : gen2) {
			// get permutations
			List<Module> combination = subset.getVector();
			List<SubstrateSet> permutations = getPermutations(combination);
			for (SubstrateSet s : permutations)
				substrates.add(s);
		}

		for (List<Module> subset : subsets) {
			SubstrateSet substrate = new SubstrateSet(subset);
			substrates.add(substrate);
		}

		return substrates;
	}

	public static List<SubstrateSet> getPermutations(List<Module> combination) {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		List<int[]> permutations = Permutations.permutations(
				combination.size(), combination.size(), 500);
		for (int[] p : permutations) {
			List<Module> modules = new ArrayList<Module>();
			for (int i = 0; i < p.length; i++)
				modules.add(combination.get(p[i] - 1));
			substrates.add(new SubstrateSet(modules));
		}
		return substrates;
	}

}