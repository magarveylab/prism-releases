package ca.mcmaster.magarveylab.prism.cluster.annotation.ribosomal;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.analysis.RibosomalClusterAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;
import ca.mcmaster.magarveylab.prism.util.exception.RibosomalScaffoldGenerationException;

/**
 * Get all potential sites of sulfhydryl-alpha carbon bond formation in
 * sactipeptides, where the molecule must form a hairpin, at least two residues
 * must separate each cysteine-alpha carbon bond, and the final cysteine-alpha
 * carbon bond must be at least two residues from the end.
 * 
 * @author skinnider
 *
 */
public class SactipeptideAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.AlbA };
	}

	public List<SubstrateSet> findSubstrates(Domain domain,
			List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		List<Module> cysteines = RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.CYSTEINE, permutation);
		if (cysteines.size() == 0) {
			System.out.println("Error: could not find any cysteine "
					+ "residues in this permutation!");
			return substrates;
		}

		List<List<Module>> results = new ArrayList<List<Module>>();
		if (cluster.contains(RibosomalDomains.SkfH)) {
			// if there is a disulfide bond, consider only combinations w/ 1
			// cysteine
			for (Module cysteine : cysteines) {
				int idx = permutation.indexOf(cysteine);
				for (int i = idx + 6; i < permutation.size(); i++) {
					Module module = permutation.get(i);
					if (cysteines.indexOf(module) != -1)
						continue;
					List<Module> modules = new ArrayList<Module>();
					modules.add(cysteine);
					modules.add(module);
					results.add(modules);
				}
			}
		} else {
			// if there is no disulfide bond, consider combinations w/ 3
			// cysteines
			int lastCysteineIdx = permutation.indexOf(cysteines.get(cysteines
					.size() - 1));
			annotate(permutation, results, lastCysteineIdx + 3, 0,
					new ArrayList<Module>());
		}

		for (List<Module> result : results) {
			SubstrateSet substrate = new SubstrateSet(result);
			substrates.add(substrate);
		}

		System.out.println("Found " + substrates.size()
				+ " substrate combinations: ");
		for (SubstrateSet substrate : substrates) {
			StringBuffer sb = new StringBuffer();
			for (Module module : substrate.getAllModules())
				sb.append(module.scaffold().name().replace("_", "") + " ");
			System.out.println("Combination: " + sb.toString());
		}

		return substrates;
	}

	public static void annotate(List<Module> permutation,
			List<List<Module>> results, int start, int depth,
			List<Module> current) {
		List<Module> cysteines = RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.CYSTEINE, permutation);

		if (depth == cysteines.size()) {
			List<Module> result = new ArrayList<Module>();
			for (int i = 0; i < current.size(); i++)
				result.add(current.get(i));
			results.add(result);
			current.clear();
			return;
		} else {
			if (permutation.size() - 5 - start <= 0)
				current.clear();

			for (int i = start + 3; i < permutation.size() - 2; i++) {
				if (results.size() > 0 && current.size() == 0) {
					int lastIdx = results.size() - 1;
					List<Module> last = results.get(lastIdx);
					for (int j = 0; j < depth * 2; j++)
						current.add(last.get(j));
				}

				Module previousCysteine = cysteines.get(cysteines.size()
						- (depth + 1));
				current.add(previousCysteine);
				Module alphaCarbonModule = permutation.get(i);
				current.add(alphaCarbonModule);

				annotate(permutation, results, i, depth + 1, current);
			}
		}
	}

	public static void main(String[] args) throws InvalidSmilesException,
			IOException, RibosomalScaffoldGenerationException {
		String sactipeptide = "CMGCWASKSIAMTRVCALPHPAMRAI";
		List<List<Module>> lists = RibosomalClusterAnalyzer
				.generateRibosomalModulePermutations(sactipeptide, null);
		List<Module> permutation = lists.get(0);

		List<SubstrateSet> substrates = new SactipeptideAnnotator()
				.findSubstrates(null, permutation, null);
		System.out.println("Found " + substrates.size()
				+ " substrate combinations: ");
		for (SubstrateSet substrate : substrates) {
			StringBuffer sb = new StringBuffer();
			for (Module module : substrate.getAllModules())
				sb.append(module.scaffold().name().replace("_", "") + " ");
			System.out.println("Combination: " + sb.toString());
		}
	}

}
