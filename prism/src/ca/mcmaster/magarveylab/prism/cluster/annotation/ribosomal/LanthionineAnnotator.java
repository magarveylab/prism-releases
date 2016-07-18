package ca.mcmaster.magarveylab.prism.cluster.annotation.ribosomal;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

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
 * Find potential combinations of serine/threonine and cysteine domains, at
 * which the LanC cyclase can catalyze lanthionine bond formation. This
 * annotator also detects Dha/Dhb-Cys combinations for GarO, the actagardine
 * S-oxidase.
 * 
 * @author skinnider
 *
 */
public class LanthionineAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.LanC, RibosomalDomains.GarO };
	}

	public List<SubstrateSet> findSubstrates(Domain domain,
			List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		// get all serine/threonine
		List<Module> serthr = new ArrayList<Module>();
		serthr.addAll(RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.SERINE, permutation));
		serthr.addAll(RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.THREONINE, permutation));

		// get all cysteine
		List<Module> cysteines = RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.CYSTEINE, permutation);

		// empty results check
		if (serthr.size() == 0 || cysteines.size() == 0)
			return substrates;

		// remove 1st ser/thr if ElxO is present
		if (cluster.contains(RibosomalDomains.ElxO)) {
			Module m0 = permutation.get(0);
			if (serthr.contains(m0))
				serthr.remove(m0);
		}

		int n = serthr.size();
		int k = cysteines.size();
		if (cluster.contains(RibosomalDomains.LanD))
			k -= 1;
		if (n >= k) {
			// get all permutations of serine/threonine of size [# cysteines]
			List<int[]> serThrPermutations = Permutations.permutations(n, k, 500);
			for (int[] serThrPermutation : serThrPermutations) {
				SubstrateSet substrate = new SubstrateSet();
				for (int i = 0; i < serThrPermutation.length; i++) {
					Module lan = serthr.get(serThrPermutation[i] - 1);
					Module cys = cysteines.get(i);
					substrate.add(lan);
					substrate.add(cys);
				}
				substrates.add(substrate);
			}
		} else {
			// get all permutations of cysteines of size [# ser/thr]
			List<int[]> cysteinePermutations = Permutations.permutations(k, n, 500);
			for (int[] cysteinePermutation : cysteinePermutations) {
				SubstrateSet substrate = new SubstrateSet();
				for (int i = 0; i < cysteinePermutation.length; i++) {
					Module lan = serthr.get(i);
					Module cys = cysteines.get(cysteinePermutation[i] - 1);
					substrate.add(lan);
					substrate.add(cys);
				}
				substrates.add(substrate);
			}
		}

		return substrates;
	}

}
