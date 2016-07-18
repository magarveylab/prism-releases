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
 * All known glycocins contain two disulfide bonds, so get all combinations of
 * four cysteines to form two disulfide bonds.
 * 
 * @author skinnider
 *
 */
public class SunAAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.SunA };
	}

	public List<SubstrateSet> findSubstrates(Domain domain,
			List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		List<Module> cys = RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.CYSTEINE, permutation);

		if (cys.size() < 2)
			return substrates;
		
		List<int[]> permutations = Permutations.permutations(cys.size(), 4, 500);
		for (int[] p : permutations) {
			// need to form hairpins
			if (p[0] > p[2] || p[0] > p[1] || p[2] > p[3] || p[3] > p[1])
				continue;

			List<Module> modules = new ArrayList<Module>();
			for (int i : p)
				modules.add(cys.get(i-1));
			
			SubstrateSet substrate = new SubstrateSet(modules);
			substrates.add(substrate);
		}

		return substrates;
	}

}
