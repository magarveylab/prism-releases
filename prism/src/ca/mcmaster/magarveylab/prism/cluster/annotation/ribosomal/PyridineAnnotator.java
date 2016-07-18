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
 * Get potential points of reaction for the thiopeptide cycloaddition enzyme
 * implicated in formation of the core pyridine ring, defined here as
 * combinations of two serine residues (plus the residue immediately before the
 * second serine).<br>
 * <br>
 * When the hit is to the LazC model, and the first residue is a serine, then
 * only pyridines (not tetrahydropyridines) are formed. Conversely, when the hit
 * is to the tetrahydropyridine-forming LazC_b model, and there are more than
 * two serines, then only tetrahydropyridines are formed.
 * 
 * @author skinnider
 *
 */
public class PyridineAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.LazC,
				RibosomalDomains.LazC_b, RibosomalDomains.NosC };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		
		if (domain.type() == RibosomalDomains.LazC
				&& cluster.contains(RibosomalDomains.LazC_b))
			return substrates;
		
		if (permutation.size() == 0)
			return substrates;
		Module m0 = permutation.get(0);
		
		List<Module> serine = RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.SERINE, permutation);
		List<List<Module>> subsets = RibosomalUtil.getSubsets(serine, 2);
		for (List<Module> subset : subsets) {
			Module first = subset.get(0);
			Module second = subset.get(1);
			int idx1 = permutation.indexOf(first);
			int idx2 = permutation.indexOf(second);
			if (idx2 - idx1 <= 7 || idx2 - idx1 >= 14)
				// thiopeptide macrocycles are 9-12 aa, so only consider 8-13 aa 
				continue;
			if (idx2 == 0)
				// need residue (n-1) to form 6-membered ring
				continue;
			if (domain.type() == RibosomalDomains.LazC
					&& serine.contains(m0) && idx1 > 0)
				// must form pyridine (not tetrahydropyridine) when LazC and
				// first residue serine
				continue;
			if (domain.type() == RibosomalDomains.LazC_b && idx1 == 0
					&& serine.size() > 2)
				// must form tetrahydropyridine when LazC and >2 serines
				continue;
			
			Module third = permutation.get(idx2 - 1);

			SubstrateSet substrate = new SubstrateSet();
			substrate.addAll(subset);
			substrate.add(third);
			substrates.add(substrate);
		}
		
		System.out.println("Found " + substrates.size() + " LazC substrates");
		
		return substrates;
	}

}
