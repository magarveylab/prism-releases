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
 * Get all combinations of lysine residues with a glutamate or aspartate
 * residue, for microviridin omega-amide bond formation as catalyzed by the ATP
 * grasp enzyme MdnB.
 * 
 * @author skinnider
 *
 */
public class MdnBAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.MdnB };
	}

	public List<SubstrateSet> findSubstrates(Domain domain,
			List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		List<Module> lys = RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.LYSINE, permutation);

		List<Module> de = new ArrayList<Module>();
		de.addAll(RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.ASPARTATE, permutation));
		de.addAll(RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.GLUTAMATE, permutation));

		if (lys.size() == 0)
			return substrates;
		if (de.size() == 0)
			return substrates;

		List<List<Module>> input = new ArrayList<List<Module>>();
		input.add(lys);
		input.add(de);
		List<List<Module>> output = Combinations.getCombinations(input);

		for (List<Module> combination : output) {
			SubstrateSet substrate = new SubstrateSet(combination);
			substrates.add(substrate);
		}

		return substrates;
	}

}
