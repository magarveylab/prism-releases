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
 * Get potential substrates for modified indolic acid attachment in the
 * nosiheptide and nocathiacin clusters. The ketone is assumed to be esterified
 * or thioesterified at free serines and cysteines, while the oxygen is assumed
 * to react at glutamate or aspartate.
 * 
 * @author skinnider
 *
 */
public class NosIAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.NosI };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		
		List<Module> cs = new ArrayList<Module>();
		cs.addAll(RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.CYSTEINE, permutation));
		cs.addAll(RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.SERINE, permutation));
		
		List<Module> de = new ArrayList<Module>();
		de.addAll(RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.GLUTAMATE, permutation));
		de.addAll(RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.ASPARTATE, permutation));
		
		if (de.size() == 0 || cs.size() == 0)
			return substrates;
		
		List<List<Module>> input = new ArrayList<List<Module>>();
		input.add(cs);
		input.add(de);
		List<List<Module>> output = Combinations.getCombinations(input);
		
		for (List<Module> combination : output) {
			SubstrateSet substrate = new SubstrateSet(combination);
			substrates.add(substrate);
		}

		return substrates;
	}

}
