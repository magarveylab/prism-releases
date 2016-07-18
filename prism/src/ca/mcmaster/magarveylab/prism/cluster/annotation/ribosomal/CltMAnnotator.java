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
 * Get all combinations of cysteine substrates with dehydroalanine residues for
 * tertiary thioether formation, as catalyzed by the cyclothiazolomycin enzmye
 * CltM.
 * 
 * @author skinnider
 *
 */
public class CltMAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.CltM };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		
		List<Module> cys = RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.CYSTEINE, permutation);
		List<Module> ser = RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.SERINE, permutation);
		
		List<List<Module>> input = new ArrayList<List<Module>>();
		input.add(cys);
		input.add(ser);
		List<List<Module>> output = Combinations.getCombinations(input);
		
		for (List<Module> combination : output) {
			int cysIdx = permutation.indexOf(combination.get(0));
			int serIdx = permutation.indexOf(combination.get(1));
			
			// can't be N-terminal serine -- this will be pyruvate 
			if (serIdx == 0)
				continue;
			// reject strained rings  
			if (Math.abs(serIdx - cysIdx) <= 10)
				continue; 
			
			SubstrateSet substrate = new SubstrateSet(combination);
			substrates.add(substrate);
		}

		return substrates;
	}

}
