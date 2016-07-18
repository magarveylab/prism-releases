package ca.mcmaster.magarveylab.prism.cluster.annotation.ribosomal;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;
import org.paukov.combinatorics.Factory;
import org.paukov.combinatorics.Generator;
import org.paukov.combinatorics.ICombinatoricsVector;

import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.enums.interfaces.SubstrateType;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;

/**
 * Get potential sites for prenylation as catalyzed by the cyanobactin
 * prenyltransferase PatF: threonines, serines, tyrosines, arginines, and the
 * N-terminal residue.
 * 
 * @author skinnider
 *
 */
public class PatFAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.PatF };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		if (permutation.size() == 0) {
			System.out.println("Error: could not get N-terminus of a permutation with size 0");
			return substrates;
		}

		List<Module> modules = new ArrayList<Module>();
		for (int i = 0; i < permutation.size(); i++) {
			Module module = permutation.get(i);
			if (module.type() == ModuleTypes.RIBOSOMAL
					&& module.scaffold() != null
					&& module.scaffold().topSubstrate() != null) {
				SubstrateType aa = module.scaffold().topSubstrate().type();
				if (i == 0) {
					modules.add(module);
				} else if (aa == ProteinogenicAminoAcids.THREONINE
						|| aa == ProteinogenicAminoAcids.SERINE 
						|| aa == ProteinogenicAminoAcids.TYROSINE
						|| aa == ProteinogenicAminoAcids.ARGININE) {
					modules.add(module);
				}
			}			
		}

		// now get all combinations of size 0, 1, and 2 
		List<List<Module>> subsets = new ArrayList<List<Module>>();
		ICombinatoricsVector<Module> vector = Factory.createVector(modules);
//		Generator<Module> gen0 = Factory.createSimpleCombinationGenerator(vector, 0);
		Generator<Module> gen1 = Factory.createSimpleCombinationGenerator(vector, 1);
		Generator<Module> gen2 = Factory.createSimpleCombinationGenerator(vector, 2);
//		for (ICombinatoricsVector<Module> subset : gen0) 
//			subsets.add(subset.getVector());
		for (ICombinatoricsVector<Module> subset : gen1) 
			subsets.add(subset.getVector());
		for (ICombinatoricsVector<Module> subset : gen2) 
			subsets.add(subset.getVector());

		for (List<Module> subset : subsets) {
			SubstrateSet substrate = new SubstrateSet(subset);
			substrates.add(substrate);
		}
		
		return substrates;
	}

}
