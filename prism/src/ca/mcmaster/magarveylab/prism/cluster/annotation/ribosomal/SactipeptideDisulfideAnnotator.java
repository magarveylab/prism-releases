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
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;

/**
 * Get all potential sites of disulfide bond formation in sactipeptides (where
 * only a single disulfide bond is present, at least in known compounds).
 * 
 * @author skinnider
 *
 */
public class SactipeptideDisulfideAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.SkfH };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		
		List<Module> cys = RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.CYSTEINE, permutation);
		
		ICombinatoricsVector<Module> vector = Factory.createVector(cys);
		Generator<Module> gen0 = Factory.createSimpleCombinationGenerator(vector, 2);
		for (ICombinatoricsVector<Module> subset : gen0) {
			SubstrateSet substrate = new SubstrateSet(subset.getVector());
			substrates.add(substrate);
		}
		
		return substrates;
	}

}
