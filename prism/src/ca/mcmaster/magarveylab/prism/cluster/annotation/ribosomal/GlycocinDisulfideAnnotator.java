package ca.mcmaster.magarveylab.prism.cluster.annotation.ribosomal;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;
import org.paukov.combinatorics.Factory;
import org.paukov.combinatorics.Generator;
import org.paukov.combinatorics.ICombinatoricsVector;
import org.paukov.combinatorics.util.ComplexCombinationGenerator;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.prism.cluster.analysis.RibosomalClusterAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;

/**
 * Get all potential sites of disulfide bond formation in glycocins (where only
 * 2 disulfide bonds are considered).<br>
 * <br>
 * BdbB does not catalyze disulfide bond formation in the glycocin F cluster, so
 * this has been replaced with a SunA annotator.
 * 
 * @author skinnider
 *
 */
@Deprecated	
public class GlycocinDisulfideAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		List<Module> cys = RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.CYSTEINE, permutation);
		if (cys.size() != 5) {
			List<List<Module>> subsets = new ArrayList<List<Module>>();
			ICombinatoricsVector<Module> vector = Factory.createVector(cys);
			Generator<Module> gen1 = Factory.createSimpleCombinationGenerator(vector, 4);
			for (ICombinatoricsVector<Module> subset : gen1) 
				subsets.add(subset.getVector());
			
			for (List<Module> subset : subsets) {
				SubstrateSet substrate = new SubstrateSet(subset);
				substrates.add(substrate);
			}
		} else {
			List<Module> c1234 = new ArrayList<Module>();
			List<Module> c1245 = new ArrayList<Module>();
			c1234.add(cys.get(0));
			c1234.add(cys.get(1));
			c1234.add(cys.get(2));
			c1234.add(cys.get(3));
			c1245.add(cys.get(0));
			c1245.add(cys.get(1));
			c1245.add(cys.get(3));
			c1245.add(cys.get(4));
			
			ICombinatoricsVector<Module> v1 = Factory.createVector(c1234);
			ICombinatoricsVector<Module> v2 = Factory.createVector(c1245);
			Generator<ICombinatoricsVector<Module>> gen1 = new ComplexCombinationGenerator<Module>(v1, 2);
			Generator<ICombinatoricsVector<Module>> gen2 = new ComplexCombinationGenerator<Module>(v2, 2);
			for (ICombinatoricsVector<ICombinatoricsVector<Module>> subset : gen1)
				if (subset.getVector().get(0).getVector().size() == 2) {
					SubstrateSet substrate = new SubstrateSet();
					substrate.addAll(subset.getVector().get(0).getVector());
					substrate.addAll(subset.getVector().get(1).getVector());
				}
			for (ICombinatoricsVector<ICombinatoricsVector<Module>> subset : gen2)
				if (subset.getVector().get(0).getVector().size() == 2) {
					SubstrateSet substrate = new SubstrateSet();
					substrate.addAll(subset.getVector().get(0).getVector());
					substrate.addAll(subset.getVector().get(1).getVector());
				}
		}
		
		return substrates;
	}

}
