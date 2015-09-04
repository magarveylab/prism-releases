package ca.mcmaster.magarveylab.prism.cluster.annotation.typeII;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.TypeIIPolyketideDomains;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;

/**
 * Find potential sites at which an aureolic acid Baeyer-Villiger monooxygenase
 * could catalyze the ring opening reaction that gives the characteristic
 * aureolic acid scaffold.
 * 
 * @author skinnider
 *
 */
public class BVMOAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { TypeIIPolyketideDomains.BVMO };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		if (permutation.size() < 10)
			return substrates;
		
		List<Module> modules = new ArrayList<Module>();
		Module c1_2 = permutation.get(0);
		Module c3_4 = permutation.get(1);
		Module c15_16 = permutation.get(7);
		Module c17_18 = permutation.get(8);
		modules.add(c1_2);
		modules.add(c3_4);
		modules.add(c15_16);
		modules.add(c17_18);

		SubstrateSet substrate = new SubstrateSet(modules);
		substrates.add(substrate);
		return substrates;
	}

}
