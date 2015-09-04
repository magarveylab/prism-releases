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
 * Find potential sites at which a type II polyketide cyclase could catalyze the
 * linear cyclization of the A and B rings, with the A ring cyclizing on C7-C12.
 * 
 * @author skinnider
 *
 */
public class C7ABRingCyclaseAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { TypeIIPolyketideDomains.CYCLASE_CLADE_10,
				TypeIIPolyketideDomains.CYCLASE_CLADE_8a };
	}

	public List<SubstrateSet> findSubstrates(Domain domain,
			List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		if (permutation.size() < 7) {
			System.out.println("[C7ABRingCyclaseAnnotator] Error: fewer than 7 modules detected!");
			return substrates;
		}

		List<Module> modules = new ArrayList<Module>();
		Module c5_6 = permutation.get(2);
		Module c7_8 = permutation.get(3);
		Module c9_10 = permutation.get(4);
		Module c11_12 = permutation.get(5);
		Module c13_14 = permutation.get(6);
		modules.add(c5_6);
		modules.add(c7_8);
		modules.add(c9_10);
		modules.add(c11_12);
		modules.add(c13_14);

		SubstrateSet substrate = new SubstrateSet(modules);
		substrates.add(substrate);
		return substrates;
	}

}
