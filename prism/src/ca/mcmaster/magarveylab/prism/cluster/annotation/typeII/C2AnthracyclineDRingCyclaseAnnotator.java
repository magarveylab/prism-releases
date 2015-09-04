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

public class C2AnthracyclineDRingCyclaseAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { TypeIIPolyketideDomains.CYCLASE_CLADE_6a, 
				TypeIIPolyketideDomains.CYCLASE_CLADE_6b_SUBTYPE_1 };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster) 
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>(); 
		if (permutation.size() < 10) {
			System.out.println("[AnthracyclineFourthRingCyclaseAnnotator] Error: fewer than 10 modules detected!");
			return substrates;
		}
		
		List<Module> modules = new ArrayList<Module>();
		Module c1_2 = permutation.get(0);
		Module c17_18 = permutation.get(8);
		Module c19_20 = permutation.get(9);
		modules.add(c1_2);
		modules.add(c17_18);
		modules.add(c19_20);
		
		SubstrateSet substrate = new SubstrateSet(modules);
		substrates.add(substrate);
		return substrates;
	}

}
