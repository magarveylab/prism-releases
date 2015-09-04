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

public class C9C14FirstRingCyclaseAnnotator implements Annotator {
	
	public DomainType[] domains() {
		return new DomainType[] { TypeIIPolyketideDomains.CYCLASE_CLADE_7 };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster) 
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>(); 
		if (permutation.size() < 9) {
			System.out.println("[C9C14FirstRingCyclaseAnnotator] Error: fewer than 9 modules detected!");
			return substrates; // need at least 6 acetate units to have C12
		}
		
		List<Module> modules = new ArrayList<Module>();
		Module c5_6 = permutation.get(2);
		Module c7_8 = permutation.get(3);
		Module c9_10 = permutation.get(4);
		Module c11_12 = permutation.get(5);
		Module c13_14 = permutation.get(6);
		Module c15_16 = permutation.get(7);
		Module c16_17 = permutation.get(8);
		modules.add(c5_6);
		modules.add(c7_8);
		modules.add(c9_10);
		modules.add(c11_12);
		modules.add(c13_14);
		modules.add(c15_16);
		modules.add(c16_17);
		
		SubstrateSet substrate = new SubstrateSet(modules);
		substrates.add(substrate);
		return substrates;
	}
	
}
