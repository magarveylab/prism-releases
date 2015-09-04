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

public class C5CRingCyclaseAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { TypeIIPolyketideDomains.CYCLASE_CLADE_1 };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster) 
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>(); 
		if (permutation.size() < 9) {
			System.out.println("[C5C18ThirdRingCyclaseAnnotator] Error: fewer than 9 modules detected!");
			return substrates;
		}
		
		List<Module> modules = new ArrayList<Module>();
		Module c5_6 = permutation.get(2);
		Module c17_c18 = permutation.get(8);
		modules.add(c5_6);
		modules.add(c17_c18);
		
		SubstrateSet substrate = new SubstrateSet(modules);
		substrates.add(substrate);
		return substrates;
	}

}
