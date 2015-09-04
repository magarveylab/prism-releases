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

public class C3CRingCyclaseAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { TypeIIPolyketideDomains.CYCLASE_CLADE_2, TypeIIPolyketideDomains.CYCLASE_CLADE_3 };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster) 
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>(); 
		if (permutation.size() < 8) {
			System.out.println("[C3C16ThirdRingCyclaseAnnotator] Error: fewer than 8 modules detected!");
			return substrates;
		}
		
		List<Module> modules = new ArrayList<Module>();
		Module c3_4 = permutation.get(1);
		Module c15_c16 = permutation.get(7);
		modules.add(c3_4);
		modules.add(c15_c16);
		
		SubstrateSet substrate = new SubstrateSet(modules);
		substrates.add(substrate);
		return substrates;
	}


}
