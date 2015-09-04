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

public class C2DERingCyclaseAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { TypeIIPolyketideDomains.CYCLASE_CLADE_5a };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster) 
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>(); 
		if (permutation.size() < 12) {
			System.out.println("[C2C23AngularCyclaseAnnotator] Error: fewer than 12 modules detected!");
			return substrates;
		}
		
		List<Module> modules = new ArrayList<Module>();
		Module c1_2 = permutation.get(0);
		Module c3_4 = permutation.get(1);
		Module c19_20 = permutation.get(9);
		Module c21_22 = permutation.get(10);
		Module c23_24 = permutation.get(11);
		modules.add(c1_2);
		modules.add(c3_4);
		modules.add(c19_20);
		modules.add(c21_22);
		modules.add(c23_24);
		
		SubstrateSet substrate = new SubstrateSet(modules);
		substrates.add(substrate);
		return substrates;
	}
	
}
