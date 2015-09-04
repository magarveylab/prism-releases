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
 * Get the C1 and C2 carbons of a type II polyketide scaffold.
 * 
 * @author skinnider
 *
 */
public class C1C2Annotator implements Annotator {
	
	public DomainType[] domains() {
		return new DomainType[] { TypeIIPolyketideDomains.C2AMT, 
				TypeIIPolyketideDomains.C2MT,
				TypeIIPolyketideDomains.CARBOXY_OMT };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster) 
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>(); 
		if (permutation.size() < 1)
			return substrates;
		Module c1 = permutation.get(0);
		SubstrateSet substrate = new SubstrateSet(c1);
		substrates.add(substrate);
		return substrates;
	}
	
}
