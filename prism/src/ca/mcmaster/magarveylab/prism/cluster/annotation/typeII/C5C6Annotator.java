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
 * Get the C5 and C6 carbons of a type II polyketide scaffold.
 * 
 * @author skinnider
 *
 */
public class C5C6Annotator implements Annotator {
	
	public DomainType[] domains() {
		return new DomainType[] { TypeIIPolyketideDomains.C6CMT };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster) 
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>(); 
		if (permutation.size() < 3)
			return substrates;
		Module c5 = permutation.get(2);
		SubstrateSet substrate = new SubstrateSet(c5);
		substrates.add(substrate);
		return substrates;
	}
	
}
