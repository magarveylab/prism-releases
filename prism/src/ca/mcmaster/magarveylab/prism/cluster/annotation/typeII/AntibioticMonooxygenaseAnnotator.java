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

public class AntibioticMonooxygenaseAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { TypeIIPolyketideDomains.ABM,
				TypeIIPolyketideDomains.C6OX };
	}

	public List<SubstrateSet> findSubstrates(Domain domain,
			List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		
		if (cluster.contains(TypeIIPolyketideDomains.CYCLASE_CLADE_7)) {
			// C9 first ring cyclization
			if (permutation.size() < 9)
				return substrates;
			Module first = permutation.get(3);
			Module second = permutation.get(7);
			SubstrateSet substrate = new SubstrateSet(first, second);
			substrates.add(substrate);
		} else if (cluster.contains(TypeIIPolyketideDomains.CYCLASE_CLADE_8_9)
				|| cluster.contains(TypeIIPolyketideDomains.CYCLASE_CLADE_8a)
				|| cluster.contains(TypeIIPolyketideDomains.CYCLASE_CLADE_10)) {
			// C7 first ring cyclization 
			if (permutation.size() < 7)
				return substrates;
			Module first = permutation.get(2);
			Module second = permutation.get(6);
			SubstrateSet substrate = new SubstrateSet(first, second);
			substrates.add(substrate);
		}
		
		return substrates;
	}
	
}
