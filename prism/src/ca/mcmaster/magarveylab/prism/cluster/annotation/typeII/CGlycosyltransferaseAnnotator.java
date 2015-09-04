package ca.mcmaster.magarveylab.prism.cluster.annotation.typeII;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;

/**
 * Get all type II polyketide synthase or type II starter modules for 
 * type II-specific tailoring reactions.
 * @author skinnider
 *
 */
public class CGlycosyltransferaseAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { TailoringDomains.C_GLYCOSYLTRANSFERASE };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		for (Module module : permutation) {
			if (module.type() == ModuleTypes.TYPE_II_PKS 
					|| module.type() == ModuleTypes.TYPE_II_PKS_STARTER) {
				SubstrateSet substrate = new SubstrateSet(module);
				substrates.add(substrate);
			}
		}
		return substrates;
	}

}
