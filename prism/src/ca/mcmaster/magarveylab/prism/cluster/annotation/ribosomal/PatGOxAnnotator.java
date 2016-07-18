package ca.mcmaster.magarveylab.prism.cluster.annotation.ribosomal;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;

/**
 * Find reaction sites for the bifunctional PatG protease with an oxidation
 * domain built in, which catalyzes both macrocyclization and thiazole/oxazole
 * formation (oxidation of thiazolines/oxazolines).
 * 
 * @author skinnider
 *
 */
public class PatGOxAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.PatG_ox };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		List<SubstrateSet> azolines = new AzoleAnnotator().findSubstrates(
				domain, permutation, cluster);
		List<SubstrateSet> macrocyclization = new MacrocyclizationAnnotator()
				.findSubstrates(domain, permutation, cluster);

		if (macrocyclization.size() == 0)
			return substrates;
		
		SubstrateSet macrocycle = macrocyclization.get(0);
		if (azolines.size() > 0) {
			for (SubstrateSet azoline : azolines) {
				SubstrateSet substrate = new SubstrateSet();
				substrate.addAll(macrocycle.getAllModules());
				substrate.addAll(azoline.getAllModules());
				substrates.add(substrate);
			}
		} else {
			SubstrateSet substrate = new SubstrateSet();
			substrate.addAll(macrocycle.getAllModules());
			substrates.add(substrate);
		}
		
		return substrates;
	}

}
