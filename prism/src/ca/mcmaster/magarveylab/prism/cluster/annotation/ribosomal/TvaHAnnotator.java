package ca.mcmaster.magarveylab.prism.cluster.annotation.ribosomal;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;

/**
 * Get the residues at which TvaH, the thioviridamide thioamide-forming enzyme,
 * reacts (i.e. residues 2-6).
 * 
 * @author skinnider
 *
 */
public class TvaHAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.TvaH };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		if (permutation.size() == 0) {
			System.out.println("Error: could not get TvaH substrates in a permutation with size 0");
			return substrates;
		}
		
		SubstrateSet substrate = new SubstrateSet();
		for (int i = 1; i < 6; i++) {
			Module m = permutation.get(i);
			if (m.type() == ModuleTypes.RIBOSOMAL
					&& m.scaffold() != null
					&& m.scaffold().topSubstrate() != null)  {
				substrate.add(m);
			}
		}
		substrates.add(substrate);
		
		return substrates;
	}

}
