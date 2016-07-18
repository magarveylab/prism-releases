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
 * Get the C-terminal residue of a ribosomal peptide. Used for C-terminal amide
 * formation (NosA/NocA/BerI, TsrC/SioC via a different mechanism) and amino
 * acetone formation (TpaJ).
 *
 * @author skinnider
 *
 */
public class CTerminusAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { 
				RibosomalDomains.NosA, 
				RibosomalDomains.TsrC,
				RibosomalDomains.TsrB,
				RibosomalDomains.TpaJ };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		if (permutation.size() == 0) {
			System.out.println("Error: could not get C-terminus of a permutation with size 0");
			return substrates;
		}

		Module last = permutation.get(permutation.size() - 1);
		if (last.type() == ModuleTypes.RIBOSOMAL
				&& last.scaffold() != null
				&& last.scaffold().topSubstrate() != null) {
			SubstrateSet substrate = new SubstrateSet();
			substrate.add(last);
			substrates.add(substrate);
		}
		
		return substrates;
	}

}
