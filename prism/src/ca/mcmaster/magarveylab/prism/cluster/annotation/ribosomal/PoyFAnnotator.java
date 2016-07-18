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
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;

/**
 * Get potential substrates for the LanM-like proteusin C-terminal dehydratase PoyF.
 * 
 * @author skinnider
 *
 */
public class PoyFAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.PoyF };
	}

	public List<SubstrateSet> findSubstrates(Domain domain,
			List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		
		if (permutation.size() == 0) 
			return substrates;
		
		if (permutation.size() == 0) {
			System.out.println("Error: could not get N-terminus of a permutation with size 0");
			return substrates;
		}

		// PoyF acts at N-terminal threonine (and, maybe, serine) residues 
		Module first = permutation.get(0);
		if (first.type() == ModuleTypes.RIBOSOMAL
				&& first.scaffold() != null
				&& first.scaffold().topSubstrate() != null
				&& first.scaffold().topSubstrate().type() == ProteinogenicAminoAcids.SERINE
				|| first.scaffold().topSubstrate().type() == ProteinogenicAminoAcids.THREONINE) {
			SubstrateSet substrate = new SubstrateSet();
			substrate.add(first);
			substrates.add(substrate);
		}

		return substrates;
	}

}
