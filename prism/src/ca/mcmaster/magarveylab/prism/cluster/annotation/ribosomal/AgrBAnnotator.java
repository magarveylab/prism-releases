package ca.mcmaster.magarveylab.prism.cluster.annotation.ribosomal;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;

/**
 * Get the C-terminal residue and the conserved serine or cysteine residue four
 * amino acids from the C-terminus in auto-inducing peptides for
 * macrocyclization.
 * 
 * @author skinnider
 *
 */
public class AgrBAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.AgrB };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		SubstrateSet substrate = new SubstrateSet();

		if (permutation.size() < 5) {
			System.out.println("Error: could not get substrates for AgrB: "
					+ "permutation has fewer than 5 residues");
			return substrates;
		}

		Module last = permutation.get(permutation.size() - 1);
		if (last.type() == ModuleTypes.RIBOSOMAL
				&& last.scaffold() != null
				&& last.scaffold().topSubstrate() != null) {
			substrate.add(last);
		} else {
			System.out.println("Could not get last residue for AgrB");
		}

		Module negative5 = permutation.get(permutation.size() - 5);
		if (negative5.type() == ModuleTypes.RIBOSOMAL
				&& negative5.scaffold() != null
				&& negative5.scaffold().topSubstrate() != null
				&& (negative5.scaffold().topSubstrate().type() == ProteinogenicAminoAcids.SERINE 
				|| negative5.scaffold().topSubstrate().type() == ProteinogenicAminoAcids.CYSTEINE)) {
			substrate.add(negative5);
		} else {
			System.out.println("Could not get Cys/Ser residue at position -5 for AgrB");
		}
		
		if (substrate.size() == 2) 
			substrates.add(substrate);
		
		return substrates;
	}

}
