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
 * Get the N-terminal residue and the 4th residue for macrocyclodehydration in
 * the bottromycin family of ribosomal peptides.
 * 
 * @author skinnider
 *
 */
public class BotCAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.BotC };
	}

	public List<SubstrateSet> findSubstrates(Domain domain,
			List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		SubstrateSet substrate = new SubstrateSet();

		if (permutation.size() < 5) {
			System.out.println("Error: could not get substrates for BotC: "
					+ "permutation has fewer than 4 residues");
			return substrates;
		}

		Module first = permutation.get(0);
		if (first.type() == ModuleTypes.RIBOSOMAL
				&& first.scaffold() != null
				&& first.scaffold().topSubstrate() != null) {
			substrate.add(first);
		} else {
			System.out.println("Could not get N-terminal residue for BotC");
		}

		Module fourth = permutation.get(3);
		if (fourth.type() == ModuleTypes.RIBOSOMAL
				&& fourth.scaffold() != null
				&& fourth.scaffold().topSubstrate() != null) {
			substrate.add(fourth);
		} else {
			System.out.println("Could not get fourth residue for BotC");
		}
		
		Module fifth = permutation.get(4);
		if (fifth.type() == ModuleTypes.RIBOSOMAL
				&& fifth.scaffold() != null
				&& fifth.scaffold().topSubstrate() != null) {
			substrate.add(fifth);
		} else {
			System.out.println("Could not get fifth residue for BotC");
		}
		
		if (substrate.size() == 3) 
			substrates.add(substrate);
		
		return substrates;
	}

}
