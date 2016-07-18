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
 * Get the C-terminal cysteine residue of a ribosomal peptide. Used for
 * C-terminal oxidative decarboxylation and thiazole formation in the
 * bottromycin family (as presumably catalyzed by BotP).
 *
 * @author skinnider
 *
 */
public class BotPAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.BotP };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		if (permutation.size() == 0) {
			System.out.println("Error: could not get C-terminus of a permutation with size 0");
			return substrates;
		}

		for (Module module : permutation) {
			if (module.type() == ModuleTypes.RIBOSOMAL
					&& module.scaffold() != null
					&& module.scaffold().topSubstrate() != null
					&& module.scaffold().topSubstrate().type() == ProteinogenicAminoAcids.CYSTEINE) {
				SubstrateSet substrate = new SubstrateSet();
				substrate.add(module);
				substrates.add(substrate);
			}
		}

		return substrates;
	}

}
