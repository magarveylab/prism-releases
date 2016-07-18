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
 * Get potential substrates for thioviridamide N-terminal serine-derived monomer
 * biosynthesis. This class finds the N-terminal residue if and only if the
 * cluster also contains other domains putatively required for the biosynthesis
 * of this monomer (TvaC and TvaE).
 * 
 * @author skinnider
 *
 */
public class TvaDAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.TvaD };
	}

	public List<SubstrateSet> findSubstrates(Domain domain,
			List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		if (permutation.size() == 0) {
			System.out.println("Error: could not get N-terminus of a permutation with size 0");
			return substrates;
		}

		if (!cluster.contains(RibosomalDomains.TvaC)
				|| !cluster.contains(RibosomalDomains.TvaE)) {
			System.out.println("Error: could not get N-terminus for TvaD: "
					+ "cluster does not contain TvaC or TvaE");
			return substrates;
		}

		Module first = permutation.get(0);
		if (first.type() == ModuleTypes.RIBOSOMAL
				&& first.scaffold() != null
				&& first.scaffold().topSubstrate() != null
				&& first.scaffold().topSubstrate().type() == ProteinogenicAminoAcids.SERINE) {
			SubstrateSet substrate = new SubstrateSet();
			substrate.add(first);
			substrates.add(substrate);
		}
		
		return substrates;
	}

}
