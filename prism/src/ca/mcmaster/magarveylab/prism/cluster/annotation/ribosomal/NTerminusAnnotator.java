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
 * Get the N-terminal residue of a ribosomal peptide. Used for linardin and
 * linear azole-containing peptide N,N-dimethylation (CypM and PznL,
 * respectively), pyruvate/2-oxobutyrate dehydrogenation (ElxO), and N-terminal
 * acetylation (MdnD, GodH, PaeN).
 * 
 * @author skinnider
 *
 */
public class NTerminusAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { 
				RibosomalDomains.CypM,
				RibosomalDomains.ElxO, 
				RibosomalDomains.PznL,
				RibosomalDomains.MdnD,
				RibosomalDomains.GodH,
				RibosomalDomains.PaeN };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		if (permutation.size() == 0) {
			System.out.println("Error: could not get N-terminus of a permutation with size 0");
			return substrates;
		}

		Module first = permutation.get(0);
		if (first.type() == ModuleTypes.RIBOSOMAL
				&& first.scaffold() != null
				&& first.scaffold().topSubstrate() != null) {
			SubstrateSet substrate = new SubstrateSet();
			substrate.add(first);
			substrates.add(substrate);
		}
		
		return substrates;
	}

}
