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
 * Get the C- and N-terminal modules for macrocyclization, as catalyzed by the
 * cyanobactin protease PatG, the bacterial head-to-tail cyclized peptide
 * DUF95-family protein, the SkfC-type and split AlbE-type sactipeptide
 * proteases, and the YM-216391 family cyclase YmF. This annotator is also used
 * to find substrates for the bifunctional
 * O-methyltransferase/N-prenyltransferase AgeF1.
 * 
 * @author skinnider
 *
 */
public class MacrocyclizationAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.PatG,
				RibosomalDomains.AgeF1, RibosomalDomains.DUF95,
				RibosomalDomains.YmF, RibosomalDomains.SkfC,
				RibosomalDomains.AlbE };
	}

	public List<SubstrateSet> findSubstrates(Domain domain,
			List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		SubstrateSet substrate = new SubstrateSet();

		if (permutation.size() == 0) 
			return substrates;
		
		// AlbE requires the AlbF domain as well 
		if (domain.type() == RibosomalDomains.AlbE
				&& !cluster.contains(RibosomalDomains.AlbF))
			return substrates;
		
		Module first = permutation.get(0);
		if (first.type() == ModuleTypes.RIBOSOMAL
				&& first.scaffold() != null
				&& first.scaffold().topSubstrate() != null) {
			substrate.add(first);
		} else {
			System.out.println("Error: could not get N-terminus for module permutation");
		}
		
		Module last = permutation.get(permutation.size() - 1);
		if (last.type() == ModuleTypes.RIBOSOMAL
				&& last.scaffold() != null
				&& last.scaffold().topSubstrate() != null
				&& last != first) {
			substrate.add(last);
		} else {
			System.out.println("Error: could not get C-terminus for module permutation");
		}
		
		if (substrate.size() == 2) 
			substrates.add(substrate);
		
		return substrates;
	}

}
