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
 * Find potential sites at which YmB1 and homologs can catalyze phenyloxazoline
 * formation, and YmC1 further catalyzes phenyloxazole formation.
 * 
 * @author skinnider
 *
 */
public class YmB1Annotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.YmB1,
				RibosomalDomains.YmC1 };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		
		// get first and last modules
		Module phe = permutation.get(0);
		if (phe.type() != ModuleTypes.RIBOSOMAL
				|| phe.scaffold() == null
				|| phe.scaffold().topSubstrate() == null
				|| phe.scaffold().topSubstrate().type() != ProteinogenicAminoAcids.PHENYLALANINE) 
			return substrates;
		
		Module last = permutation.get(permutation.size() - 1);
		if (last.type() != ModuleTypes.RIBOSOMAL
				|| last.scaffold() == null
				|| last.scaffold().topSubstrate() == null)
			return substrates;
		
		SubstrateSet substrate = new SubstrateSet(phe, last);
		substrates.add(substrate);

		return substrates;
	}

}
