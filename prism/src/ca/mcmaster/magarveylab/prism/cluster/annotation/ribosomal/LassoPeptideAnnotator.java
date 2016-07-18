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
 * Find possible sites at which a lasso peptide asparagine synthase can cyclize
 * a lasso peptide precursor. Typically, this occurs between the N-terminal
 * residue (usually a glycine) and an aspartate or glutamate at position 8 or 9,
 * but may also occur at positions 6 or 10.
 * 
 * @author skinnider
 *
 */
public class LassoPeptideAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.Asparagine_synthase };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		
		if (permutation.size() < 9) {
			System.out.println("Error: could not get substrates for lasso peptide "
					+ "asparagine synthase in module permutation with < 9 amino acids!");
			return substrates;
		}

		Module first = permutation.get(0);
		
		for (int i = 0; i < permutation.size(); i++) {
			Module module = permutation.get(i);
			if (i == 8 || i == 7 || i == 6 || i == 9) {
				if (module.type() == ModuleTypes.RIBOSOMAL) {
					Domain scaffold = module.scaffold();
					if (scaffold == null) {
						System.out.println("Error: could not get scaffold domain for ribosomal module");
						return substrates;					
					} else if (scaffold.topSubstrate() == null) {
						System.out.println("Error: could not get substrate for ribosomal module");
						return substrates;										
					}
					if (scaffold.topSubstrate().type() == ProteinogenicAminoAcids.ASPARTATE 
							|| scaffold.topSubstrate().type() == ProteinogenicAminoAcids.GLUTAMATE) {
						SubstrateSet substrate = new SubstrateSet();
						substrate.add(first);
						substrate.add(module);
						substrates.add(substrate);
					}
				}
			}
		}
		
		return substrates;
	}

}
