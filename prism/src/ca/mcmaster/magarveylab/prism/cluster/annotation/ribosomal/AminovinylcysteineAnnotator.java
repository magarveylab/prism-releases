package ca.mcmaster.magarveylab.prism.cluster.annotation.ribosomal;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.clusters.RibosomalClusterTypes;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;

/**
 * Find potential substrates at which the aminovinylcysteine biosynthesis
 * flavoprotein from the linardin class of RiPPs can react, defined as a
 * C-terminal cysteine with any other cysteine (in linardins) or a C-terminal
 * cysteine with a Dha residue (in any cluster).
 * 
 * @author skinnider
 *
 */
public class AminovinylcysteineAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.LanD,
				RibosomalDomains.TvaF };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		
		if (permutation.size() == 0) {
			System.out.println("Error: could not get substrate for AviCys biosynthesis protein");
			return substrates;
		}
		
		Module last = permutation.get(permutation.size() - 1);
		if (last.type() == ModuleTypes.RIBOSOMAL
				&& last.scaffold() != null
				&& last.scaffold().topSubstrate() != null
				&& last.scaffold().topSubstrate().type() == ProteinogenicAminoAcids.CYSTEINE) {
			for (Module module : permutation) {
				if (module.type() == ModuleTypes.RIBOSOMAL && module != last) {
					Domain scaffold = module.scaffold();
					if (scaffold == null) {
						System.out.println("Error: could not get scaffold domain for ribosomal module");
						continue;					
					} else if (scaffold.topSubstrate() == null) {
						System.out.println("Error: could not get substrate for ribosomal module");
						continue;
					}
					if (scaffold.topSubstrate().type() == ProteinogenicAminoAcids.SERINE
							|| (cluster.types().contains(RibosomalClusterTypes.LINARIDIN) 
								&& scaffold.topSubstrate().type() == ProteinogenicAminoAcids.CYSTEINE)) {
						SubstrateSet substrate = new SubstrateSet();
						substrate.add(module);
						substrate.add(last);
						substrates.add(substrate);
					}
				}
			}
		}

		return substrates;
	}

}
