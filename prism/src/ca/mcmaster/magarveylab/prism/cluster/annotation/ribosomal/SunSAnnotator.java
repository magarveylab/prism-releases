package ca.mcmaster.magarveylab.prism.cluster.annotation.ribosomal;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.analysis.RibosomalClusterAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;

/**
 * Get potential sites of S-glycosylation as catalyzed by the glycocin family
 * enzyme SunS. Note that if the tail (C-terminal) cysteine is glycosylated, it
 * is assumed that a hairpin serine may also be glycosylated.
 * 
 * @author skinnider
 *
 */
public class SunSAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.SunS };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		List<Module> cysteines = RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.CYSTEINE, permutation);
		List<Module> serines = RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.SERINE, permutation);
		
		for (Module cysteine : cysteines) {
			SubstrateSet substrate = new SubstrateSet(cysteine);
			substrates.add(substrate);

			// if the tail cysteine is glycosylated, so too will be a hairpin
			// serine
			int cysIdx = cysteines.indexOf(cysteine);
			if (cysIdx == cysteines.size() - 1) {
				if (cysIdx == 0)
					continue;
				Module lastCys = cysteines.get(cysIdx - 1);
				
				for (Module serine : serines) {
					// needs to be before the last cysteine
					int serIdx = permutation.indexOf(serine);
					if (serIdx > permutation.indexOf(lastCys))
						continue;
					
					SubstrateSet substrate2 = new SubstrateSet(cysteine, serine);
					substrates.add(substrate2);
				}
			}
		}
		
		return substrates;
	}

}
