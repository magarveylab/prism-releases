package ca.mcmaster.magarveylab.prism.cluster.annotation.ribosomal;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.ModuleTypes;
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
 * Get potential substrates for modified quinaldic acid attachment in the
 * thiostrepton and siomycin clusters. The ketone is assumed to be esterified or
 * thioesterified at free threonines and serines (although only the former is
 * known to occur), while the N-terminal amide is assumed to attach to a
 * quinaldic acid carbon.
 * 
 * @author skinnider
 *
 */
public class TsrIAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.TsrI };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		
		if (permutation.size() == 0) {
			System.out.println("Error: could not get N-terminus of a permutation with size 0");
			return substrates;
		}

		Module first = permutation.get(0);
		if (first.type() != ModuleTypes.RIBOSOMAL
				|| first.scaffold() == null
				|| first.scaffold().topSubstrate() == null) 
			return substrates;

		List<Module> serthr = new ArrayList<Module>();
		serthr.addAll(RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.THREONINE, permutation));
		serthr.addAll(RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.SERINE, permutation));
				
		if (serthr.size() == 0)
			return substrates;

		for (Module module : serthr) {
			if (permutation.indexOf(module) > 0) {
				SubstrateSet substrate = new SubstrateSet(first, module);
				substrates.add(substrate);
			}
		}
		
		return substrates;
	}

}
