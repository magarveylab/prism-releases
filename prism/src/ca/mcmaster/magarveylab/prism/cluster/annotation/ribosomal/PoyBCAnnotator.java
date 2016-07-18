package ca.mcmaster.magarveylab.prism.cluster.annotation.ribosomal;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.analysis.RibosomalClusterAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.cluster.reactions.RibosomalUtil;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;

/**
 * Get potential substrates for the proteusin iterative radical SAM superfamily
 * methyltransferases PoyB/PoyC.
 * 
 * @author skinnider
 *
 */
public class PoyBCAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.PoyBC };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		// PoyBC is known to act at I,V,T,Q,M beta carbons  
		List<Module> modules = new ArrayList<Module>();
		modules.addAll(RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.ISOLEUCINE, permutation));
		modules.addAll(RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.VALINE, permutation));
		modules.addAll(RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.THREONINE, permutation));
		modules.addAll(RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.GLUTAMINE, permutation));
		modules.addAll(RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.METHIONINE, permutation));
		
		// get subsets of possibel modules 
		List<List<Module>> subsets = RibosomalUtil.getSubsetsWithSizeRange(modules, 10, 15);

		// convert to modules
		for (List<Module> subset : subsets) {
			SubstrateSet substrate = new SubstrateSet(subset);
			substrates.add(substrate);
		}

		return substrates;
	}

}
