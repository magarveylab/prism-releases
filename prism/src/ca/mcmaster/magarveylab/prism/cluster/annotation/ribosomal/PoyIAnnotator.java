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
 * Get potential substrates for the proteusin iterative oxidoreductase PoyI,
 * which catalyzes beta-hydroxylation of asparagine and valine residues in
 * polytheonamide.
 * 
 * @author skinnider
 *
 */
public class PoyIAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.PoyI };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		// PoyI is known to act at V,N beta carbons  
		List<Module> modules = new ArrayList<Module>();
		modules.addAll(RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.ASPARAGINE, permutation));
		modules.addAll(RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.VALINE, permutation));
		
		// get subsets of possible modules 
		List<List<Module>> subsets = RibosomalUtil.getSubsetsWithSizeRange(modules, 3, 5);

		// convert to modules
		for (List<Module> subset : subsets) {
			SubstrateSet substrate = new SubstrateSet(subset);
			substrates.add(substrate);
		}

		return substrates;
	}

}
