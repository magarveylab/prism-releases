package ca.mcmaster.magarveylab.prism.cluster.annotation.impl;

import java.io.IOException;
import java.util.List;

import org.openscience.cdk.exception.CDKException;

import ca.mcmaster.magarveylab.enums.ClusterFamilies;
import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.domains.BetaLactamDomains;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;
import ca.mcmaster.magarveylab.prism.cluster.analysis.ModuleAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.analysis.RibosomalClusterAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.AnnotatorUtil;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.prism.enums.hmms.AdenylationHmms;

/**
 * Find sites at which tailoring reactions which occur at hydroxyls (e.g.
 * carbamoylation, glycosylation, sulfonation) could occur.
 * 
 * @author skinnider
 *
 */
public class HydroxylAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { TailoringDomains.SULFOTRANSFERASE,
				TailoringDomains.CARBAMOYLTRANSFERASE,
				TailoringDomains.GLYCOSYLTRANSFERASE };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster) 
			throws IOException, CDKException {
		List<Module> modules = ModuleAnalyzer.hydroxyls(permutation);
		
		if (cluster.contains(BetaLactamDomains.DOACS)) {
			modules.addAll(ModuleAnalyzer.modules(permutation, AdenylationHmms.VALINE_1));
			modules.addAll(ModuleAnalyzer.modules(permutation, AdenylationHmms.VALINE_2));
			modules.addAll(ModuleAnalyzer.modules(permutation, AdenylationHmms.VALINE_3));
		}
		
		if (cluster.families().contains(ClusterFamilies.TYPE_II_POLYKETIDE))
			for (Module module : permutation)
				if (module.type() == ModuleTypes.TYPE_II_PKS 
						|| module.type() == ModuleTypes.TYPE_II_PKS_STARTER)
					modules.add(module);
		if (cluster.families().contains(ClusterFamilies.RIBOSOMAL)) {
			modules.addAll(RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.THREONINE, permutation));
			modules.addAll(RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.SERINE, permutation));
		}

		List<SubstrateSet> substrates = AnnotatorUtil.convertModulesToSubstrateSets(modules);
		return substrates;
	}
	
}
