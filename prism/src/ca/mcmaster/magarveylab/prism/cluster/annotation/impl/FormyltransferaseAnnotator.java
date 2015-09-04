package ca.mcmaster.magarveylab.prism.cluster.annotation.impl;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.enums.substrates.AdenylationSubstrates;
import ca.mcmaster.magarveylab.prism.cluster.analysis.ModuleAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;

/**
 * Get potential sites of N-formylation.<br>
 * If a cluster contains ornithine or hydroxy-ornithine, assume it can get
 * formylated.<br>
 * Also assume the first residue can be formylated.
 * 
 * @author skinnider
 *
 */
public class FormyltransferaseAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { TailoringDomains.FORMYLTRANSFERASE };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		if (permutation.size() == 0)
			return substrates;
		
		Module first = permutation.get(0);
		if (first.isAdenylationModule()) {
			SubstrateSet substrate = new SubstrateSet(first);
			substrates.add(substrate);
		}
		
		List<Module> orn = new ArrayList<Module>();
		orn.addAll(ModuleAnalyzer.modules(permutation, AdenylationSubstrates.ORNITHINE));
		orn.addAll(ModuleAnalyzer.modules(permutation, AdenylationSubstrates.N5_HYDROXYORNITHINE));
		
		for (Module module : orn)
			if (permutation.indexOf(module) > 0) {
				SubstrateSet s = new SubstrateSet(module);
				substrates.add(s);
			}

		return substrates;
	}
	
}
