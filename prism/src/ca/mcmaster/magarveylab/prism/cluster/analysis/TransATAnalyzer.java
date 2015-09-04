package ca.mcmaster.magarveylab.prism.cluster.analysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import ca.mcmaster.magarveylab.enums.DomainFamilies;
import ca.mcmaster.magarveylab.enums.Frames;
import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.enums.substrates.AcyltransferaseSubstrates;
import ca.mcmaster.magarveylab.prism.combinatorialization.OrfPermuter;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.CombinatorialData;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.data.Substrate;

public class TransATAnalyzer {

	public static List<Orf> getNonTransATOrfs(Cluster cluster) {
		List<Orf> copy = new ArrayList<Orf>();
		for (Orf orf : cluster.moduleOrfs())
			copy.add(new Orf(orf));
		if (cluster.frame() == Frames.NEGATIVE)
			Collections.reverse(copy);
		return copy;
	}

	public static List<Orf> constructArtificialOrfs(List<Orf> orfs) {
		ArtificialOrfCreator aoc = new ArtificialOrfCreator(orfs);
		aoc.initiate();
	
		for (Orf orf : orfs) {
			StringBuffer sb = new StringBuffer();
			sb.append("[TransATAnalyzer] Created artificial orf with domains: ");
			for (Domain domain : orf.domains(DomainFamilies.THIOTEMPLATED))
				sb.append(domain.type().abbreviation() + " - ");
			System.out.println(sb.toString());
		}

		return orfs;
	}

	public static void insertTransAT(List<Module> modules, Cluster cluster) {
		List<Module> transAT = cluster.modules(ModuleTypes.TRANS_AT);
		Domain domain = new Domain(transAT.get(0).first());
		
		// if there are more than one trans-AT domains, simply insert malonyl-CoA
		if (transAT.size() > 1) {
			System.out.println("[TransATAnalyzer] Found more than one trans-AT module; inserting malonyl-CoA instead");
			Substrate s = domain.topSubstrate();
			if (s.type() != AcyltransferaseSubstrates.MALONYL_COA_1 && 
					s.type() != AcyltransferaseSubstrates.MALONYL_COA_2) 
				s.setType(AcyltransferaseSubstrates.MALONYL_COA_1);
		}
	
		for (Module module : modules)
			if (module.type() == ModuleTypes.TRANS_AT_INSERTION)
				insertTransAT(domain, module);
	}

	public static void insertTransAT(Domain transAT, Module insertion) {
		if (insertion.contains(ThiotemplatedDomains.KETOSYNTHASE)) {
			// if KS, insert after KS
			int index = insertion.firstIndexOf(ThiotemplatedDomains.KETOSYNTHASE);
			insertion.domains().add(index + 1, transAT);
		} else {
			// else, insert at beginning
			insertion.domains().add(0, transAT);
		}
		System.out.println("[TransATAnalyzer] Inserted trans-AT module with substrate "
						+ transAT.topSubstrate().type().toString());
		insertion.setType(ModuleTypes.ACYLTRANSFERASE);
	}

	public static List<List<Module>> getTransATModulePermutations(Cluster cluster) {
		List<List<Module>> permutations = new ArrayList<List<Module>>();
		CombinatorialData cd = cluster.combinatorialData();

		if (cluster.isColinear()) {
			List<Orf> copy = getNonTransATOrfs(cluster);
			List<Orf> artificial = constructArtificialOrfs(copy);
			cd.setNumOrfPermutations(artificial.size());
			List<Module> modules = OrfAnalyzer.getModules(artificial);
			insertTransAT(modules, cluster);
			permutations.add(modules);
		} else {
			List<Orf> copy = getNonTransATOrfs(cluster);
			Orf start = ClusterAnalyzer.start(copy);
			Orf end = ClusterAnalyzer.end(copy);
			List<List<Orf>> orfPermutations = OrfPermuter.permuteOrfs(start, end, copy, cluster);
	
			int artificialSize = 0;
			for (List<Orf> orfPermutation : orfPermutations) {
				List<Orf> artificial = constructArtificialOrfs(orfPermutation);
				artificialSize += artificial.size();
				List<Module> modules = OrfAnalyzer.getModules(artificial);
				insertTransAT(modules, cluster);
				permutations.add(modules);
			}
			cd.setNumOrfPermutations(artificialSize);
		}
		
		return permutations;
	}

}
