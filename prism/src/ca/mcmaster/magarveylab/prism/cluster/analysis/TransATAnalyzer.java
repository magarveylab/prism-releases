package ca.mcmaster.magarveylab.prism.cluster.analysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import ca.mcmaster.magarveylab.enums.DomainFamilies;
import ca.mcmaster.magarveylab.enums.Frames;
import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.prism.combinatorialization.OrfPermuter;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.CombinatorialData;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.data.Substrate;
import ca.mcmaster.magarveylab.prism.enums.hmms.AcyltransferaseHmms;

/**
 * Analyze type I polyketide biosynthetic gene clusters with trans-acting
 * acyltransferase domains.
 * 
 * @author skinnider
 *
 */
public class TransATAnalyzer {

	/**
	 * Construct copies of all modular open reading frames in a biosynthetic
	 * gene cluster.
	 * 
	 * @param cluster
	 *            cluster in question
	 * @return copies of all orfs
	 */
	public static List<Orf> getNonTransATOrfs(Cluster cluster) {
		List<Orf> copy = new ArrayList<Orf>();
		for (Orf orf : cluster.moduleOrfs())
			copy.add(new Orf(orf));
		if (cluster.frame() == Frames.NEGATIVE)
			Collections.reverse(copy);
		return copy;
	}

	/**
	 * Construct artificial orfs to find modules that span multiple open reading
	 * frames in trans-AT clusters.
	 * 
	 * @param orfs
	 *            list of open reading frames to construct artificial orfs from
	 * @return concatenated open reading frames
	 */
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

	/**
	 * Insert trans-acting acyltransferase domains into all trans-AT insertion
	 * modules within a trans-AT polyketide gene cluster.
	 * 
	 * @param modules
	 *            module permutation to insert trans-AT domains into
	 * @param cluster
	 *            parent trans-AT cluster
	 */
	public static void insertTransAT(List<Module> modules, Cluster cluster) {
		List<Module> transAT = cluster.modules(ModuleTypes.TRANS_AT);
		Domain domain = new Domain(transAT.get(0).first());
		
		// if there are more than one trans-AT domains, simply insert malonyl-CoA
		if (transAT.size() > 1) {
			System.out.println("[TransATAnalyzer] Found more than one trans-AT module; inserting malonyl-CoA instead");
			Substrate s = domain.topSubstrate();
			if (s.type() != AcyltransferaseHmms.MALONYL_COA_1 && 
					s.type() != AcyltransferaseHmms.MALONYL_COA_2) 
				s.setType(AcyltransferaseHmms.MALONYL_COA_1);
		}
	
		for (Module module : modules)
			if (module.type() == ModuleTypes.TRANS_AT_INSERTION)
				insertTransAT(domain, module);
	}

	/**
	 * Insert a trans-acting acyltransferase domain into an insertion module
	 * (KS-[KR]-[DH]-[ER]-T with no AT).
	 * 
	 * @param transAT
	 *            the trans-acyltransferase domain to insert
	 * @param insertion
	 *            the insertion module
	 */
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
