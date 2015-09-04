package ca.mcmaster.magarveylab.prism.cluster.analysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import ca.mcmaster.magarveylab.enums.DomainFamilies;
import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.OrfTypes;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.util.Strings;

/**
 * Analyzes orfs for combinatorial library generation.
 * @author skinnider
 *
 */
public class OrfAnalyzer {

	/**
	 * Determine whether this orf is a possible cluster end. 'End' is defined as a thioesterase domain at the end of
	 * an orf with at least one active biosynthetic module.
	 * @param orf	orf to analyze
	 * @return		true if orf is a possible cluster end
	 */
	public static boolean isEnd(Orf orf) {
		boolean flag = false;
		if (orf.activeModules().size() > 0) {
			List<Domain> domains = orf.domains(DomainFamilies.THIOTEMPLATED);
			int lastIdx = domains.size() - 1;
			Domain last = domains.get(lastIdx);
			if (last.type() == ThiotemplatedDomains.THIOESTERASE)
				flag = true;
		}
		if (flag != false) 
			System.out.println("[OrfAnalyzer] Found end orf: " + orf.name());
		return flag;
	}

	/**
	 * Determine whether this orf contains a C-starter domain within its first adenylation module.
	 * @param orf	orf to analyze
	 * @return		true if orf is a C starter
	 */
	public static boolean isCStarter(Orf orf) {
		boolean flag = false;
		List<Module> modules = orf.activeModules();
		if (modules.size() > 0) {
			Module first = modules.get(0);
			if (first.type() == ModuleTypes.C_STARTER)
				flag = true;
		}
		if (flag != false) 
			System.out.println("[OrfAnalyzer] Found C-starter orf: " + orf.name());
		return flag;
	}
		
	/**
	 * Determine whether this orf's first module contains a KS1 starter domain.
	 * @param orf	orf to analyze
	 * @return		true if orf is KS1
	 */
	public static boolean isStarterKS(Orf orf) {
		boolean flag = false;
		List<Module> modules = orf.activeModules();
		if (modules.size() > 0) {
			Module first = modules.get(0);
			if (first.type() == ModuleTypes.ACYLTRANSFERASE)
				for (Domain domain : first.domains())
					if (DomainAnalyzer.isStarterKS(domain))
						flag = true;
		}
		if (flag != false) 
			System.out.println("[OrfAnalyzer] Found KS1 orf: " + orf.name());
		return flag;
	}
	
	/**
	 * Determine whether this orf's first module is a starter acyl-adenylating module.
	 * @param orf	orf to analyze
	 * @return		true if orf is a starter acyl-adenylating module
	 */
	public static boolean isStarterAcylAdenylate(Orf orf) {
		boolean flag = false;
		List<Module> modules = orf.activeModules();
		if (modules.size() > 0) {
			Module first = modules.get(0);
			if (first.type() == ModuleTypes.ACYL_ADENYLATE && !first.canExtend())
				flag = true;
		}
		if (flag != false) 
			System.out.println("[OrfAnalyzer] Found non-extending acyl-adenylating orf: " + orf.name());
		return flag;
	}
	
	/**
	 * Determine whether this orf is a trans-adenylation orf.
	 * @param orf	orf to analyze
	 * @return		true if trans A
	 */
	public static boolean isTransAdenylation(Orf orf) {
		boolean flag = false;
		List<Module> modules = orf.activeModules();
		if (modules.size() == 1 && modules.get(0).type() == ModuleTypes.TRANS_ADENYLATION)
			flag = true;
		return flag;
	}
	
	/**
	 * Detect the type (NRPS, PKS, or hybrid) of a biosynthetic orf.
	 * Order in which types are preferentially detected: ribosomal > type II PKS > (PKS/NRPS/hybrid) > sugar > modification > 
	 * inactive > tailoring > resistance > regulator. 
	 * @param orf		the orf to search
	 * @param config	the current Prism configuration
	 */
	public static void detectType(Orf orf) {
		boolean pks = false, nrps = false;
		OrfTypes type = null;
		List<Module> modules = orf.activeModules();
		
		for (Module module : modules) {
			if (module.isAcyltransferaseModule())
				pks = true;
			if (module.isAdenylationModule())
				nrps = true;
		}
		
		if (modules.size() == 0) {
			if (orf.domains(DomainFamilies.REGULATOR).size() > 0) 
				type = OrfTypes.REGULATOR;
			if (orf.domains(DomainFamilies.RESISTANCE).size() > 0) 
				type = OrfTypes.RESISTANCE;
			if (orf.domains(DomainFamilies.PREREQUISITE).size() > 0) 
				type = OrfTypes.PREREQUISITE;
			if (orf.domains(DomainFamilies.RIBOSOMAL).size() > 0) 
				type = OrfTypes.RIBOSOMAL;
			if (orf.domains(DomainFamilies.TAILORING).size() > 0) 
				type = OrfTypes.TAILORING;
			if (orf.domains(DomainFamilies.AMINOGLYCOSIDE).size() > 0) 
				type = OrfTypes.AMINOGLYCOSIDE;
			if (orf.domains(DomainFamilies.BETA_LACTAM).size() > 0) 
				type = OrfTypes.BETA_LACTAM;
			if (orf.domains(DomainFamilies.TYPE_II_POLYKETIDE).size() > 0) 
				type = OrfTypes.TYPE_II_PKS;
			if (orf.domains(DomainFamilies.NUCLEOSIDE).size() > 0) 
				type = OrfTypes.NUCLEOSIDE;
		}
		
		if (pks == true || orf.domains(ThiotemplatedDomains.ENEDIYNE_PPTASE).size() > 0)
			type = OrfTypes.PKS;
		if (nrps == true) 
			type = OrfTypes.NRPS;
		if (pks == true && nrps == true) 
			type = OrfTypes.HYBRID;
		if (pks == false && nrps == false && type == null) {
			// if the orf isn't pks or nrps or tailoring, check whether it is inactive
			if (orf.domains(DomainFamilies.SUGAR).size() > 0) {
				type = OrfTypes.SUGAR;
			} else if (orf.hasBiosyntheticDomains()) {
				type = OrfTypes.MODIFICATION;
			} else {
				type = OrfTypes.INACTIVE;
			}
		}
		
		if (type == null)
			type = OrfTypes.NULL;
		
		System.out.println("[OrfAnalyzer] Classified " + orf.name() + " as " + type.toString());
		orf.setType(type); 
	}

	/**
	 * Get the parent open frame of a domain within a cluster.
	 * 
	 * @param domain
	 *            the domain in question
	 * @param cluster
	 *            the parent cluster
	 * @return the parent open reading frame, or null if none is identified
	 */
	public static Orf getParentOrf(Domain domain, Cluster cluster) {
		Orf parent = null;
		for (Orf orf : cluster.orfs())
			if (orf.contains(domain)) {
				parent = orf;
				break;
			}
		if (parent == null)
			System.out.println("Error: no parent orf for domain "
					+ domain.name() + " in cluster " + cluster.index());
		return parent;
	}

	/**
	 * Test whether a list of open reading frames are colinear: i.e., whether
	 * they are all in the same frame on the chromosome.
	 * 
	 * @param orfs
	 *            list of open reading frames to analyze
	 * @return true if all open reading frames are in the same frame
	 */
	public static boolean areColinear(List<Orf> orfs) {
		boolean flag = false;
		if (orfs.size() > 1) {
			List<String> frames = new ArrayList<String>();
			for (Orf orf : orfs)
				frames.add(orf.frame());
			if (!(Strings.contains("+", frames) && Strings.contains("-", frames)))
				flag = true;
		}
		return flag;
	}

	/**
	 * Determine whether a subset of orfs within a cluster are adjacent. This
	 * method works by getting the index of each orf, sorting the list in
	 * ascending order, and checking if last - first = length - 1.
	 * 
	 * @param orfs
	 *            orfs to analyze
	 * @param cluster
	 *            cluster in question
	 * @return true if the orfs are adjacent within the cluster
	 */
	public static boolean areAdjacent(List<Orf> orfs, Cluster cluster) {
		boolean flag = false;
		if (orfs.size() > 1) {
			List<Integer> indices = new ArrayList<Integer>();
			for (Orf orf : orfs)
				indices.add(cluster.orfs().indexOf(orf));
			Collections.sort(indices);
			int length = indices.size();
			int last = indices.get(length - 1);
			int first = indices.get(0);
			if (last - first == length -1)
				flag = true;
		}
		return flag;
	}
	
	/**
	 * Convert an ordered list of orfs to an ordered list of modules.
	 * 
	 * @param orfs
	 *            list of orfs
	 * @return list of modules contained within these orfs
	 */
	public static List<Module> getModules(List<Orf> orfs) {
		List<Module> modules = new ArrayList<Module>();
		for (Orf orf : orfs)
			for (Module module : orf.modules())
				modules.add(module);
		return modules;
	}

}
