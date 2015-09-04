package ca.mcmaster.magarveylab.prism.cluster.analysis;

import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.enums.CStarterSubstrates;
import ca.mcmaster.magarveylab.enums.DomainFamilies;
import ca.mcmaster.magarveylab.enums.SubstrateDomainSearches;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.prism.blast.BlastSearchResult;
import ca.mcmaster.magarveylab.prism.data.Domain;

/**
 * Analyzes domains for combinatorial library generation.
 * @author skinnider
 *
 */
public class DomainAnalyzer {
	
	/**
	 * Get all domains of a given type from a list of domains.
	 * @param type	the type of domain to return
	 * @return	all domains of this type
	 */
	public static List<Domain> domains(DomainType type, List<Domain> domains) {
		List<Domain> typeDomains = new ArrayList<Domain>();
		for (Domain domain : domains) 
			if (domain.type() == type)
				typeDomains.add(domain);
		return typeDomains;
	}	

	/**
	 * Get all unique nonbiosynthetic domains within a list of untyped domains.
	 * @param domains	list of domains to analyze
	 * @return			all unique domains within this list
	 */
	public static List<Domain> getUniqueNonbiosyntheticDomains(List<Domain> domains) {
		ArrayList<Domain> nonbiosynthetic = new ArrayList<Domain>();
		for (Domain domain : domains) 
			if (!nonbiosynthetic.contains(domain)) 
				nonbiosynthetic.add(domain);
		return nonbiosynthetic;
	}
	
	/**
	 * Determine whether this domain is a C-starter domain.
	 * @param domain	domain to analyze
	 * @return			true if it is a C* domain
	 */
	public static boolean isCStarter(Domain domain) {
		boolean flag = false;
		if (domain.type() == ThiotemplatedDomains.CONDENSATION && domain.blastResults().size() > 0)
			for (String substrate : CStarterSubstrates.names()) {
				String name = domain.blastResults().get(0).subject();
				if (name.contains(substrate) || name.endsWith("_start"))
					flag = true;
			}
		return flag;
	}

	/**
	 * Determine whether this domain is an epimerization domain.
	 * @param domain	domain to analyze
	 * @return			true if epimerization
	 */
	public static boolean isEpimerase(Domain domain) {
		boolean flag = false;
		if (domain.type() == ThiotemplatedDomains.CONDENSATION)
			if (domain.blastResults().size() > 0 && (domain.blastResults().get(0).subject().endsWith("E") ||
					domain.blastResults().get(0).subject().endsWith("_dual")))
					flag = true;
		return flag;
	}

	/**
	 * Determine whether this domain is a "KSQ" or other module-starting ketosynthase domain
	 * @param domain	domain to analyze
	 * @return			true if KS starter
	 */
	public static boolean isStarterKS(Domain domain) {
		boolean flag = false;
		String[] starterKS = { "SorA_Q9ADL6_1KSB", "PikAI_Q9ZGI5_1KSB", "SP01_3145_1", "StiA_Q8RJY6_1KSB",
				"Stro1024_KS1", "TetA_BAE93722_KS1", "TylGI_O33954_KS1", "ChlA1_AAZ77693_1KSB", "FurA1_ABB88519_KSB", 
				"JamE_AAS98777_KS1", "MerB_ABJ97437_1KSB", "MtaB_Q9RFL0_1KSB", "MxaF_Q93TW6_1KSB", "NidA1_O30764_1KSB", 
				"OleAI_Q9KIV4_1KSB" };
		if (domain.type() == ThiotemplatedDomains.KETOSYNTHASE) {
			if (domain.blastResults().size() > 0) {
				String name = domain.blastResults().get(0).subject();
				for (String starter : starterKS)
					if (name.contains(starter)) 
						flag = true;
			}
		}
		if (flag == true)
			System.out.println("[DomainsUtil] Found starter KS domain matching " + domain.blastResults().get(0).subject());
		return flag;
	}

	/**
	 * Determine whether this domain is a cyclization domain.
	 * @param domain	domain to analyze
	 * @return			true if cyclization
	 */
	public static boolean isCyclization(Domain domain) {
		boolean flag = false;
		if (domain.type() == ThiotemplatedDomains.CONDENSATION)
			if (domain.blastResults().size() > 0 && domain.blastResults().get(0).subject().endsWith("_cyc"))
				flag = true;
		return flag;
	}
	
	/**
	 * Determine whether this domain type is associated with a domain-specific substrate search. 
	 * @param type	domain type to analyze
	 * @return		true if this domain type has a substrate 
	 */
	public static boolean isSubstrateDomainType(DomainType type) {
		boolean flag = false;
		for (SubstrateDomainSearches search : SubstrateDomainSearches.values())
			if (search.type() == type)
				flag = true;
		return flag;
	}

	/**
	 * Get the chemical structure of the top-scoring BLASTp result for a fatty acyl-AMP ligase. By default, returns 12 carbons.
	 * @param domain	domain to search
	 * @return			top scoring structure, as a FattyAcidsEnum
	 */
	public static CStarterSubstrates cStarterType(Domain domain) {
		BlastSearchResult top = domain.blastResults().get(0);
		CStarterSubstrates type = null;
		
		for (CStarterSubstrates substrate : CStarterSubstrates.values()) {
			if (top.subject().contains(substrate.toString()))
				type = substrate;
		}
		
		return type;
	}

	public static String getDomainTypeNameForCSS(Domain domain) {
		String name = null;
		if (domain.family() == DomainFamilies.TAILORING) {
			name = "tailoring";
		} else if (domain.family() == DomainFamilies.TYPE_II_POLYKETIDE) {
			name = "type_ii_pks";
		} else {
			name = domain.type().toString();
		} 
		return name;
	}
	
}
