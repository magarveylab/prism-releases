package ca.mcmaster.magarveylab.prism.cluster.module;

import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.Orf;

/**
 * Finds biosynthetic modules within a polyketide synthase.
 * @author skinnider
 *
 */
public class PolyketideModuleFinder {

	/**
	 * Detect acyltransferase modules in an orf. Acyltransferase module detection works by finding a ketosynthase domain
	 * and the next thiolation domain, and including everything in between in the putative module. If the module does not
	 * contain an acyltransferase domain, it is rejected. 
	 * @param orf	the orf to search
	 * @return		a list of acyltransferase modules
	 */
	public static List<Module> detectAcyltransferaseModules(Orf orf) {
		List<Module> modules = new ArrayList<Module>();
		List<Domain> domains = orf.domains();
	
		int start = -1, end = -1;
		boolean acyltransferase = false, adenylation = false;
		for (int i = 0; i < domains.size(); i++) {
			Domain domain = domains.get(i);
			if (domain.type() == ThiotemplatedDomains.KETOSYNTHASE)
				start = i;
			if (start != -1) { // don't detect thiolation first
				if (domain.type() == ThiotemplatedDomains.THIOLATION || domain.type() == ThiotemplatedDomains.THIOESTERASE) {
					end = i;
					Module module = new Module(ModuleTypes.ACYLTRANSFERASE);
					for (int j = start; j <= end; j++) {
						Domain d = domains.get(j);
						if (d.type() == ThiotemplatedDomains.ACYLTRANSFERASE)
							acyltransferase = true;
						if (d.type() == ThiotemplatedDomains.ADENYLATION)
							adenylation = true;
						module.add(d);
					}
					if (acyltransferase) {
						modules.add(module);
					} else if (!adenylation) {
						module.setType(ModuleTypes.TRANS_AT_INSERTION);
						modules.add(module);
					}
					// clear flags
					start = -1;
					end = -1;
					acyltransferase = false;
				}
			}
		}
		return modules;
	}

	/**
	 * Detect acyltransferase modules which begin an orf (i.e., which lack a KS domain).
	 * @param orf	orf to search
	 * @return		list of AT start modules
	 */
	public static List<Module> detectAcyltransferaseStarterModules(Orf orf) {
		List<Module> modules = new ArrayList<Module>();
		List<Domain> domains = orf.domains();
		
		if (domains.size() > 0 && domains.get(0).type() == ThiotemplatedDomains.ACYLTRANSFERASE) {
			for (int i = 0; i < domains.size(); i++) {
				Domain domain = domains.get(i);
				if (domain.type() == ThiotemplatedDomains.THIOLATION) {
					Module module = new Module(ModuleTypes.ACYLTRANSFERASE);
					for (int j = 0; j <= i; j++) {
						Domain d = domains.get(j);
						module.add(d);
					}
					modules.add(module);
					break;
				}
			}
		}
		
		return modules;
	}

	/**
	 * Detect acyltransferase modules which end an orf (i.e., an orf which ends with KS-(x)-AT, without a thiolation domain). 
	 * @param orf	orf to search
	 * @return		a list of AT end modules
	 */
	public static List<Module> detectAcyltransferaseEndModules(Orf orf) {
		List<Module> modules = new ArrayList<Module>();
		List<Domain> domains = orf.domains();
	
		if (domains.size() > 0) {
			// get last KS
			int KSIdx = -1;
			for (int i = 0; i < domains.size() - 1; i++) {
				Domain domain = domains.get(i);
				if (domain.type() == ThiotemplatedDomains.KETOSYNTHASE)
					KSIdx = i;
			}
			
			int lastIdx = domains.size() - 1;
			Domain last = domains.get(lastIdx);
			if (last.type() == ThiotemplatedDomains.ACYLTRANSFERASE) {
				Module module = new Module(ModuleTypes.ACYLTRANSFERASE);
				if (KSIdx != -1) { // if there is a ketosynthase on the orf before this AT
					for (int j = lastIdx; j >= KSIdx; j--) {
						Domain d = domains.get(j);
						module.add(d);
					}
					modules.add(module);
				}
			}
		}
		
		return modules;
	}
	
	/**
	 * Detect trans-AT modules on an open reading frame. These AT domains must be the only domain on their orf.
	 * @param orf	orf to search
	 * @return		a list of trans-AT modules
	 */
	public static List<Module> detectTransATModules(Orf orf) {
		List<Module> modules = new ArrayList<Module>();
		if (orf.domains().size() == orf.domains(ThiotemplatedDomains.ACYLTRANSFERASE).size()) {
			Module module = new Module(ModuleTypes.TRANS_AT);
			for (Domain domain : orf.domains())
				module.add(domain);
			modules.add(module);
		}
		return modules;
	}
	
}
