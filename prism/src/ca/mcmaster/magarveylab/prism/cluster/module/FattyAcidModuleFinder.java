package ca.mcmaster.magarveylab.prism.cluster.module;

import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.enums.CStarterSubstrates;
import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.enums.substrates.AdenylationSubstrates;
import ca.mcmaster.magarveylab.prism.cluster.analysis.DomainAnalyzer;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.Orf;

/**
 * Finds modules which add starter units or acyl adenylates. 
 * @author skinnider 
 *
 */
public class FattyAcidModuleFinder {

	/**
	 * Detect acyl-adenylating modules on their own open reading frame (i.e., not in conventional adenylation module
	 * arrangements).
	 * @param orf	orf to search
	 * @return		a list of acyl-adenylating modules. 
	 */
	public static List<Module> detectFattyAcidModules(Orf orf) {
		List<Module> modules = new ArrayList<Module>();
		if (orf.domains().size() == orf.domains(ThiotemplatedDomains.ACYL_ADENYLATING).size()) {
			Module module = new Module(ModuleTypes.ACYL_ADENYLATE);
			for (Domain domain : orf.domains(ThiotemplatedDomains.ACYL_ADENYLATING))
				module.add(domain);
			modules.add(module);
		}
		return modules;
	}

	/**
	 * Detect pyrrole ligases (prolyl-AMP ligases on their own orf).
	 * @param orf	orf to search
	 * @return		a list of prolyl-AMP ligases
	 */
	public static List<Module> findPyrroleModules(Orf orf) {
		List<Module> modules = new ArrayList<Module>();
		if (orf.domains().size() == orf.domains(ThiotemplatedDomains.ADENYLATION).size()) {
			Module module = new Module(ModuleTypes.ACYL_ADENYLATE);
			for (Domain domain : orf.domains(ThiotemplatedDomains.ADENYLATION))
				if (domain.topSubstrate().type() == AdenylationSubstrates.PROLINE_3)
					module.add(domain);
			if (module.domains().size() == 0)
				return modules;
			modules.add(module);
		}
		return modules;
	}
	
	/**
	 * Detect starter condensation domain modules which activate a fatty acid.
	 * @param orf	orf to search
	 * @return		a list of C-starter fatty acid modules
	 */
	public static List<Module> detectCStarterModules(Orf orf) {
		List<Module> modules = new ArrayList<Module>();
		for (Domain domain : orf.domains()) {
			if (DomainAnalyzer.isCStarter(domain)) {
				String result = domain.blastResults().get(0).subject();
				String[] names = CStarterSubstrates.names();
				for (String name : names) {
					if (result.indexOf(name) != -1) {
						Module module = new Module(ModuleTypes.C_STARTER);
						module.add(domain);
						modules.add(module);
					}
				}
			}
		}
		return modules;
	}

}
