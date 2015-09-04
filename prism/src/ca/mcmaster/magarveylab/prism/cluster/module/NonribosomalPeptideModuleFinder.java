package ca.mcmaster.magarveylab.prism.cluster.module;

import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.Orf;

/**
 * Finds biosynthetic modules within a nonribosomal peptide synthetase.
 * @author skinnider
 *
 */
public class NonribosomalPeptideModuleFinder {

	/**
	 * Detect adenylation and acyl-adenylating modules in an orf in the prototypical C-A-T arrangement. 
	 * Adenylation module detection works by finding a condensation (or epimerization)
	 * domain and the next thiolation domain, and including everything in between in the putative module. If the module 
	 * does not contain an adenylation domain, it is labelled a potential site of trans-adenylation insertion.
	 * @param orf	the orf to search
	 * @return		a list of adenylation modules
	 */
	public static List<Module> detectAdenylationModules(Orf orf) {
		List<Module> modules = new ArrayList<Module>();
		List<Domain> domains = orf.domains();

		int start = -1, end = -1;
		Domain adenylation = null, acylAdenylating = null;
		for (int i = 0; i < domains.size(); i++) {
			Domain domain = domains.get(i);
			if (domain.type() == ThiotemplatedDomains.CONDENSATION)
				start = i;
			if (start != -1) { // don't detect thiolation first
				if (domain.type() == ThiotemplatedDomains.THIOLATION 
						|| domain.type() == ThiotemplatedDomains.THIOESTERASE) {
					end = i;
					Module module = new Module(ModuleTypes.ADENYLATION); 
					for (int j = start; j <= end; j++) {
						Domain d = domains.get(j);
						module.add(d);
						if (d.type() == ThiotemplatedDomains.ADENYLATION)
							if (acylAdenylating == null || d.score() > acylAdenylating.score()) {
								adenylation = d;
								acylAdenylating = null;
							}
						if (d.type() == ThiotemplatedDomains.ACYL_ADENYLATING)
							if (adenylation == null || d.score() > adenylation.score()) {
								acylAdenylating = d;
								adenylation = null;
							}
					}
					if (acylAdenylating != null) {
						module.setType(ModuleTypes.ACYL_ADENYLATE);
						modules.add(module);
					} else if (adenylation != null) { 
						modules.add(module);
					} else {
						module.setType(ModuleTypes.TRANS_ADENYLATION_INSERTION);
						modules.add(module);
					}
					// clear flags
					start = -1;
					end = -1;
					adenylation = null;
					acylAdenylating = null;
				}
			}
		}
		return modules;
	}

	/**
	 * Detect potential trans-adenylation modules, which occur at the beginning
	 * of the orf (or immediately after a formylation domain), and lack a
	 * starting condensation domain. These modules must be checked later to
	 * ensure that they are the only module on their orf, and that a possible
	 * insertion site exists on the open reading frame; if they are not
	 * trans-adenylation they will be considered starter adenylation modules.
	 * 
	 * @param orf
	 *            the orf to search
	 * @return a list of potential trans-adenylation modules
	 */
	public static List<Module> detectTransAdenylationOrStarterAdenylationModules(Orf orf) {
		List<Module> modules = new ArrayList<Module>();
		List<Domain> domains = orf.domains();

		int start = -1, end = -1;
		int beginning = orf.domains(TailoringDomains.FORMYLTRANSFERASE).size() > 0 ? 1 : 0;
		Domain adenylation = null, acylAdenylating = null;
		for (int i = 0; i < domains.size(); i++) {
			Domain d = domains.get(i);
			if (d.type() == ThiotemplatedDomains.ADENYLATION && i == beginning)
				if (acylAdenylating == null || d.score() > acylAdenylating.score()) {
					start = i;
					adenylation = d;
					acylAdenylating = null;
				}
			if (d.type() == ThiotemplatedDomains.ACYL_ADENYLATING && i == beginning) 
				if (adenylation == null || d.score() > adenylation.score()) {
					start = i;
					acylAdenylating = d;
					adenylation = null;
				}
			if (start != -1) { 
				// don't detect condensation domain first
				if (d.type() == ThiotemplatedDomains.CONDENSATION) {
					start = -1;
					continue;
				}
				if (d.type() == ThiotemplatedDomains.THIOLATION) {
					end = i;
					Module module = null;
					if (acylAdenylating != null) {
						module = new Module(ModuleTypes.ACYL_ADENYLATE);
					} else if (adenylation != null) {
						module = new Module(ModuleTypes.TRANS_ADENYLATION);
					}
					if (module != null) {
						for (int j = start; j <= end; j++) {
							Domain domain = domains.get(j);
							module.add(domain);
						}
						modules.add(module);
					}
					// clear flags
					start = -1;
					end = -1;
					adenylation = null;
					acylAdenylating = null;
				}
			}
		}
		return modules;
	}
	
	/**
	 * Detect adenylation modules which end an orf (i.e., an orf which ends with C-A, without a thiolation domain). 
	 * @param orf	orf to search
	 * @return		a list of adenylation end modules
	 */
	public static List<Module> detectAdenylationEndModules(Orf orf) {
		List<Module> modules = new ArrayList<Module>();
		List<Domain> domains = orf.domains();
	
		if (domains.size() >= 2) {
			int lastIdx = domains.size() - 1;
			int secondLastIdx = domains.size() - 2;
			Domain last = domains.get(lastIdx);
			Domain secondLast = domains.get(secondLastIdx);
			if (last.type() == ThiotemplatedDomains.ADENYLATION && secondLast.type() == ThiotemplatedDomains.CONDENSATION) {
				Module module = new Module(ModuleTypes.ADENYLATION);
				module.add(secondLast);
				module.add(last);
				modules.add(module);
			}
		}
		return modules;
	}
	
}
