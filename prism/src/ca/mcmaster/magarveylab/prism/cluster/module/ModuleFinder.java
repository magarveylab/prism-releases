package ca.mcmaster.magarveylab.prism.cluster.module;

import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.util.Sorter;

/**
 * Detects biosynthetic modules.
 * @author skinnider
 *
 */
public class ModuleFinder {
	
	/**
	 * Detect all modules in an orf.
	 * @param orf	the orf to search
	 * @return		all biosynthetic modules
	 */
	public static void detectModules(Orf orf) {
		List<Module> modules = new ArrayList<Module>();
				
		List<Module> adenylation = NonribosomalPeptideModuleFinder.detectAdenylationModules(orf);
		List<Module> transAndStarter = NonribosomalPeptideModuleFinder.detectTransAdenylationOrStarterAdenylationModules(orf);
		List<Module> endAdenylation = NonribosomalPeptideModuleFinder.detectAdenylationEndModules(orf);
		
		List<Module> acyltransferase = PolyketideModuleFinder.detectAcyltransferaseModules(orf);
		List<Module> starterAcyltransferase = PolyketideModuleFinder.detectAcyltransferaseStarterModules(orf);
		List<Module> endAcyltransferase = PolyketideModuleFinder.detectAcyltransferaseEndModules(orf);
		List<Module> transAT = PolyketideModuleFinder.detectTransATModules(orf);
		
		List<Module> fattyAcid = FattyAcidModuleFinder.detectFattyAcidModules(orf);
		List<Module> pyrrole = FattyAcidModuleFinder.findPyrroleModules(orf);
		List<Module> cstarter = FattyAcidModuleFinder.detectCStarterModules(orf);
		
		modules.addAll(adenylation);
		modules.addAll(transAndStarter);
		modules.addAll(endAdenylation);
		modules.addAll(acyltransferase);
		modules.addAll(starterAcyltransferase);
		modules.addAll(endAcyltransferase);
		modules.addAll(transAT);
		modules.addAll(fattyAcid);
		modules.addAll(pyrrole);
		modules.addAll(cstarter);
		
		Sorter.sortModules(modules);
		
		// command line output
		for (Module module : modules) 
			System.out.println("[ModuleFinder] Detected " + module.type() + " module in orf " + orf.name() + " with "
					+ module.domains().size() + " domains");
		if (modules.size() == 0)
			System.out.println("[ModuleFinder] No biosynthetic modules found on " + orf.name());
		
		orf.setModules(modules);
	}
	
}
