package ca.mcmaster.magarveylab.prism.cluster.analysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.enums.CStarterSubstrates;
import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.enums.substrates.AdenylationSubstrates;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Atoms;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.Substrate;
import ca.mcmaster.magarveylab.prism.util.SmilesIO;

/**
 * Utilities class for operations on biosynthetic modules.
 * @author skinnider
 *
 */
public class ModuleAnalyzer {
	
	/**
	 * Get all modules which contain a hydroxyl group, including cyclization points, hydroxyl-containing amino acids, and 
	 * modules followed by a module with a KR domain and no DH.
	 * @param permutation	modules to analyze
	 * @return				all hydroxyl-containing modules
	 * @throws IOException 
	 * @throws CDKException 
	 */
	public static List<Module> hydroxyls(List<Module> permutation) throws IOException, CDKException {
		List<Module> hydroxyls = new ArrayList<Module>();
		
		// get cyclization -OH
		List<Module> cyclizations = cyclizationHydroxyls(permutation);
		hydroxyls.addAll(cyclizations);

		// get -OH containing amino acids
		hydroxyls.addAll(ModuleAnalyzer.modules(permutation, AdenylationSubstrates._3_5_DIHYDROXYPHENYLGLYCINE));
		hydroxyls.addAll(ModuleAnalyzer.modules(permutation, AdenylationSubstrates._4_HYDROXY_PHENYLGLYCINE));
		hydroxyls.addAll(ModuleAnalyzer.modules(permutation, AdenylationSubstrates.TYROSINE_1));
		hydroxyls.addAll(ModuleAnalyzer.modules(permutation, AdenylationSubstrates.TYROSINE_2));
		hydroxyls.addAll(ModuleAnalyzer.modules(permutation, AdenylationSubstrates.BETA_HYDROXYTYROSINE));
		hydroxyls.addAll(ModuleAnalyzer.modules(permutation, AdenylationSubstrates.HYDROXYVALINE));
		
		return hydroxyls;
	}

	/**
	 * Get all cyclization points which are also hydroxyl groups.
	 * @param permutation		module permutation to analyze
	 * @return					potential cyclization hydroxyls
	 * @throws IOException
	 * @throws CDKException 
	 */
	public static List<Module> cyclizationHydroxyls(List<Module> permutation) throws IOException, CDKException {
		List<Module> cyclizations = new ArrayList<Module>();
		
		cyclizations.addAll(ModuleAnalyzer.modules(permutation, AdenylationSubstrates.SERINE));
		cyclizations.addAll(ModuleAnalyzer.modules(permutation, AdenylationSubstrates.THREONINE_1));
		cyclizations.addAll(ModuleAnalyzer.modules(permutation, AdenylationSubstrates.THREONINE_2));
		cyclizations.addAll(ModuleAnalyzer.modules(permutation, AdenylationSubstrates.ALLO_THREONINE));
		
		// Serine/threonine modules cannot be cyclized already 
		Iterator<Module> itr = cyclizations.iterator();
		while (itr.hasNext()) {
			Module next = itr.next();
			for (Domain domain : next.domains())
				if (DomainAnalyzer.isCyclization(domain))
					itr.remove();
		}
		
		// next module has KR & no DH
		for (int i = 0; i < permutation.size() - 1; i++) {
			Module module = permutation.get(i);
			Module next = permutation.get(i+1);
			if (!module.contains(ThiotemplatedDomains.O_METHYLTRANSFERASE) && 
					next.contains(ThiotemplatedDomains.KETOREDUCTASE) && !next.contains(ThiotemplatedDomains.DEHYDRATASE) 
					&& !next.contains(ThiotemplatedDomains.ENOYLREDUCTASE) && cyclizations.indexOf(module) == -1)
				cyclizations.add(module);
		}

		// get C-starter/FAAL with free -OH, -NH2
		moduleLoop:
		for (int i = 0; i < permutation.size() - 1; i++) {
			Module module = permutation.get(i);
			Module next = permutation.get(i+1);
			if (module.type() == ModuleTypes.C_STARTER || module.type() == ModuleTypes.ACYL_ADENYLATE) {
				for (Domain domain : module.domains()) {
					String smiles = null;
					if (domain.type() == ThiotemplatedDomains.ACYL_ADENYLATING) {
						Substrate top = domain.topSubstrate();
						smiles = top.smiles();
					} else if (DomainAnalyzer.isCStarter(domain)) {
						CStarterSubstrates type = DomainAnalyzer.cStarterType(domain);
						smiles = type.smiles();
					}
					if (smiles != null) {
						IAtomContainer structure = SmilesIO.molecule(smiles);
						IAtom oxygen = Atoms.getHydroxyl(structure);
						if ((oxygen != null && structure.getAtomNumber(oxygen) != -1)
								|| next.contains(ThiotemplatedDomains.KETOREDUCTASE)) {
							cyclizations.add(module);
							continue moduleLoop;
						}
						IAtom nitrogen = Atoms.getNitrogen(structure);
						if (nitrogen!= null && structure.getAtomNumber(nitrogen) != -1
								&& structure.getBondOrderSum(nitrogen) == 1) {
							cyclizations.add(module);
							continue moduleLoop;
						}
					}
				}
			}	
		}
		
		// remove cyclization points less than/equal to 3 modules from end
		itr = cyclizations.iterator();
		while (itr.hasNext()) {
			Module next = itr.next();
			int idx = permutation.indexOf(next);
			int size = permutation.size();
			if (size - idx < 4)
				itr.remove();
		}
		
		return cyclizations;
	}
	
	/**
	 * Get all branched-chain amino and alpha-keto acids (isoleucine, leucine, valine).
	 * @param permutation	modules to analyze	
	 * @return				all branched-chain amino/alpha-keto acid modules
	 */
	public static List<Module> bcaa(List<Module> permutation) {
		List<Module> bcaa = new ArrayList<Module>();
		bcaa.addAll(modules(permutation, AdenylationSubstrates.ISOLEUCINE));
		bcaa.addAll(modules(permutation, AdenylationSubstrates.LEUCINE_1));
		bcaa.addAll(modules(permutation, AdenylationSubstrates.LEUCINE_2));
		bcaa.addAll(modules(permutation, AdenylationSubstrates.LEUCINE_3));
		bcaa.addAll(modules(permutation, AdenylationSubstrates.VALINE_1));
		bcaa.addAll(modules(permutation, AdenylationSubstrates.VALINE_2));
		bcaa.addAll(modules(permutation, AdenylationSubstrates.VALINE_3));
		return bcaa;
	}
	
	/**
	 * Get all six-membered aromatic amino acids (tyrosine, phenylalanine, and phenylglycine).
	 * @param permutation	modules to analyze
	 * @return				all 6-membered amino acid modules
	 */
	public static List<Module> sixMemberedAromatic(List<Module> permutation) {
		List<Module> aromatic = new ArrayList<Module>();
		aromatic.addAll(modules(permutation, AdenylationSubstrates.PHENYLALANINE));
		aromatic.addAll(modules(permutation, AdenylationSubstrates.BETA_HYDROXY_PHENYLALANINE));
		aromatic.addAll(modules(permutation, AdenylationSubstrates.BETA_METHYL_PHENYLALANINE));
		aromatic.addAll(modules(permutation, AdenylationSubstrates._3_5_DIHYDROXYPHENYLGLYCINE));
		aromatic.addAll(modules(permutation, AdenylationSubstrates._4_HYDROXY_PHENYLGLYCINE));
		aromatic.addAll(modules(permutation, AdenylationSubstrates.TYROSINE_1));
		aromatic.addAll(modules(permutation, AdenylationSubstrates.TYROSINE_2));
		aromatic.addAll(modules(permutation, AdenylationSubstrates.BETA_HYDROXYTYROSINE));
		return aromatic;
	}

	/**
	 * Get adenylation modules within a module permutation with a given adenylation domain substrate.
	 * @param permutation	modules to analyze
	 * @param type			substrate type to find
	 * @return				
	 */
	public static List<Module> modules(List<Module> permutation, AdenylationSubstrates type) {
		List<Module> modules = new ArrayList<Module>();
		for (Module module : permutation) 
			if (module.isAdenylationModule()) {
				Domain scaffold = module.scaffold();
				if (scaffold != null) {
					Substrate substrate = scaffold.topSubstrate();
						if (substrate.type() == type)
							modules.add(module);
					}
				}
		return modules;
	}

}
