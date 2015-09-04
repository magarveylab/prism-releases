package ca.mcmaster.magarveylab.prism.cluster.scaffold;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.exception.CDKException;

import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.domains.TypeIIPolyketideDomains;
import ca.mcmaster.magarveylab.prism.cluster.analysis.TypeIIPolyketideAnalyzer;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;

public class TypeIIPolyketideGenerator {
	
	/**
	 * Generate a linear chain of unreduced ketide units correponding to the
	 * minimal type II PKS scaffold, and set the starter unit. Ketide units are
	 * referred to by their distance from the thioester/carboxyl group, in order
	 * to be consistent with a wealth of biosynthetic research--i.e., in the
	 * opposite manner as scaffolds are constructed (so residue 1 is closest to
	 * the carboxy terminus). In order to do this, the module permutation is
	 * reversed before and after constructing the scaffold.
	 * 
	 * @param modules
	 *            list of type II PKS modules
	 * @param cluster
	 *            cluster to analyze
	 * @return linear ketide scaffold with starter unit
	 * @throws IOException
	 * @throws NoResidueException
	 * @throws ScaffoldGenerationException
	 * @throws CDKException 
	 */
	public static Scaffold generateTypeIIPolyketideScaffold(List<Module> modules, Cluster cluster)
			throws IOException, NoResidueException, ScaffoldGenerationException, CDKException {
		Collections.reverse(modules);
		Scaffold scaffold = StructureGenerator.generateLinearScaffold(modules);
		reverseScaffoldResidues(scaffold);
		Collections.reverse(modules);
		return scaffold;
	}

	/**
	 * Generate a list of type II PKS modules whose size corresponds to the
	 * chain length factor of the PKS.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return list of generated type II PKS modules
	 * @throws NoResidueException
	 */
	public static List<List<Module>> generateTypeIIPolyketideModulePermutation(Cluster cluster) 
			throws NoResidueException {
		List<List<Module>> permutation = new ArrayList<List<Module>>();
		if (TypeIIPolyketideAnalyzer.isTypeIIPolyketideCluster(cluster)) {
			List<Module> modules = new ArrayList<Module>();
			int chainLength = TypeIIPolyketideAnalyzer.getChainLength(cluster);
			for (int i = 0; i < chainLength; i++) {
				Module module = new Module(ModuleTypes.TYPE_II_PKS);
				modules.add(module);
			}
			Module starter = createStarterModule(cluster);
			modules.add(starter);
						
			permutation.add(modules);
			System.out.println("[TypeIIPolyketideGenerator] Generated " 
					+ modules.size() + "-module type II polyketide permutation");
		}
		return permutation;
	}

	/**
	 * Reverse the order in which residues access acetate units from a type II
	 * polyketide scaffold, for conceptual simplicity and consistency with
	 * biosynthesis papers.
	 * 
	 * @param scaffold
	 *            type II PKS scaffold
	 */
	public static void reverseScaffoldResidues(Scaffold scaffold) {
		Map<Module, Residue> newResidues = new LinkedHashMap<Module, Residue>();
		Map<Module, Residue> residues = scaffold.residues();
		Object[] modules = residues.keySet().toArray();
		for (int i = modules.length - 1; i >= 0; i--) {
			Module module = (Module) modules[i];
			Residue residue = residues.get(module);
			newResidues.put(module, residue);
		}
		scaffold.setResidues(newResidues);
	}
	
	/**
	 * Configure the first module in a list to a type II PKS-specific starter
	 * module by adding to the module all domains that imply the presence of a
	 * specific, non-acetate stater unit.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return the polyketide starter module
	 */
	public static Module createStarterModule(Cluster cluster) {
		Module starter = new Module(ModuleTypes.TYPE_II_PKS_STARTER);
		starter.add(cluster.domains(TypeIIPolyketideDomains.PRIMING_AT));
		starter.add(cluster.domains(TypeIIPolyketideDomains.KSIII));
		starter.add(cluster.domains(TypeIIPolyketideDomains.AMIDOTRANSFERASE));
		starter.add(cluster.domains(TypeIIPolyketideDomains.PAL));
		return starter;
	}
	
}
