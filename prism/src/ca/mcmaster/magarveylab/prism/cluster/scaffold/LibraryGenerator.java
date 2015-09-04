package ca.mcmaster.magarveylab.prism.cluster.scaffold;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import org.openscience.cdk.exception.CDKException;

import ca.mcmaster.magarveylab.prism.Prism;
import ca.mcmaster.magarveylab.prism.combinatorialization.ClusterPermuter;
import ca.mcmaster.magarveylab.prism.combinatorialization.CombinatorialPlanGenerator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Library;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.structure.CombinatorialPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.PrismFileWriter;
import ca.mcmaster.magarveylab.prism.util.exception.ClassInstantiationException;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.web.PrismConfig;
import ca.mcmaster.magarveylab.wasp.session.Session;
import ca.mcmaster.magarveylab.prism.util.SmilesIO;

/**
 * Generates combinatorial scaffold libraries.
 * 
 * @author skinnider
 *
 */
public class LibraryGenerator {

	private Library library = new Library();
	private Cluster cluster;
	private Session session;
	private PrismConfig config;

	/**
	 * Instantiate a new scaffold library generator for a single cluster.
	 * 
	 * @param cluster
	 *            cluster for which to generate scaffold library
	 * @param session
	 *            current session
	 */
	public LibraryGenerator(Cluster cluster, Session session) {
		this.cluster = cluster;
		this.session = session;
		Prism prism = (Prism) session.webapp();
		this.config = prism.config();
	}

	/**
	 * Generate a combinatorial libary of cluster products.
	 * @return			the combinatorial library
	 * @throws NoResidueException 
	 * @throws IOException 
	 * @throws ClassInstantiationException 
	 * @throws CDKException 
	 * @throws ClassNotFoundException 
	 */
	public void generate() throws IOException, NoResidueException, ClassInstantiationException, CDKException {
		session.listener().updateLastDetail("Generating predicted scaffolds for cluster " + cluster.index() + "...");

		if (config.scaffoldLimit <= 0)
			return;
		
		List<List<Module>> permutations = ClusterPermuter.getModulePermutations(cluster);		
		List<CombinatorialPlan> plans = CombinatorialPlanGenerator.generateAllPlans(cluster, permutations);
		combinatorializeScaffold(plans, config);

		cluster.setLibrary(library);
		write(library, cluster, session);
	}

	/**
	 * Generate a combinatorial library of chemical structures for each generated permutation of biosynthetic modules.	
	 * @param modulePermutations		list of module permutations
	 */
	private void combinatorializeScaffold(List<CombinatorialPlan> plans, PrismConfig config) {
		boolean[] generated = new boolean[plans.size()];
		Arrays.fill(generated, false);

		Random random = new Random(0);
		while (library.size() < config.scaffoldLimit) {
			int i = random.nextInt(plans.size());
			if (!areRemainingPlans(generated)) {
				System.out.println("[Library] Exhausted all combinatorial plans, exiting scaffold generation loop");
				break;
			}
			if (generated[i])
				continue;
			
			session.listener().updateLastDetail("Generating scaffold " + (library.size()+1) + " of a maximum " 
					+ config.scaffoldLimit + " for cluster " + cluster.index() + "...");
			System.out.println("[LibraryGenerator] Executing combinatorial plan " + i);
			
			try {
				CombinatorialPlan scheme = plans.get(i);
				Scaffold scaffold = StructureGenerator.generateScaffold(cluster, scheme);
				String smiles = SmilesIO.smiles(scaffold.molecule());
				System.out.println("[LibraryGenerator] Generated scaffold with SMILES " + smiles);

				if (!library.contains(smiles)) {
					if (scaffold.reactionCount() > library.reactionCount()) {
						library.clear();
						library.add(smiles);
						library.setReactionCount(scaffold.reactionCount());
					} else if (scaffold.reactionCount() == library.reactionCount()) {
						library.add(smiles);
					}
				}
			} catch (Exception e) {
				System.out.println(e.getMessage());
				e.printStackTrace();
			}

			generated[i] = true;
		}
		
		System.out.println("[LibraryGenerator] Generated new library of " + library.size() + " scaffolds"
				+ " with reaction count " + library.reactionCount());
	}

	/**
	 * Write the combinatorial scaffold library as a library file and an iSNAP
	 * file.
	 * 
	 * @throws IOException
	 */
	public static void write(Library library, Cluster cluster, Session session) throws IOException {
		String libraryPath = session.dir() + "cluster_" + cluster.index() + "_library.txt";
		String isnapPath = session.dir() + "cluster_" + cluster.index() + "_isnap.txt";
		PrismFileWriter.writeLibraryFile(library, libraryPath);
		PrismFileWriter.writeiSNAPFile(cluster, isnapPath);
		cluster.addFile("library", libraryPath);
		cluster.addFile("isnap", isnapPath);
	}
	
	/**
	 * Determine whether there are combinatorial plans remaining.
	 * 
	 * @param generated
	 *            boolean array of generated plans
	 * @return true if there are any un-executed plans remaining
	 */
	public boolean areRemainingPlans(boolean[] generated) {
		boolean flag = false;
		for (boolean b : generated)
			if (b == false)
				flag = true;
		return flag;
	}
	
}
