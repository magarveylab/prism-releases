package ca.mcmaster.magarveylab.prism.cluster;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ca.mcmaster.magarveylab.prism.cluster.analysis.ClusterAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.analysis.OrfAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.analysis.SugarAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.analysis.TypeIIPolyketideAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.module.ModuleFinder;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Contig;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.web.PrismConfig;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Group orfs into biosynthetic clusters. 
 * @author skinnider
 *
 */
public class ClusterFinder {

	private int counter = 1;

	/**
	 * Cluster biosynthetic orfs.
	 * Step 1: Get a list of orfs using a greedy algorithm
	 * Step 2: Detect modules, orf types within cluster
	 * Step 3: Determine if this is a cluster, and assign subtype.
	 * @param contig	a contig to cluster
	 * @param config	current PRISM configuration
	 * @param session	current session
	 * @throws IOException 
	 */
	public void cluster(Contig contig, PrismConfig config, Session session) throws IOException {
		List<Orf> orfs = contig.orfs();
		boolean[] clustered = new boolean[orfs.size()];
		Arrays.fill(clustered, false);
		List<Cluster> clusters = new ArrayList<Cluster>();

		System.out.println("[ClusterFinder] Clustering " + orfs.size() + " orfs with window " + config.window);

		int start, end;
		// For each gene, try to cluster with any genes within window and mark as such
		for (int i = 0; i < orfs.size(); i++) {
			Orf orf = orfs.get(i);
			start = orf.start();
			end = orf.end();
			if (orf.hasBiosyntheticDomains() && clustered[i] == false) {
				Cluster cluster = new Cluster();
				cluster.setIndex(counter);
				cluster.orfs().add(orf);
				clustered[i] = true;
				for (int j = 1; j < orfs.size(); j++) {
					Orf orf2 = orfs.get(j);
					if (orf2.hasDomains() && clustered[j] == false) {
						if (Math.abs(start - orf2.end()) < config.window || 
								Math.abs(orf2.start() - end) < config.window) {
							cluster.orfs().add(orf2);
							clustered[j] = true;
							if (orf2.hasBiosyntheticDomains()) {
								if (orf2.start() < start)
									start = orf2.start();
								if (orf2.end() > end)
									end = orf2.end();
							}
						}
					}
				}
				
				// detect all modules in the cluster
				for (Orf o : cluster.orfs()) 
					ModuleFinder.detectModules(o);
				
				// perform module/substrate checks
				ClusterAnalyzer.setExtendability(cluster);
				ClusterAnalyzer.checkPyrroleModules(cluster); 
				ClusterAnalyzer.checkTransAdenylationModules(cluster);
				ClusterAnalyzer.checkTransAcyltransferaseModules(cluster); 
				ClusterAnalyzer.checkStarterModules(cluster);
				ClusterAnalyzer.checkDockingModules(cluster);
				ClusterAnalyzer.checkSubstratePrerequisites(cluster);
				TypeIIPolyketideAnalyzer.checkC6Methyltransferase(cluster);
				TypeIIPolyketideAnalyzer.checkCyclaseClade6B(cluster);
				
				// detect each orf type, cluster type, and cluster frame
				for (Orf o : cluster.orfs())
					OrfAnalyzer.detectType(o);
				ClusterAnalyzer.detectFrame(cluster);
				
				// detect sugar types
				session.listener().updateLastDetail("Detecting possible sugar combinations...");
				cluster.setSugars(SugarAnalyzer.getSugars(cluster));

				// try to assign type
				ClusterAnalyzer.setClusterType(cluster);
				if (cluster.types().size() > 0) {
					clusters.add(cluster);
					counter++;
				} else {
					continue;
				}
			}
		}

		contig.clusters().addAll(clusters);
	}

}
