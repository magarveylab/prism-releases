package ca.mcmaster.magarveylab.prism.homology;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.prism.blast.BlastSearchResult;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.homology.data.HomologousCluster;
import ca.mcmaster.magarveylab.prism.util.Files;
import ca.mcmaster.magarveylab.prism.util.Sorter;
import ca.mcmaster.magarveylab.prism.util.Strings;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Search for homologous clusters.
 * @author skinnider
 *
 */
public class HomologousClusterSearch {

	private final String extension = ".phr";
	
	private Cluster cluster;
	private Session session;
	
	/**
	 * Instantiate a new database search for homologous clusters.
	 * @param cluster	query cluster
	 * @param session	current session
	 */
	public HomologousClusterSearch(Cluster cluster, Session session) {
		this.cluster = cluster;
		this.session = session;
	}
	
	/**
	 * Search PRISM's database for known biosynthetic gene clusters with homology to this cluster.
	 * @throws NumberFormatException 
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public void run() throws NumberFormatException, IOException, InterruptedException {
		instantiate();
		runIdentitySearch();
		runDomainSearch();
		
		Sorter.sortHomologousClustersByIdentityScore(cluster.homologs()); // do this first as a tiebreaker
		Sorter.sortHomologousClustersByAverageDomainScore(cluster.homologs());
		Sorter.sortHomologousClustersByDomainScore(cluster.homologs());
	}
	
	public void instantiate() {
		List<HomologousCluster> homologs = new ArrayList<HomologousCluster>();
		String dir = session.subDir("identityFasta");
		File[] files = Files.getDirectoryFiles(dir, extension);
		if (files != null)
			for (File file : files) {
				String database = file.getAbsolutePath().replace(extension, "");
				String name = Strings.name(database);
				HomologousCluster hc = new HomologousCluster(name);
				homologs.add(hc);
			}
		cluster.setHomologs(homologs);
	}
	
	public void runIdentitySearch() throws IOException, InterruptedException {
		String query = cluster.file("scaffoldOrfs");
		String dir = session.subDir("identityFasta");
		
		File[] files = Files.getDirectoryFiles(dir, extension);
		if (files != null)
			for (File file : files) {
				String database = file.getAbsolutePath().replace(extension, "");
				session.listener().updateLastDetail("Calculating cluster " + cluster.index() + " identity to " 
						+ Strings.name(database) + "...");
				IdentityScorer search = new IdentityScorer(database, query, cluster, session);
				search.run();
			}
	}
	
	public void runDomainSearch() throws IOException, InterruptedException {
		session.listener().updateLastDetail("Calculating cluster homology scores...");

		DomainScorer ds = new DomainScorer(cluster, session);
		List<BlastSearchResult> results = ds.score();

		for (HomologousCluster homolog : cluster.homologs())
			for (BlastSearchResult result : results) {
				String homologName = homolog.name();
				String resultName = result.subject();
				if (resultName.contains(homologName)) {
					homolog.setCoverage(1.0d * result.coverage());
					homolog.setDomainScore(result.score());
				}
			}
	}
	
	public static int getDomainSize(Cluster cluster) {
		int size = 0;
		for (Orf orf : cluster.orfs())
			for (Domain domain : orf.domains())
				size += domain.size();
		return size;
	}

}
