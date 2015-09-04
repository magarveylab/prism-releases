package ca.mcmaster.magarveylab.prism.web.html.graph;

import java.io.IOException;

import ca.mcmaster.magarveylab.prism.homology.data.HomologousCluster;
import ca.mcmaster.magarveylab.prism.util.Files;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * A graph of all the biosynthetic open reading frames in a homologous natural product gene cluster.
 * @author skinnider
 *
 */
public class ClusterHomologGraph {
	
	/**
	 * Get the cluster graph for this homologous cluster by reading the appropriate HTML text file. 
	 * @param homolog	homologous cluster to graph
	 * @param session	current session
	 * @return			cluster graph as HTML string 
	 * @throws IOException
	 */
	public static String html(HomologousCluster homolog, Session session) throws IOException {
		String name = homolog.name();
		String file = session.subDir("clusterGraphs") + name + ".txt";
		String contents = Files.readFile(file);
		return contents;
	}
	
}
