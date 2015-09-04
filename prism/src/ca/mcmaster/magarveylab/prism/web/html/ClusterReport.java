package ca.mcmaster.magarveylab.prism.web.html;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.openscience.cdk.exception.CDKException;

import ca.mcmaster.magarveylab.prism.Prism;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.CombinatorialData;
import ca.mcmaster.magarveylab.prism.data.Genome;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.util.PrismStringBuffer;
import ca.mcmaster.magarveylab.prism.web.PrismConfig;
import ca.mcmaster.magarveylab.prism.web.html.graph.CircularGenomeGraph;
import ca.mcmaster.magarveylab.prism.web.html.graph.ClusterAssemblyGraph;
import ca.mcmaster.magarveylab.prism.web.html.graph.ClusterGraph;
import ca.mcmaster.magarveylab.prism.web.html.graph.Legend;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Get the HTML output for a single biosynthetic gene cluster.  
 * @author skinnider
 *
 */
public class ClusterReport {
	
	private Session session;
	private Prism prism;
	private Cluster cluster;
	private PrismReport report;
	private PrismConfig config;
	
	/**
	 * Instantiate a new cluster HTML report.
	 * @param cluster	cluster to generate HTML for
	 * @param prism		parent PRISM search 
	 */
	public ClusterReport(Cluster cluster, Prism prism) {
		this.cluster = cluster;
		this.prism = prism;
		this.config = prism.config();
		this.session = prism.session();
		this.report = (PrismReport) session.report();
	}
	
	/**
	 * Write a cluster report to a HTML file.
	 * @throws IOException
	 */
	public void write() throws IOException {
		File html = new File(session.dir() + "cluster_" + cluster.index() + ".html");
		if (!html.exists())
			html.createNewFile();
		FileWriter fw = new FileWriter(html);
		BufferedWriter bw = new BufferedWriter(fw);
		
		StringBuffer sb = new StringBuffer();
		PrismStringBuffer psb = new PrismStringBuffer(sb);
		
		String scaffoldFastaLink = "<a href='" + session.webDir() + "cluster_" + cluster.index() + ".fasta'>FASTA</a>";
		String clusterFastaLink = "<a href='" + session.webDir() + "cluster_" + cluster.index() + "_full.fasta'>FASTA</a>";
		String scaffoldLibraryLink = "<a href='" + session.webDir() + "cluster_" + cluster.index() + "_library.txt'>TXT</a>";
		String isnapLibraryLink = "<a href='" + session.webDir() + "cluster_" + cluster.index() + "_isnap.txt'>TXT</a>";
		String clusterGenomicLink = "<a href='" + session.webDir() + "cluster_" + cluster.index() + "_genomic.fasta'>FASTA</a>";
		String clusterJsonLink = "<a href='" + session.webDir() + "cluster_" + cluster.index() + ".json'>JSON</a>";

		try {			
			// write header
			String header = report.writeClusterHeader();
			psb.append(header);
			psb.appendLine("<div class='container'>");
			 
			// embed highlighted genome graph
			Genome genome = prism.genome();
			if (genome.contigs().size() == 1) {
				CircularGenomeGraph highlightedGenomeGraph = cluster.graph();
				if (cluster.graph() != null) {
					String highlightedGenomeGraphHtml = highlightedGenomeGraph.html();
					psb.append(highlightedGenomeGraphHtml); 
				} else {
					System.out.println("[ClusterReport] Error: no genome graph object");
				}
			}
			
			psb.appendLine("<header class='sub'>");
			psb.appendLine("<h3>Cluster " + cluster.index() + "</h3>");
			psb.appendLine("<p>From search #" + session.id());
			if (prism.genome().filename() != null)
				psb.appendLine("<br>Sequence: " + prism.genome().filename());
//			psb.appendLine("<br>Fasta header: " + genome.organism().rawfastaheader());
//			psb.appendLine("<br>Accession: " + genome.organism().accession());
//			psb.appendLine("<br>Genus: " + genome.organism().genus());
//			psb.appendLine("<br>Species: " + genome.organism().species());
//			psb.appendLine("<br>Strain: " + genome.organism().strain());
			psb.appendLine("</p>");
			psb.appendLine("</header>");

			psb.appendLine("<section class='cluster cf'>");
			psb.appendLine("<div class='clusterHeader'>");

			psb.appendLine("<h3>Biosynthetic assembly:</h3>");
			String assembly = ClusterAssemblyGraph.html(cluster);
			psb.append(assembly);

			psb.appendLine("<h3>Cluster:</h3>");
			String clusterGraph = ClusterGraph.html(cluster);
			psb.append(clusterGraph);
			
			psb.appendLine("<h3>Legend:</h3>");
			String legend = Legend.getOrfGraphLegend(config);
			psb.appendLine(legend);

			psb.appendLine("</div>");

			psb.appendLine("<section class='reportHeaderRight'>");
			psb.appendLine("<h3>Predicted cluster product:</h3>");
			if (cluster.library() != null 
					&& cluster.library().size() > 0) {
				psb.appendLine("<p>Combinatorial structure library: " + scaffoldLibraryLink + "<br>");
				psb.appendLine("iSNAP library: " + isnapLibraryLink + "</p>");
				
				CombinatorialData cd = cluster.combinatorialData();
				if (cd != null) {
					psb.appendLine("<h3>Combinatorial data evaluated:</h3>");
					psb.appendLine("<ul>");
					int orfPermutations = cd.getNumOrfPermutations();
					String orfPermutationsSize = orfPermutations >= 500 ? "500+"
							: orfPermutations + "";
					psb.appendLine("<li>Open reading frame permutations: " 
							+ orfPermutationsSize + "<br>(maximum 500)</li>");
					int sugars = cd.getNumSugars();
					String sugarSize = sugars >= 100 ? "100+" : sugars + "";
					psb.appendLine("<li>Sugar combinations: " 
							+ sugarSize + " (maximum 100)</li>");
					psb.appendLine("<li>Cyclization patterns: " 
							+ cd.getNumCyclizations() + "</li>");
					psb.appendLine("<li>Tailoring reaction plans: " 
							+ cd.getNumReactions() + "</li>");
					int combinatorialPlans = cd.getNumCombinatorialPlans();
					String combinatorialPlanSize = combinatorialPlans == 1001 ? "1,000+"
							: combinatorialPlans + "";
					psb.appendLine("<li>Total combinatorial plans: " 
							+ combinatorialPlanSize + " (maximum 1,000)</li>");
					int librarySize = cluster.library().size();
					psb.appendLine("<li>Scaffolds generated: "
							+ librarySize + " (maximum " + config.scaffoldLimit
							+ ")</li>");
					psb.appendLine("</ul>");
				}
			} else {
				psb.appendLine("<strong>Error:</strong> Unable to generate predicted structure!");
			}
			
			psb.appendLine("<h3>Downloads:</h3>");
			if (cluster.file("sequence") != null) 
				psb.appendLine("<p>Genomic cluster sequence: " + clusterGenomicLink + "<br>");
			if (cluster.file("scaffoldOrfs") != null)
				psb.appendLine("Biosynthetic open reading frames: " + scaffoldFastaLink + "<br>");
			if (cluster.file("clusterOrfs") != null)
				psb.appendLine("All cluster open reading frames: " + clusterFastaLink + "<br>");
			if (cluster.file("clusterJson") != null)
				psb.appendLine("JSON report: " + clusterJsonLink + "<br>");
			psb.appendLine("</p>");
			psb.appendLine("</section>");

			if (cluster.sugars().size() > 0)
				psb.append(SugarReport.getHTML(cluster));
			
			if (config.score)
				psb.append(DereplicationReport.getHTML(cluster, session));
			psb.appendLine("</section>");

			psb.appendLine("<header class='sub cf'>");
			psb.appendLine("<h3>Analysis</h3>");
			psb.appendLine("</header>");
			
			for (Orf orf : cluster.orfs()) {
				String orfHtml = OrfReport.getHTML(orf, session);
				psb.appendLine(orfHtml);
			}
			
			psb.appendLine("</div>"); 
			String footer = report.writeFooter();
			psb.append(footer);
		} catch (CDKException e) {
			session.exceptionHandler().throwException(e);
		}

		bw.write(psb.toString());
		bw.close();
	}

}
