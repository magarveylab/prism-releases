package ca.mcmaster.magarveylab.prism.web.html;

import java.io.IOException;
import java.util.List;

import org.openscience.cdk.exception.CDKException;

import ca.mcmaster.magarveylab.prism.Prism;
import ca.mcmaster.magarveylab.prism.cluster.analysis.ClusterAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.analysis.type.ClusterTypeAnalyzer;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Contig;
import ca.mcmaster.magarveylab.prism.data.Genome;
import ca.mcmaster.magarveylab.prism.util.PrismStringBuffer;
import ca.mcmaster.magarveylab.prism.web.PrismConfig;
import ca.mcmaster.magarveylab.prism.web.html.graph.CircularGenomeGraph;
import ca.mcmaster.magarveylab.prism.web.html.graph.ClusterGraph;
import ca.mcmaster.magarveylab.wasp.report.BasicReport;
import ca.mcmaster.magarveylab.wasp.report.Report;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Generates Prism HTML search reports.
 * 
 * @author skinnider
 *
 */
public class PrismReport extends BasicReport implements Report {
		
    private String context = "/prism";

	public PrismReport(Session session) {
		super(session);
	}
	
	public void writeClusterPages() throws IOException {
		if (session != null 
				&& session.webapp() != null 
				&& !session.webapp().isTerminated()) {
			Prism prism = (Prism) session.webapp();
			Genome genome = prism.genome();
			if (genome == null)
				return;
			for (Contig contig : genome.contigs()) {
				for (Cluster cluster : contig.clusters()) {
					session.listener().updateLastDetail("Writing report for cluster " 
							+ cluster.index() + " report...");
					ClusterReport cr = new ClusterReport(cluster, contig, prism);
					cr.write();
				}
			}
		}
	}
	
	/**
	 * Get the HTML for the Prism search report.
	 * 
	 * @throws IOException
	 * @throws CDKException
	 */
	public String getReportHtml() throws IOException {
		StringBuffer sb = new StringBuffer();
		PrismStringBuffer psb = new PrismStringBuffer(sb);

		Prism search = (Prism) session.webapp();
		Genome genome = search.genome();
		PrismConfig config = search.config();
		
		if (genome != null) {
			List<Cluster> clusters = genome.clusters();
			if (clusters.size() > 0) {
				psb.appendLine("<header class='sub'>");
				psb.appendLine("<h3>Clusters</h3>");
				psb.appendLine("<p>");
				if (genome.filename() != null) 
					psb.appendLine("Sequence: " + genome.filename() + "<br>");
				String configLink = session.webDir() + "config.txt";
				psb.appendLine("Search configuration: <a href='" + configLink + "' download>download</a><br>");
				String jsonLink = session.webDir() + genome.filename() + ".json";
				psb.appendLine("JSON Results: <a href='" + jsonLink + "' download>download</a><br>");
				psb.appendLine("These results will be available at <a href='" 
						+ session.webDir() + "'>" + session.webDir() + "</a> for 60 days.");
				psb.appendLine("</p>");
				psb.appendLine("</header>");

				// embed circular genome graph
				CircularGenomeGraph genomeGraph = genome.graph();
				if (genomeGraph != null) {
					String genomeGraphHtml = genomeGraph.html(config);
					psb.append(genomeGraphHtml);
				}

				psb.appendLine("<table class='clusterTable'>");
				psb.appendLine("<thead>");
				psb.appendLine("<tr>");
				psb.appendLine("<th style='width:162px;'>Results</th>");
				if (config.score) {
					psb.appendLine("<th style='width:552px;'>Cluster</th>");
					psb.appendLine("<th style='width:190px;'>Product</th>");				
				} else {
					psb.appendLine("<th style='width:754px;'>Cluster</th>");
				}
				psb.appendLine("</tr>");
				psb.appendLine("</thead>");
				psb.appendLine("<tbody>");

				for (Contig contig : genome.contigs()) {
					for (Cluster cluster : contig.clusters()) {
						String typeString = ClusterTypeAnalyzer.getTypeString(cluster);
						psb.appendLine("<tr>");
						psb.appendLine("<td><a href='" + session.webDir() + "cluster_" + cluster.index() + ".html'>Cluster " + 
								cluster.index() + "</a><br>" + typeString + "</td>");
						String clusterGraph = ClusterGraph.html(cluster);
						psb.appendLine("<td>" + clusterGraph + "</td>");
						if (config.score) {
							String known = ClusterAnalyzer.isKnown(cluster, config);
							psb.appendLine("<td>" + known + "</td>");
						}
						psb.appendLine("</tr>");
					}
				}

				psb.appendLine("</tbody>");
				psb.appendLine("</table>");
			} else {
				psb.appendLine("<div class='cf stageHolder firstStageHolder'>");
				psb.appendLine("<div class='stage'>Results</div>");
				psb.appendLine(
						"<p class='detail'>No clusters found in your sequence file!</p>");
				psb.appendLine("</div>");
			}
		} else {
			System.out
					.println("Error: could not write report: genome is null!");
		}

		return psb.toString();
	}

	public String writeHeader() {
		StringBuffer sb = new StringBuffer();
		PrismStringBuffer psb = new PrismStringBuffer(sb);

		String head = writeHead();
		psb.append(head);

		psb.appendLine("<body>");
		psb.appendLine("<nav id='header' class='cf'>");
		psb.appendLine("<div class='container cf'>");
		
		psb.appendLine("<h1><a href='" + context + "/'>&larr; PRISM</a></h1>");
		psb.appendLine("</div></nav>");
		
		Prism prism = (Prism) session.webapp();
		if (prism != null) {
			PrismConfig config = prism.config();
			if (config != null 
					&& !config.version.equals(new PrismConfig().version)) {
				psb.appendLine("<div class='versionWarning'>");
				psb.appendLine("<p>Warning: these results were generated with PRISM version " 
						+ config.version + ". The current version is " + new PrismConfig().version 
						+ ". Some inconsistencies may appear in your saved output due to compatibility issues!");
				psb.appendLine("</div>");
			}
		}
		return psb.toString();
	}

	public String writeClusterHeader() {
		StringBuffer sb = new StringBuffer();
		PrismStringBuffer psb = new PrismStringBuffer(sb);

		String head = writeHead();
		psb.append(head);

		String url = context + "/tasks/" + session.id() + "/";

		psb.appendLine("<body>");
		psb.appendLine("<nav id='header' class='cf'>");
		psb.appendLine("<div class='container cf'>");
		psb.appendLine("<h1><a href=" + url + ">&larr; Search #" + session.id() + "</a></h1>");
		psb.appendLine("</div></nav>");

		return psb.toString();
	}

	public String writeHead() {
		StringBuffer sb = new StringBuffer();
		PrismStringBuffer psb = new PrismStringBuffer(sb);

		psb.appendLine("<!DOCTYPE html>");
		psb.appendLine("<head>");
		psb.appendLine("<meta charset='utf-8'>");
		psb.appendLine("<meta http-equiv='X-UA-Compatible' content='IE=edge,chrome=1'>");
		psb.appendLine("<title>PRISM&mdash;Prediction Informatics for Secondary Metabolomes</title>");
        psb.appendLine("<link rel='shortcut icon' href='" + context + "/img/favicon.ico' />");
		psb.appendLine("<link rel='stylesheet' href='" + context
				+ "/css/style.css?version=1' type='text/css' media='screen' />");
		psb.appendLine("<link rel='stylesheet' href='" + context
				+ "/css/ChemDoodleWeb.css' type='text/css' media='screen' />");
		psb.appendLine("<script src='" + context + "/js/modernizr-2.6.2.min.js'></script>");
		psb.appendLine("</head>");

		return psb.toString();
	}

	public String writeFooter() {
		StringBuffer sb = new StringBuffer();
		PrismStringBuffer psb = new PrismStringBuffer(sb);

		psb.appendLine("<script src='//ajax.googleapis.com/ajax/libs/jquery/2.0.3/jquery.min.js'></script>");
		psb.appendLine("<script src='" + context + "/js/jquery-2.0.3.min.js'></script></script>");
		psb.appendLine("<script type='text/javascript' src='" + context + "/js/support.js'></script>");
		psb.appendLine("<script type='text/javascript' src='" + context + "/js/ChemDoodleWeb.js'></script>");
		psb.appendLine("</body>");
		psb.appendLine("</html>");

		return psb.toString();
	}

}
