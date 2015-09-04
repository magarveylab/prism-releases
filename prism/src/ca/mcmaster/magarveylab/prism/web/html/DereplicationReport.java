package ca.mcmaster.magarveylab.prism.web.html;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.prism.Prism;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Library;
import ca.mcmaster.magarveylab.prism.homology.HomologousClusterSearch;
import ca.mcmaster.magarveylab.prism.homology.data.HomologousCluster;
import ca.mcmaster.magarveylab.prism.tanimoto.data.TanimotoScore;
import ca.mcmaster.magarveylab.prism.util.PrismStringBuffer;
import ca.mcmaster.magarveylab.prism.util.Strings;
import ca.mcmaster.magarveylab.prism.web.PrismConfig;
import ca.mcmaster.magarveylab.prism.web.html.graph.ClusterHomologGraph;
import ca.mcmaster.magarveylab.wasp.session.Session;

public class DereplicationReport {

	public static String getHTML(Cluster cluster, Session session) throws IOException {		
		StringBuffer sb = new StringBuffer();
		PrismStringBuffer psb = new PrismStringBuffer(sb);
		
		Prism prism = (Prism) session.webapp();
		PrismConfig config = prism.config();
	
		List<HomologousCluster> homologs = cluster.homologs();
		
		psb.appendLine("<div class='clusterHomologs'>");
		psb.appendLine("<h3>Homologous clusters:</h3><p>");

		psb.appendLine("<table class='homology'>");
		psb.appendLine("<thead>");
		psb.appendLine("<tr>");
		psb.appendLine("<th style='width:151px;'>Name</th>");
		psb.appendLine("<th style='width:432px;'>Cluster</th>");
		psb.appendLine("<th style='width:76px;'>Domain score</th>");
		psb.appendLine("<th style='width:76px;'>Coverage</th>");
		psb.appendLine("<th style='width:76px;'>Average domain score</th>");
		psb.appendLine("<th style='width:80px;'>Identity score</th>");
		psb.appendLine("</tr>");
		psb.appendLine("</thead>");
		psb.appendLine("<tbody>");
		for (int i = 0; i < config.display && i < homologs.size(); i++) {
			HomologousCluster homolog = homologs.get(i);
		//	if (homolog.domainScore() > 0) {
				psb.append("<tr>");
				psb.append("<td>" + Strings.removeUnderscores(Strings.capitalizeFirstLetter(homolog.name())) + "</td>");
				psb.append("<td>" + ClusterHomologGraph.html(homolog, session) + "</td>");
				double domainScore = Double.parseDouble(String.format("%.2f", homolog.domainScore()));
				psb.append("<td>" + domainScore + "</td>");
				double coverage = Double.parseDouble(String.format("%.2f", 
						1.0d * homolog.coverage() / HomologousClusterSearch.getDomainSize(cluster)));
				psb.append("<td>" + coverage + "</td>");
				double weightedDomainScore = Double.parseDouble(String.format("%.2f", 
						1.0d * homolog.domainScore() / homolog.coverage()));
				psb.append("<td>" + weightedDomainScore + "</td>");
				double identityScore = Double.parseDouble(String.format("%.2f", homolog.identityScore()));
				psb.append("<td>" + identityScore + "</td>");
				psb.append("</tr>");
		//	}
		}
		psb.appendLine("</tbody>");
		psb.appendLine("</table>");
		
		Library library = cluster.library();
		if (library != null && library.scaffolds().size() > 0) {
			List<TanimotoScore> scores = cluster.scores();
			if (scores.size() > 0) {
				psb.appendLine("<h3 class='tanimotoHeader'>Similar compounds:</h3><p>");
				
				// avoid duplicates
				List<String> used = new ArrayList<String>();
				
				psb.appendLine("<table class='tanimoto'>");
				psb.appendLine("<thead>");
				psb.appendLine("<tr>");
				psb.appendLine("<th style='width:287px;'>Name(s)</th>");
				psb.appendLine("<th style='width:100px;'>Scaffold</th>");
				psb.appendLine("<th style='width:60px;'>ECFP6</th>");
				psb.appendLine("<th style='width:60px;'>FCFP6</th>");
				psb.appendLine("<th style='width:400px;'>SMILES</th>");
				psb.appendLine("</tr>");
				psb.appendLine("</thead>");
				psb.appendLine("<tbody>");
				int size = config.display;
				for (int i = 0; i < size && i < scores.size() - 1; i++) {
					TanimotoScore ts = scores.get(i);
					
					String name = ts.target().name();
					String smiles = ts.target().smiles();
					String query = ts.query().name();
					float ecfp6 = ts.score("ecfp6");
					float fcfp6 = ts.score("fcfp6");
					
					boolean isUsed = false;
					String[] split = name.split(",");
					for (String s : split)
						if (used.indexOf(s) != -1)
							isUsed = true;
					for (String s : split)
						used.add(s.trim());
					
					if (!isUsed) {
						psb.appendLine("<tr>");
						psb.appendLine("<td>" + name + "</td>");
						psb.appendLine("<td>" + query + "</td>");
						psb.appendLine("<td>" + String.format("%.2f", ecfp6) + "</td>");
						psb.appendLine("<td>" + String.format("%.2f", fcfp6) + "</td>");
						psb.appendLine("<td class='smiles'>" + smiles + "</td>");
						psb.appendLine("</tr>");
					} else {
						if (size < scores.size() - 1)
							size++;
					}
				}
				psb.appendLine("</tbody>");
				psb.appendLine("</table>");
			} else {
				System.out.println("[ClusterReport] ERROR: No tanimoto scores!");
			}
		}
		
		psb.appendLine("</div>");
		return psb.toString();
	}	
	
}
