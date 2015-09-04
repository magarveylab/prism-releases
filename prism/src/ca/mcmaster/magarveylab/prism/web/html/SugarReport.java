package ca.mcmaster.magarveylab.prism.web.html;

import java.util.List;

import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.sugar.Sugar;
import ca.mcmaster.magarveylab.prism.util.PrismStringBuffer;
import ca.mcmaster.magarveylab.prism.util.Sorter;
import ca.mcmaster.magarveylab.prism.web.html.graph.SugarAssemblyGraph;

public class SugarReport {

	public static String getHTML(Cluster cluster) {		
		StringBuffer sb = new StringBuffer();
		PrismStringBuffer psb = new PrismStringBuffer(sb);
	
		List<List<Sugar>> combinations = cluster.sugars();
		Sorter.sortSugarCombinationsByName(combinations);
		System.out.println("[SugarReport] Graphing " + combinations.size() + " sugar combinations");

		psb.appendLine("<div class='sugars'>");
		psb.appendLine("<h3>Sugars:</h3>");
		String s = combinations.size() > 1 ? "s" : "";
		psb.appendLine("<p>" + combinations.size() + " possible sugar combination" + s + " identified.</p>");

		psb.appendLine("<table>");
		psb.appendLine("<thead>");
		psb.appendLine("<tr>");
		psb.appendLine("<th style='width:331px;'>Sugars</th>");
		psb.appendLine("<th style='width:;'>Domains</th>");
		psb.appendLine("</tr>");
		psb.appendLine("</thead>");
		psb.appendLine("<tbody>");

		if (combinations.size() > 5) {
			for (int i = 0; i < 5; i++) {
				List<Sugar> combination = combinations.get(i);
				psb.appendLine(getRow(combination, cluster));
			}
			psb.appendLine("</tbody>");
			psb.appendLine("<tbody id='sugarTableHidden'>");
			for (int i = 5; i < combinations.size(); i++) {
				List<Sugar> combination = combinations.get(i);
				psb.appendLine(getRow(combination, cluster));
			}
			psb.appendLine("</tbody>");
			psb.appendLine("<tbody id='sugarTableShow'>");
			psb.appendLine("<tr><td colspan='2'><span onclick='return toggleSugars(this);'><span id='showText'>Show all "
					+ combinations.size() + "</span><span id='hideText'>Collapse</span></span></td></tr>");
			psb.appendLine("</tbody>");
		} else {
			for (List<Sugar> combination : combinations)
				if (combination != null)
					psb.appendLine(getRow(combination, cluster));
		}

		psb.appendLine("</tbody>");
		psb.appendLine("</table>");
		psb.appendLine("</div>");
		
		return psb.toString();
	}
	
	public static String getRow(List<Sugar> combination, Cluster cluster) {
		StringBuffer sb = new StringBuffer();
		PrismStringBuffer psb = new PrismStringBuffer(sb);
		
		if (combination == null)
			return psb.toString();

		psb.appendLine("<tr>");
		psb.appendLine("<td>");
		if (combination.size () > 0) {
			for (int i = 0; i < combination.size() - 1; i++)
				psb.appendLine(combination.get(i).toString() + ", ");
			psb.appendLine(combination.get(combination.size() - 1).toString());
		} 
		psb.appendLine("</td>");
		psb.appendLine("<td>");
		psb.appendLine(SugarAssemblyGraph.html(combination, cluster));
		psb.appendLine("</td>");
		psb.appendLine("</tr>");				
		return psb.toString();
	}

}
