package ca.mcmaster.magarveylab.prism.web.html.graph;

import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.util.PrismStringBuffer;

public class ClusterAssemblyGraph {

	public static String html(Cluster cluster) {
		StringBuffer sb = new StringBuffer();
		PrismStringBuffer psb = new PrismStringBuffer(sb);

		psb.appendLine("<div class='clusterAssemblyGraph cf'>");
		for (Orf orf : cluster.orfs()) 
			if (orf.domains().size() > 0){
				psb.appendLine("<div class='orfAssemblyGraph'>");
				psb.appendLine("<span class='orfAssemblyGraphTitle'>" + orf.name() + "</span>");
				String assembly = OrfAssemblyGraph.html(orf);
				psb.append(assembly);
				psb.appendLine("</div>");
			}
		psb.appendLine("</div>");

		return psb.toString();
	}

}
