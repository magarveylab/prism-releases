package ca.mcmaster.magarveylab.prism.web.html.graph;

import java.util.Iterator;
import java.util.List;

import ca.mcmaster.magarveylab.enums.OrfTypes;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.util.PrismStringBuffer;

/**
 * A graph of all the biosynthetic open reading frames in a putative natural product gene cluster.
 * @author skinnider
 *
 */
public class ClusterGraph {
	
	/**
	 * Get the HTML graph of this biosynthetic gene cluster.
	 * @param cluster	cluster in question
	 * @return			graph in HTML format
	 */
	public static String html(Cluster cluster) {
		StringBuffer sb = new StringBuffer();
		PrismStringBuffer psb = new PrismStringBuffer(sb);
		psb.appendLine("<div class='clusterGraph'>");
		
		List<Orf> orfs = cluster.orfs();
		Iterator<Orf> iterator = orfs.iterator();
		while (iterator.hasNext()) {
			Orf orf = iterator.next();
			StringBuffer type = new StringBuffer();
			boolean scaffold = false;
			if (orf.type() != null) {
				if (orf.type() == OrfTypes.INACTIVE || orf.type() == OrfTypes.NULL)
					continue;
				type.append(orf.type().toString());
				scaffold = (orf.type() == OrfTypes.HYBRID || orf.type() == OrfTypes.NRPS || 
						orf.type() == OrfTypes.PKS) ? true : false;
			} else {
				continue;
			}
			
			double width = (double) orf.length() * 0.02;
			if (width > 260) 
				width = 260;
			if (width < 10)
				width = 10;
			if (width <= 50 && scaffold)
				width = 50;
			
			psb.append("<span class='clusterGraphOrf ");
			if (!orf.frame().startsWith("+"))
				psb.append("negative ");
			psb.append(type.toString() + "' style='width:" + (int) width + "px;'>");
			if (orf.type() != null && orf.type() == OrfTypes.HYBRID)
				psb.append("<span class='hybrid'></span>");
			if (scaffold) {
				psb.append("<span class='name'>" + orf.name() + "</span>");
			} else {
				psb.append("<span class='name'>&nbsp;</span>");
			}
			psb.append("</span>");
		}
		
		psb.appendLine("</div>");
		return psb.toString();
	}
	
}
