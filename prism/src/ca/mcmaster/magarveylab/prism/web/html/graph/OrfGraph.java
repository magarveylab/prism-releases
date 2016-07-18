package ca.mcmaster.magarveylab.prism.web.html.graph;

import ca.mcmaster.magarveylab.enums.DomainFamilies;
import ca.mcmaster.magarveylab.prism.cluster.analysis.DomainAnalyzer;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.data.Propeptide;
import ca.mcmaster.magarveylab.prism.util.PrismStringBuffer;
import ca.mcmaster.magarveylab.prism.util.Sorter;

public class OrfGraph {
	
	/**
	 * Get the HTML graph of a single orf.
	 * 
	 * @param orf
	 *            the orf to graph
	 * @param config
	 *            current Prism configuration
	 * @return graph in HTML format
	 */
	public static String html(Orf orf) {
		StringBuffer sb = new StringBuffer();
		PrismStringBuffer psb = new PrismStringBuffer(sb);

		int length = orf.length() / 3;
		int divWidth = 440;
		double unit = (double) divWidth / length;

		psb.appendLine("<div class='orfGraph cf'>");

		for (Domain domain : orf.domains()) {
			String type = null;

			// assign type
			if (DomainAnalyzer.isEpimerase(domain)) {
				type = "epimerization";
			} else if (domain.family() == DomainFamilies.TAILORING) { 
				type = "tailoring";
			} else if (domain.family() == DomainFamilies.TYPE_II_POLYKETIDE) { 
				type = "type_ii_pks";
			} else {
				type = domain.type().toString().toLowerCase() + " " + domain.family().toString().toLowerCase();	
			}
			
			double left = domain.start() * unit;
			double width = (domain.end() - domain.start() - 1) * unit;
			width = (width > 0) ? width : 1.0;
			width = (left + width > divWidth) ? divWidth - left : width;
			psb.appendLine("<span class='orfGraphDomain " + type + "' style='" +
					"width: " + (int) width + "px; " + "left: " + (int) left + "px; " +
					"top: " + 0 + "px; " + "'></span>");

			// sort propeptides
			Sorter.sortPropeptides(orf.propeptides());
			
			// create propeptides 
			for (Propeptide propeptide : orf.propeptides()) {
				String propeptideType = "propeptide";
				double propeptideLeft = propeptide.getStart() * unit;
				double propeptideWidth = propeptide.getLength() * unit; 
				psb.appendLine("<span class='orfGraphDomain " + propeptideType + "' style='" +
						"width: " + (int) propeptideWidth + "px; " + "left: " + (int) propeptideLeft + "px; " +
						"top: " + 0 + "px; " + "'></span>");

			}
		}

		psb.appendLine("</div>");

		return psb.toString();
	}

}
