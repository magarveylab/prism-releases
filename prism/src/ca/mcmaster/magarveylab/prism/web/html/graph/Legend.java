package ca.mcmaster.magarveylab.prism.web.html.graph;

import ca.mcmaster.magarveylab.prism.util.PrismStringBuffer;
import ca.mcmaster.magarveylab.prism.web.PrismConfig;

/**
 * Create legends for graphical output.
 * 
 * @author skinnider
 *
 */
public class Legend {

	/**
	 * Get the legend for an open reading frame graph, as HTML.
	 * 
	 * @param config
	 *            the current PRISM configuration
	 * @return the orf graph legend, as HTML
	 */
	public static String getOrfGraphLegend(PrismConfig config) {
		StringBuffer sb = new StringBuffer();
		PrismStringBuffer psb = new PrismStringBuffer(sb);

		psb.append("<div class='legend cf'>");
		psb.append("<p><span class='nrps'></span> Nonribosomal peptide</p>");
		psb.append("<p><span class='pks'></span> Polyketide</p>");
		psb.append("<p class='hybrid'><span><span></span></span> Hybrid</p>");
		psb.append("<p><span class='type_ii_pks'></span> Type II polyketide</p>");
		psb.append("<p><span class='sugar'></span> Sugar</p>");
		psb.append("<p><span class='tailoring'></span> Tailoring</p>");
		psb.append("<p><span class='beta_lactam'></span> Beta-lactam</p>");
		if (config.ribosomal)
			psb.append("<p><span class='ribosomal'></span> Ribosomal</p>");
		if (config.resistance)
			psb.append("<p><span class='resistance'></span> Resistance</p>");
		if (config.regulation)
			psb.append("<p><span class='regulator'></span> Regulator</p>");
		if (config.aminoglycoside)
			psb.append("<p><span class='aminoglycoside'></span> Aminoglycoside</p>");
		if (config.nucleoside)
			psb.append("<p><span class='nucleoside'></span> Nucleoside</p>");
		psb.append("<p><span class='modification'></span> Other</p>");
		psb.append("</div>");

		return psb.toString();
	}

	
	/**
	 * Get the legend for this circular genome graph, as HTML.
	 * 
	 * @return graph legend as HTML
	 */
	public static String getGenomeGraphLegend() {
		StringBuffer sb = new StringBuffer();
		PrismStringBuffer psb = new PrismStringBuffer(sb);
		
		psb.appendLine("<div class='legend'>");
		psb.appendLine("<p><strong>Legend:</strong></p>");
		psb.appendLine("<p><span class='nrps'></span> Nonribosomal peptide</p>");
		psb.appendLine("<p><span class='pks'></span> Type I polyketide</p>");
		psb.appendLine("<p><span class='type_ii_pks'></span> Type II polyketide</p>");
		psb.appendLine("<p><span class='beta_lactam'></span> Beta-lactam</p>");
		psb.appendLine("<p><span class='ribosomal'></span> Ribosomal peptide</p>");
		psb.appendLine("<p><span class='aminoglycoside'></span> Aminoglycoside</p>");
		psb.appendLine("<p><span class='nucleoside'></span> Nucleoside</p>");
		psb.appendLine("<p><span class='current'></span> Selected cluster</p>");
		psb.appendLine("</div>");
		
		return psb.toString();
	}
	
}
