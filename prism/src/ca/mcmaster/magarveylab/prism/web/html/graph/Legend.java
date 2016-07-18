package ca.mcmaster.magarveylab.prism.web.html.graph;

import ca.mcmaster.magarveylab.prism.data.Contig;
import ca.mcmaster.magarveylab.prism.util.Numbers;
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
		if (config.thiotemplated) {
			psb.append("<p><span class='nrps'></span> Nonribosomal peptide</p>");
			psb.append("<p><span class='pks'></span> Polyketide</p>");
			psb.append("<p class='hybrid'><span><span></span></span> Hybrid</p>");
			psb.append("<p><span class='type_ii_pks'></span> Type II polyketide</p>");
		}
		psb.append("<p><span class='tailoring'></span> Tailoring</p>");
		if (config.sugar)
			psb.append("<p><span class='sugar'></span> Sugar</p>");
		if (config.ribosomal) {
			psb.append("<p><span class='ribosomal'></span> Ribosomal</p>");
			psb.appendLine("<p><span class='propeptide'></span> Cleaved propeptide</p>");
		}
		if (config.resistance)
			psb.append("<p><span class='resistance'></span> Resistance</p>");
		if (config.regulation)
			psb.append("<p><span class='regulator'></span> Regulator</p>");
		psb.append("<p><span class='modification'></span> Other</p>");
		psb.append("</div>");

		return psb.toString();
	}

	
	/**
	 * Get the legend for this circular genome graph, as HTML.
	 * 
	 * @return graph legend as HTML
	 */
	public static String getGenomeGraphLegend(PrismConfig config) {
		StringBuffer sb = new StringBuffer();
		PrismStringBuffer psb = new PrismStringBuffer(sb);
		
		// calculate top margin
		int boxes = 1; // must have "selected cluster"
		if (config.ribosomal)
			boxes++;
		if (config.thiotemplated)
			boxes += 3;
		int top = -231 - (boxes * 11);
		
		psb.appendLine("<div class='legend' style='margin-top:" + top + "px'>");
		psb.appendLine("<p><strong>Legend:</strong></p>");
		if (config.thiotemplated) {
			psb.appendLine("<p><span class='nrps'></span> Nonribosomal peptide</p>");
			psb.appendLine("<p><span class='pks'></span> Type I polyketide</p>");
			psb.appendLine("<p><span class='type_ii_pks'></span> Type II polyketide</p>");
		}
		if (config.ribosomal)
			psb.appendLine("<p><span class='ribosomal'></span> Ribosomal peptide</p>");
		psb.appendLine("<p><span class='current'></span> Selected cluster</p>");
		psb.appendLine("</div>");
		
		return psb.toString();
	}
	
	public static String getContigGraphLegend(Contig contig, PrismConfig config) {
		StringBuffer sb = new StringBuffer();
		PrismStringBuffer psb = new PrismStringBuffer(sb);

		// calculate top margin
		int boxes = 1; // must have "selected cluster"
		if (config.ribosomal)
			boxes++;
		if (config.thiotemplated)
			boxes += 3;
		int top = -231 - (boxes * 11);

		psb.appendLine("<div class='legend' style='margin-top:" + top + "px;'>");
		psb.appendLine("<p>Contig length: " + Numbers.humanLength(contig.length()));
		psb.appendLine("<p><strong>Legend:</strong></p>");
		if (config.thiotemplated) {
			psb.appendLine("<p><span class='nrps'></span> Nonribosomal peptide</p>");
			psb.appendLine("<p><span class='pks'></span> Type I polyketide</p>");
			psb.appendLine("<p><span class='type_ii_pks'></span> Type II polyketide</p>");
		}
		if (config.ribosomal)
			psb.appendLine("<p><span class='ribosomal'></span> Ribosomal peptide</p>");
		psb.appendLine("<p><span class='current'></span> Selected cluster</p>");
		psb.appendLine("</div>");
		
		return psb.toString();
	}
	
	
}
