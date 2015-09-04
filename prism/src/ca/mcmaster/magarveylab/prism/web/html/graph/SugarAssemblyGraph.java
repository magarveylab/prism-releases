package ca.mcmaster.magarveylab.prism.web.html.graph;

import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.enums.DomainFamilies;
import ca.mcmaster.magarveylab.enums.domains.DeoxySugarDomains;
import ca.mcmaster.magarveylab.prism.cluster.analysis.SugarAnalyzer;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.sugar.DeoxySugar;
import ca.mcmaster.magarveylab.prism.data.sugar.Sugar;
import ca.mcmaster.magarveylab.prism.util.PrismStringBuffer;
import ca.mcmaster.magarveylab.prism.util.Strings;

public class SugarAssemblyGraph {
	
	public static String html(List<Sugar> sugars, Cluster cluster) {
		StringBuffer sb = new StringBuffer();
		PrismStringBuffer psb = new PrismStringBuffer(sb);

		for (Sugar sugar : sugars) {
			psb.append("<span class='sugarAssemblyGraph'>");
			psb.append("<span style='font-size:" + Strings.calculateFontSize("GTr") + "px;' "
					+ "class='orfAssemblyGraphDomain glycosyltransferase'>GTr</span>");
			
			if (sugar.isDeoxygenated()) {
				DeoxySugar d = (DeoxySugar) sugar;
				for (DeoxySugarDomains gene : d.genes()) {
					String abbreviation = gene.abbreviation();
					String type = gene.toString().toLowerCase();
					String missing = (SugarAnalyzer.contains(gene, cluster)) ? "" : "missing";
					psb.append("<span style='font-size:" + Strings.calculateFontSize(abbreviation.toString()) + "px;' "
							+ "class='orfAssemblyGraphDomain " + type + " " + missing + "'>" 
							+ abbreviation.toString() + "</span>");
				}
			}
			
			psb.append("<span class='sugarName'>" + sugar.toString() + "</span>");
			psb.append("</span>");
		}
		
		// get unused domains HTML
		StringBuffer remaining = new StringBuffer();
		List<DeoxySugarDomains> functions = new ArrayList<DeoxySugarDomains>();
		for (Sugar sugar : sugars)
			if (sugar.isDeoxygenated()) {
				DeoxySugar d = (DeoxySugar) sugar;
				for (DeoxySugarDomains function : d.genes())
					functions.add(function);
			}
		for (Domain domain : cluster.domains(DomainFamilies.SUGAR)) {
			if (!SugarAnalyzer.contains(domain, functions)) {
				StringBuffer typeName = new StringBuffer();
				StringBuffer abbreviation = new StringBuffer();

				typeName.append(domain.type().toString().toLowerCase());
				abbreviation.append(domain.type().abbreviation());

				remaining.append("<span style='font-size:" + Strings.calculateFontSize(abbreviation.toString()) + "px;' "
						+ "class='orfAssemblyGraphDomain " + typeName.toString() + "'>" + abbreviation.toString() + "</span>");
			}
		}
		if (remaining.length() > 0) {
			psb.append("<span class='sugarAssemblyGraph remaining'>");
			psb.append(remaining.toString());
			psb.append("<span class='sugarName'>Unused</span>");
			psb.append("</span>");
		}

		return psb.toString();
	}

}
