package ca.mcmaster.magarveylab.prism.web.html.graph;

import java.util.List;

import ca.mcmaster.magarveylab.enums.DomainFamilies;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.prism.cluster.analysis.DomainAnalyzer;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.data.Substrate;
import ca.mcmaster.magarveylab.prism.util.PrismStringBuffer;
import ca.mcmaster.magarveylab.prism.util.Strings;

public class OrfAssemblyGraph {

	public static String html(Orf orf) {
		StringBuffer sb = new StringBuffer();
		PrismStringBuffer psb = new PrismStringBuffer(sb);

		List<Domain> domains = orf.domains();
		for (Domain domain : domains) {
			StringBuffer typeName = new StringBuffer();
			StringBuffer abbreviation = new StringBuffer();

			DomainType type = domain.type();
			if (type == ThiotemplatedDomains.CONDENSATION) {
				if (DomainAnalyzer.isEpimerase(domain)) {
					abbreviation.append("E");
					typeName.append("epimerization");
				} else if (DomainAnalyzer.isCStarter(domain)) {
					abbreviation.append("C*");
					typeName.append("condensation starter");
				} else if (DomainAnalyzer.isCyclization(domain)) {
					abbreviation.append("Cy");
					typeName.append("condensation cyclization");
				} else {
					abbreviation.append(type.abbreviation());
					typeName.append("condensation");
				}
			} else if (DomainAnalyzer.isSubstrateDomainType(type)) {
				typeName.append(type.toString().toLowerCase());
				Substrate top = domain.topSubstrate();
				abbreviation.append(top.type().abbreviation());
			} else if (type == ThiotemplatedDomains.KETOSYNTHASE) {
				if (DomainAnalyzer.isStarterKS(domain)) {
					abbreviation.append("KS<sub>Q</sub>");
				} else {
					abbreviation.append(type.abbreviation());
				}
				typeName.append("ketosynthase");
			} else if (domain.family() == DomainFamilies.TYPE_II_POLYKETIDE) { 
				abbreviation.append(type.abbreviation());
				typeName.append(type.toString());
				typeName.append(" type_ii_pks");
			} else if (domain.family() == DomainFamilies.TAILORING) { 
				typeName.append("tailoring");
				abbreviation.append(type.abbreviation());
			} else {
				typeName.append(type.toString().toLowerCase() + " " + type.family().toString().toLowerCase());
				abbreviation.append(type.abbreviation());
			}

			psb.append("<span style='font-size:" + Strings.calculateFontSize(abbreviation.toString()) + "px;' "
					+ "class='orfAssemblyGraphDomain " + typeName.toString() + "'>" + abbreviation.toString() + 
					"</span>");
		}


		return psb.toString();
	}

}
