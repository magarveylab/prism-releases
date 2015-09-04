package ca.mcmaster.magarveylab.prism.database;

import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.prism.Prism;
import ca.mcmaster.magarveylab.prism.blast.BlastSearchResult;
import ca.mcmaster.magarveylab.prism.cluster.analysis.DomainAnalyzer;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Substrate;

public class Utils {

	public static String[] getDomainData(Domain domain) {
		int domainStart = 0;
		int domainEnd = 0;
		double domainScore = 0.0;
		String starterKS = "false";

		String domainAbbreviation = domain.type() == null ? "" 
				: DomainAnalyzer.isEpimerase(domain) ? "E" : domain.type().abbreviation();
		String domainName = domain.type() == null ? "" : domain.type().name();
		domainStart = domain.start();
		domainEnd = domain.end();
		domainScore = domain.score();

		if (DomainAnalyzer.isStarterKS(domain)) {
			starterKS = "true";
		}
		
		String[] data = { domainAbbreviation, String.valueOf(domainStart),
				String.valueOf(domainEnd), String.valueOf(domainScore),
				starterKS, domainName };
		return data;
	}

	public static ArrayList<ArrayList<String[]>> getSubstrateData(Prism prism,
			Domain domain) {
		ArrayList<ArrayList<String[]>> data = new ArrayList<ArrayList<String[]>>();
		for (Substrate substrate : domain.substrates()) {
			ArrayList<String[]> substrateData = new ArrayList<String[]>();
			String[] set = new String[4];
			set[0] = substrate.type().abbreviation();
			set[1] = String.valueOf(substrate.start());
			set[2] = String.valueOf(substrate.end());
			set[3] = String.valueOf(substrate.score());
			substrateData.add(set);
			data.add(substrateData);
		}
		return data;
	}

	public static List<List<String[]>> getSourceData(Prism prism,
			Domain domain) {
		List<List<String[]>> data = new ArrayList<List<String[]>>();
		for (BlastSearchResult result : domain.blastResults()) {
			List<String[]> substrateData = new ArrayList<String[]>();
			String[] set = new String[4];
			set[0] = result.subject();
			set[1] = String.valueOf(result.score());
			substrateData.add(set);
			data.add(substrateData);
		}
		return data;
	}

}
