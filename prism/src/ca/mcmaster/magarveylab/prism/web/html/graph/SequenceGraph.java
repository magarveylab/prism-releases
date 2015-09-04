package ca.mcmaster.magarveylab.prism.web.html.graph;

import java.util.List;

import ca.mcmaster.magarveylab.prism.cluster.analysis.DomainAnalyzer;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.util.PrismStringBuffer;

public class SequenceGraph {
	
	public static String getHTML(Orf orf) {
		StringBuffer sb = new StringBuffer();
		PrismStringBuffer psb = new PrismStringBuffer(sb);
		
		if (orf.sequence().isEmpty() || orf.sequence().equals(""))
			return "";
		
		psb.appendLine("<p><strong>Sequence:</strong> " 
				+ "<a class='showLink' href='#' onClick='return toggleSequence(this, event)'>show</a> "
				+ "<span class='sequence'>");
		
		String sequence = orf.sequence();
		List<Domain> domains = orf.domains();
		if (domains.size() > 0) {
			Domain first = domains.get(0);
			int firstStart = first.start();
			if (firstStart > 0)
				psb.append(sequence.substring(0, firstStart));
			for (int i = 0; i < domains.size() - 1; i++) {
				Domain domain = domains.get(i);
				Domain next = domains.get(i+1);
				
				String type = null;
				type = DomainAnalyzer.getDomainTypeNameForCSS(domain).toLowerCase();
				
				int start = domain.start();
				int end = domain.end();
				String subsequence = sequence.substring(start, end);
 				psb.append("<span class='" + type + " " + domain.family().toString().toLowerCase() 
 						+ "'>" + subsequence + "</span>");

				int nextStart = next.start();
				
				if (nextStart > end && nextStart > 0) 
					psb.append(sequence.substring(end, nextStart));
			}
			
			Domain last = domains.get(domains.size() - 1);
			int start = last.start();
			int end = last.end();
			String subsequence = sequence.substring(start, end); 
			String type = DomainAnalyzer.getDomainTypeNameForCSS(last).toLowerCase();
			psb.append("<span class='" + type + " " + last.family().toString().toLowerCase()
					+ "'>" + subsequence + "</span>");

			int lastEnd = last.end();
			if (lastEnd < sequence.length())
				psb.append(sequence.substring(lastEnd, sequence.length()));
		}
		
		psb.appendLine("</span> <a class='hideLink' href='#' onClick='return toggleSequence(this, event)'>(hide)</a></p>");
		return psb.toString();
	}

}
