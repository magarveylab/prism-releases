package ca.mcmaster.magarveylab.prism.web.html;

import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.List;

import org.openscience.cdk.exception.CDKException;

import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.prism.Prism;
import ca.mcmaster.magarveylab.prism.blast.BlastSearchResult;
import ca.mcmaster.magarveylab.prism.cluster.analysis.DomainAnalyzer;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.data.Propeptide;
import ca.mcmaster.magarveylab.prism.data.Substrate;
import ca.mcmaster.magarveylab.prism.util.PrismStringBuffer;
import ca.mcmaster.magarveylab.prism.web.PrismConfig;
import ca.mcmaster.magarveylab.prism.web.html.graph.OrfAssemblyGraph;
import ca.mcmaster.magarveylab.prism.web.html.graph.OrfGraph;
import ca.mcmaster.magarveylab.prism.web.html.graph.SequenceGraph;
import ca.mcmaster.magarveylab.wasp.session.Session;
import ca.mcmaster.magarveylab.prism.motif.Motif;
import ca.mcmaster.magarveylab.prism.motif.MotifList;
import ca.mcmaster.magarveylab.prism.orfs.GenePredictionModes;

/**
 * Generates HTML output for a single orf. 
 * @author skinnider
 *
 */
public class OrfReport {

	public static String getHTML(Orf orf, Session session) throws IOException, CDKException {
		StringBuffer sb = new StringBuffer();
		PrismStringBuffer psb = new PrismStringBuffer(sb);
		Prism prism = (Prism) session.webapp();
		PrismConfig config = prism.config();

		psb.appendLine("<article class='orf cf'>");

		// append orf header
		psb.appendLine("<h1>" + orf.name() + "</h1>");

		// append domains table header
		//	psb.appendLine("<h2>Domains:</h2>");
		psb.appendLine("<table>");
		psb.appendLine("<thead>");
		psb.appendLine("<tr>");
		psb.appendLine("<th style='width:217px;'>Domain</th>");
		psb.appendLine("<th style='width:63px;'>Start</th>");
		psb.appendLine("<th style='width:63px;'>End</th>");
		psb.appendLine("<th style='width:70px;'>Score</th>");
		psb.appendLine("</tr>");
		psb.appendLine("</thead>");
		psb.appendLine("<tbody>");

		// sort domains for tabular output
		List<Domain> domains = orf.domains();
		for (Domain domain : domains) {
			String name = null;

			if (DomainAnalyzer.isEpimerase(domain)) {
				name = "Epimerization";
			} else {
				name = domain.type().fullName();
			}

			String details = (domain.blastResults().size() > 0 
					|| domain.substrates().size() > 0 
					|| orf.hasPropeptides()) 
					? "<a href onclick='return showDetails(this, event)'>+</a>": "";

			String score = domain.score() < 0 ? "N/A" : domain.score() + "";
			
			psb.appendLine("<tr>");
			psb.appendLine("<td>" + name + details + "</td>");
			psb.appendLine("<td>" + domain.start() + "</td>");
			psb.appendLine("<td>" + domain.end() + "</td>");
			psb.appendLine("<td>" + score + "</td>");
			psb.appendLine("</tr>");
			psb.appendLine("<tbody class='details'>");

			if (domain.blastResults().size() > 0) {
				boolean gtr = false;
				if (domain.type() == TailoringDomains.GLYCOSYLTRANSFERASE)
					gtr = true;

				// BLASTp results
				for (int i = 0; i < config.display && domain.blastResults().size() > i; i++) {
					BlastSearchResult result = domain.blastResults().get(i);
					psb.appendLine("<tr class='substrate'>");
					String subject = (gtr) ? result.subject().substring(0, result.subject().lastIndexOf("_")) 
							: result.subject();
					psb.appendLine("<td>" + subject + "</td>");
					psb.appendLine("<td></td>");
					psb.appendLine("<td></td>");
					psb.appendLine("<td>" + result.score() + "</td>");
					psb.appendLine("</tr>");
				}
			} else if (orf.hasPropeptides()) {
				for (Propeptide propeptide : orf.propeptides()) {
					psb.appendLine("<tr class='substrate'>");
					psb.appendLine("<td>" +"<strong>" + "Cleaved propeptide" + "</strong>" + "</td>");
					psb.appendLine("<td>" + propeptide.getStart() + "</td>");
					psb.appendLine("<td>" + propeptide.getEnd() + "</td>");
					psb.appendLine("<td>" + " " + "</td>");
					psb.appendLine("</tr>");
					psb.appendLine("<tr class='substrate'>");
					psb.appendLine("<td colspan=4  style='word-wrap:break-word;'><strong>Sequence:</strong> "
							+ propeptide.getSequence() + "</td>");
					psb.appendLine("</tr>");
					
					MotifList motifs = propeptide.getMotifs();
					for (Motif motif : motifs) {
						String pValue;
						double p = motif.getPValue();
						if (p < 0.005) {
							String pString = p + "";
							String[] split = pString.split("E");
							String num = split[0].length() >= 4 ? split[0].substring(0,4) : split[0];
							pValue = num + "E" + split[1];
						} else {
							NumberFormat formatter = new DecimalFormat("#0.00");
							pValue = formatter.format(motif.getPValue());
						}
						
						psb.appendLine("<tr class='substrate'>");
						psb.appendLine("<td>" + motif.getType().getFullName() + "</td>");
						psb.appendLine("<td>" + motif.getStart() + "</td>");
						psb.appendLine("<td>" + motif.getEnd() + "</td>");
						psb.appendLine("<td>" + pValue + "</td>");
						psb.appendLine("</tr>");
					}
				}
		} else if (domain.substrates().size() > 0) {
				// substrates
				List<Substrate> substrates = domain.substrates();
				for (int i = 0; substrates.size() > i && i < config.display; i++) {
					Substrate substrate = substrates.get(i);
					psb.appendLine("<tr class='substrate'>");
					
					if(substrate.type() == null){
						psb.appendLine("<td> null </td>");
					}
					else{
						psb.appendLine("<td>" + substrate.type().fullName() + "</td>");
					}
					
					psb.appendLine("<td>" + substrate.start() + "</td>");
					psb.appendLine("<td>" + substrate.end() + "</td>");
					psb.appendLine("<td>" + substrate.score() + "</td>");
					psb.appendLine("</tr>");
				} 
			} 
			
			psb.appendLine("</tbody>");

		}
		psb.appendLine("</tbody>");
		psb.appendLine("</table>");			

		// get detection mdoe
		String detection = null;
		if (orf.getMode() == GenePredictionModes.PRODIGAL) {
			detection = "predicted by Prodigal";
		} else {
			detection = "potential coding sequence";
		}
		
		psb.appendLine("<p>");
		if (detection != null)
			psb.appendLine("<strong>Detection: </strong>" + detection + "<br>");
		psb.appendLine("<strong>Start: </strong>" + orf.start() + " | ");
		psb.appendLine("<strong>End: </strong>" + orf.end() + " | ");
		psb.appendLine("<strong>Frame: </strong>" + orf.frame());
		psb.appendLine("</p>");

		// append sequence
		psb.appendLine(SequenceGraph.getHTML(orf));

		psb.appendLine("<p><strong>Domain analysis:</strong><br>");
		String graph = OrfGraph.html(orf);
		psb.append(graph);
		psb.appendLine("<p>");

		psb.appendLine("<p class='cf'><strong>Biosynthetic assembly:</strong><br></p>");
		psb.appendLine("<div class='orfAssemblyGraph'>");
		String assembly = OrfAssemblyGraph.html(orf);
		psb.append(assembly);
		psb.appendLine("</div>");

		psb.appendLine("</article>");
		return psb.toString();
	}

}
