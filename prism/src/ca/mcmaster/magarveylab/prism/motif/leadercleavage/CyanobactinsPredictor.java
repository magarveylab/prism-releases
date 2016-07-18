package ca.mcmaster.magarveylab.prism.motif.leadercleavage;

import java.util.List;
import java.util.ArrayList;

import ca.mcmaster.magarveylab.enums.RibosomalPrecursorMotifs;
import ca.mcmaster.magarveylab.enums.clusters.RibosomalClusterTypes;
import ca.mcmaster.magarveylab.prism.cluster.analysis.OrfAnalyzer;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.data.Propeptide;
import ca.mcmaster.magarveylab.prism.motif.*;

/**
 * Predicts cleavage sites for cyanobactins.
 * 
 * @author skinnider
 */
public class CyanobactinsPredictor extends AbstractLeaderPredictor implements
		LeaderPredictor {

	public List<Propeptide> cleaveLeader(Domain domain, Cluster cluster) {
		List<Propeptide> propeptides = new ArrayList<Propeptide>();
		String leader = OrfAnalyzer.getParentOrf(domain, cluster).sequence();
		MotifList motifs = domain.motifs();
		
		Motif clade1_n1 = motifs.getBestMotif(RibosomalPrecursorMotifs.Cyanobactin_clade_1_N1);
		Motif clade1_n2 = motifs.getBestMotif(RibosomalPrecursorMotifs.Cyanobactin_clade_1_N2);
		Motif clade1_c = motifs.getBestMotif(RibosomalPrecursorMotifs.Cyanobactin_clade_1_C);
		Motif clade2_n1 = motifs.getBestMotif(RibosomalPrecursorMotifs.Cyanobactin_clade_2_N1);
		Motif clade2_n2 = motifs.getBestMotif(RibosomalPrecursorMotifs.Cyanobactin_clade_2_N2);
		Motif clade2_r = motifs.getFirstMotif(RibosomalPrecursorMotifs.Cyanobactin_clade_2_R);
		Motif clade2_r2 = motifs.getFirstMotif(RibosomalPrecursorMotifs.Cyanobactin_clade_2_R2);
		// only care about the first repeating motif 
		
		boolean clade1 = (clade1_n1 != null && clade1_c != null)
				&& (clade2_r == null || clade1_c.getPValue() < clade2_r.getPValue());
		boolean clade2 = (clade2_n1 != null && clade2_n2 != null && clade2_r != null)
				&& (clade1_c == null || clade2_r.getPValue() < clade1_c.getPValue());
		boolean clade2r2 = (clade2_n1 != null && clade2_n2 != null && clade2_r2 != null);

		int start = -1, end = -1;
		int length = leader.length();
		MotifList precursorMotifs = new MotifList();
		if (clade1) {
			System.out.println("[CyanobactinsPredictor] "
					+ "Cleaving clade 1 cyanobactin");
			precursorMotifs.addMotif(clade1_n1);
			precursorMotifs.addMotif(clade1_c);
			end = (length > clade1_c.getEnd() - 7 && clade1_c.getEnd() - 7 >= 0) ? 
					clade1_c.getEnd() - 7 : 0;
			if (clade1_n2 != null) {
				precursorMotifs.addMotif(clade1_n2);
				start = clade1_n2.getEnd();
			} else {
				start = (length > clade1_n1.getEnd() + 12) ? 
						clade1_n1.getEnd() + 12 : length;
			}
		} else if (clade2) {
			precursorMotifs.addMotif(clade2_n1);
			precursorMotifs.addMotif(clade2_n2);
			precursorMotifs.addMotif(clade2_r);
			System.out.println("[CyanobactinsPredictor] "
					+ "Cleaving clade 2 cyanobactin");
			start = (length > clade2_n2.getEnd() + 4) ? 
					clade2_n2.getEnd() + 4 : 0;
			end = (clade2_r.getEnd() - 4 >= 0) ? clade2_r.getEnd() - 4 : 0;
		} else if (clade2r2) {
			precursorMotifs.addMotif(clade2_n1);
			precursorMotifs.addMotif(clade2_n2);
			precursorMotifs.addMotif(clade2_r2);
			System.out.println("[CyanobactinsPredictor] "
					+ "Cleaving clade 2 cyanobactin with repeating motif 2");
			start = (length > clade2_n2.getEnd() + 4) ? 
					clade2_n2.getEnd() + 4 : 0;
			end = clade2_r2.getStart() + 3;
		} else {
			System.out.println("[CyanobactinsPredictor] "
					+ "Couldn't identify cyanobactin precursor as clade 1 or "
					+ "clade 2. Motifs: ");
			for (Motif motif : motifs)
				System.out.println("\t" + motif.toString());
			return propeptides;
		}
		
		// max size 15
		if (length > start && length > end && end - start <= 15) {
			String sequence = leader.substring(start, end);
			Propeptide propeptide = new Propeptide(sequence, start, end);
			propeptide.setMotifs(precursorMotifs);
			propeptides.add(propeptide);
		} else {
			System.out.println("[CyanobactinsPredictor] "
					+ "Error: propeptide is out of leader peptide bounds."
					+ "Start: " + start + ", end: " + end + ", length: "
					+ length);
		}

		// save information for html output
		Orf orf = OrfAnalyzer.getParentOrf(domain, cluster);
		orf.addPropeptides(propeptides);

		return propeptides;
	}

	public RibosomalClusterTypes type() {
		return RibosomalClusterTypes.CYANOBACTIN;
	}

}
