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
 * Predicts cleavage sites for sactipeptides.
 * 
 * @author skinnider
 */
public class SactipeptidesPredictor extends AbstractLeaderPredictor implements
		LeaderPredictor {

	public List<Propeptide> cleaveLeader(Domain domain, Cluster cluster) {
		List<Propeptide> propeptides = new ArrayList<Propeptide>();
		String leader = OrfAnalyzer.getParentOrf(domain, cluster).sequence();

		MotifList motifs = domain.motifs();
		Motif motif = motifs.getBestMotif();

		int start = -1;
		if (motif == null) {
			System.out.println("[SactipeptidesPredictor] "
					+ "No sactipeptide motifs identified in precursor "
					+ "with sequence " + leader);
		} else if (motif.getType() == RibosomalPrecursorMotifs.Sactipeptide_A) {
			System.out.println("[SactipeptidesPredictor] "
					+ "Cleaving clade A sactipeptide");
			start = motif.getEnd() - 1;
		} else if (motif.getType() == RibosomalPrecursorMotifs.Sactipeptide_B) {
			System.out.println("[SactipeptidesPredictor] "
					+ "Cleaving clade B sactipeptide");
			start = motif.getStart() + 5;
		} else if (motif.getType() == RibosomalPrecursorMotifs.Sactipeptide_C) {
			System.out.println("[SactipeptidesPredictor] "
					+ "Cleaving clade C sactipeptide");
			start = motif.getEnd() - 2;
		} else if (motif.getType() == RibosomalPrecursorMotifs.Sactipeptide_D) {
			System.out.println("[SactipeptidesPredictor] "
					+ "Cleaving clade D sactipeptide");
			start = motif.getEnd() - 7;
		} else if (motif.getType() == RibosomalPrecursorMotifs.Sactipeptide_E) {
			System.out.println("[SactipeptidesPredictor] "
					+ "Cleaving clade E sactipeptide");
			start = motif.getStart() + 9;
		}

		if (start > -1 && leader.length() > start) {
			String sequence = leader.substring(start);
			Propeptide propeptide = new Propeptide(sequence, start,
					leader.length());
			propeptide.addMotif(motif);
			propeptides.add(propeptide);
		}

		// save information for html output
		Orf orf = OrfAnalyzer.getParentOrf(domain, cluster);
		orf.addPropeptides(propeptides);

		return propeptides;
	}

	public RibosomalClusterTypes type() {
		return RibosomalClusterTypes.SACTIPEPTIDE;
	}

}
