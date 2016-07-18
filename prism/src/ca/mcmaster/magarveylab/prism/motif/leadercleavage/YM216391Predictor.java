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
 * Predicts cleavage sites for YM-family peptides.
 * 
 * @author Robyn Edgar, skinnider
 */
public class YM216391Predictor extends AbstractLeaderPredictor implements
		LeaderPredictor {

	public List<Propeptide> cleaveLeader(Domain domain, Cluster cluster) {
		List<Propeptide> propeptides = new ArrayList<Propeptide>();
		String leader = OrfAnalyzer.getParentOrf(domain, cluster).sequence();

		MotifList motifs = domain.motifs();
		Motif lead = motifs.getBestMotif(RibosomalPrecursorMotifs.YM216391_N);
		Motif tail = motifs.getBestMotif(RibosomalPrecursorMotifs.YM216391_C);

		if (lead == null) {
			System.out.println("[YM216391] "
					+ "No leader motifs identified in precursor with sequence "
					+ leader);
		} else if (tail == null) {
			System.out.println("[YM216391] "
					+ "No tail motifs identified in precursor with sequence "
					+ leader);
		} else if (tail.getStart() - 8 > lead.getEnd()) {
			if (leader.length() > tail.getStart() - 8) {
				String sequence = leader.substring(lead.getEnd(),
						tail.getStart() - 8);
				Propeptide propeptide = new Propeptide(sequence, lead.getEnd(),
						tail.getStart() - 8);
				propeptide.addMotif(lead);
				propeptide.addMotif(tail);
				propeptides.add(propeptide);

				// save information for html output
				Orf orf = OrfAnalyzer.getParentOrf(domain, cluster);
				orf.addPropeptides(propeptides);
			}
		}

		return propeptides;
	}

	public RibosomalClusterTypes type() {
		return RibosomalClusterTypes.YM216391;
	}

}