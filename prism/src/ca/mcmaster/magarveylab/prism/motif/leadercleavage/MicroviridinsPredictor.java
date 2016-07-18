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
 * Predicts cleavage sites for microviridins.
 * 
 * @author skinnider
 */
public class MicroviridinsPredictor extends AbstractLeaderPredictor implements
		LeaderPredictor {

	public List<Propeptide> cleaveLeader(Domain domain, Cluster cluster) {
		List<Propeptide> propeptides = new ArrayList<Propeptide>();
		String leader = OrfAnalyzer.getParentOrf(domain, cluster).sequence();

		MotifList motifs = domain.motifs();
		Motif motif = motifs
				.getFirstMotif(RibosomalPrecursorMotifs.Microviridin);

		if (motif == null) {
			System.out.println("[MicroviridinsPrecursor] "
					+ "No microviridin motifs identified in precursor "
					+ "with sequence " + leader);
		} else {
			int start = motif.getStart() - 3;
			int end = motif.getEnd() + 1 > leader.length() ? leader.length()
					: motif.getEnd() + 1;
			if (start > 0) {
				String sequence = leader.substring(start, end);
				Propeptide propeptide = new Propeptide(sequence, start, end);
				propeptide.addMotif(motif);
				propeptides.add(propeptide);
			}
		}

		// save information for html output
		Orf orf = OrfAnalyzer.getParentOrf(domain, cluster);
		orf.addPropeptides(propeptides);

		return propeptides;
	}

	public RibosomalClusterTypes type() {
		return RibosomalClusterTypes.MICROVIRIDIN;
	}

}
