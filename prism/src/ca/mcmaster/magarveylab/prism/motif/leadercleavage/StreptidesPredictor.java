package ca.mcmaster.magarveylab.prism.motif.leadercleavage;

import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.enums.RibosomalPrecursorMotifs;
import ca.mcmaster.magarveylab.enums.clusters.RibosomalClusterTypes;
import ca.mcmaster.magarveylab.prism.cluster.analysis.OrfAnalyzer;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.data.Propeptide;
import ca.mcmaster.magarveylab.prism.motif.AbstractLeaderPredictor;
import ca.mcmaster.magarveylab.prism.motif.LeaderPredictor;
import ca.mcmaster.magarveylab.prism.motif.Motif;
import ca.mcmaster.magarveylab.prism.motif.MotifList;

/**
 * Predicts cleavage sites for streptides.
 * 
 * @author skinnider
 *
 */
public class StreptidesPredictor extends AbstractLeaderPredictor implements
		LeaderPredictor {

	public List<Propeptide> cleaveLeader(Domain domain, Cluster cluster) {
		List<Propeptide> propeptides = new ArrayList<Propeptide>();
		String leader = OrfAnalyzer.getParentOrf(domain, cluster).sequence();

		MotifList motifs = domain.motifs();
		Motif head = motifs.getBestMotif(RibosomalPrecursorMotifs.Streptide_N);
		Motif tail = motifs.getBestMotif(RibosomalPrecursorMotifs.Streptide_C);
		
		if (head == null) {
			System.out.println("[StreptidesPredictor] "
					+ "No head motif identified in precursor with sequence "
					+ leader);
		} else if (tail == null) {
			// use the head only for precursor cleavage 
			int start = head.getStart() + 13;
			
			if (leader.length() > start) {
				String sequence = leader.substring(start);
				Propeptide propeptide = new Propeptide(sequence, start,
						leader.length());
				propeptide.addMotif(head);
				propeptides.add(propeptide);
			}
		} else {
			// use both head and tail for precursor cleavage 
			int start = head.getStart() + 13;
			int end = tail.getStart() + 1;
			if (leader.length() > start) {
				if (leader.length() > end && start < end && end - start >= 5) {
					String sequence = leader.substring(start, end);
					Propeptide propeptide = new Propeptide(sequence, start, end);
					propeptide.addMotif(head);
					propeptide.addMotif(tail);
					propeptides.add(propeptide);
				} else {
					String sequence = leader.substring(start);
					Propeptide propeptide = new Propeptide(sequence, start,
							leader.length());
					propeptide.addMotif(head);
					propeptide.addMotif(tail);
					propeptides.add(propeptide);
				}
			}
		}

		// save information for html output
		Orf orf = OrfAnalyzer.getParentOrf(domain, cluster);
		orf.addPropeptides(propeptides);

		return propeptides;
	}

	public RibosomalClusterTypes type() {
		return RibosomalClusterTypes.STREPTIDE;
	}

}