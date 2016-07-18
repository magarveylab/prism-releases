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
 * Predicts cleavage sites for linaridins.
 * 
 * @author skinnider, Robyn Edgar
 */
public class LinaridinsPredictor extends AbstractLeaderPredictor implements
		LeaderPredictor {

	public List<Propeptide> cleaveLeader(Domain domain, Cluster cluster) {
		List<Propeptide> propeptides = new ArrayList<Propeptide>();
		String leader = OrfAnalyzer.getParentOrf(domain, cluster).sequence();

		MotifList motifs = domain.motifs();
		Motif motif = motifs.getBestMotif();
		if (motif == null) {
			System.out.println("[LinaridinsPredictor] "
					+ "No motifs identified in precursor "
					+ "with sequence " + leader);
		} else {
			int start = -1;
			if (motif.getType() == RibosomalPrecursorMotifs.Linaridin) {
				start = motif.getEnd() - 3;
			} else if (motif.getType() == RibosomalPrecursorMotifs.Linaridin_2) {
				start = motif.getEnd() - 3;
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
		}

		return propeptides;
	}
	
	public RibosomalClusterTypes type() {
		return RibosomalClusterTypes.LINARIDIN;
	}

}