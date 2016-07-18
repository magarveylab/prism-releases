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
 * Predicts cleavage sites for lantipeptides.
 * 
 * @author skinnider
 */
public class LantipeptidesPredictor extends AbstractLeaderPredictor implements
		LeaderPredictor {

	public List<Propeptide> cleaveLeader(Domain domain, Cluster cluster) {
		List<Propeptide> propeptides = new ArrayList<Propeptide>();
		String leader = OrfAnalyzer.getParentOrf(domain, cluster).sequence();

		MotifList motifs = domain.motifs();
		Motif motif = motifs.getBestMotif();

		if (motif == null) {
			System.out.println("[LantipeptidesPredictor] "
					+ "No motifs identified in precursor "
					+ "with sequence " + leader);
		} else {
			// this is the most reliable motif for predicting cleavage 
			if (motifs.contains(RibosomalPrecursorMotifs.LanA_B)
					&& !motifs.contains(RibosomalPrecursorMotifs.LanA_K)) 
				motif = motifs.getBestMotif(RibosomalPrecursorMotifs.LanA_B);
			
			int start = -1;
			RibosomalPrecursorMotifs type = motif.getType();
			System.out.println("[LantipeptidesPredictor] " + 
						"Cleaving " + type.toString() + " lantipeptide");
			switch (type) {
				case LanA_A: 
					start = motif.getEnd() - 4;
					break;
				case LanA_B:
					start = motif.getEnd();
					break;
				case LanA_C:
					start = motif.getStart() + 7;
					break;
				case LanA_D:
					start = motif.getEnd() + 9;
					break; 
				case LanA_E: 
					start = motif.getStart() + 8;
					break;
				case LanA_F:
					start = motif.getEnd() + 4;
					break;
				case LanA_G:
					start = motif.getEnd() - 9;
					break;
				case LanA_H:
					start = motif.getEnd() + 2;
					break;
				case LanA_I:
					start = motif.getEnd() + 4;
					break;
				case LanA_J1:
					start = motif.getEnd() - 7;
					break;
				case LanA_J2:
					start = motif.getStart() + 18;
					break;
				case LanA_K:
					start = motif.getEnd() - 5;
					break;
				case LanA_L:
					start = motif.getEnd() - 3;
					break;
				default:
					break;
			}
			
			if (start >= 0 && leader.length() > start) {
				String sequence = leader.substring(start);
				Propeptide propeptide = new Propeptide(sequence, start,
						leader.length());
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
		return RibosomalClusterTypes.LANTIPEPTIDE;
	}
	
}
