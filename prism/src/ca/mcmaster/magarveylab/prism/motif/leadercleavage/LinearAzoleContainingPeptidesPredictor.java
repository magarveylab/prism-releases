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
 * Predicts cleavage sites for linear azole-containing peptides.
 * 
 * @author skinnider
 */
public class LinearAzoleContainingPeptidesPredictor extends
		AbstractLeaderPredictor implements LeaderPredictor {

	public List<Propeptide> cleaveLeader(Domain domain, Cluster cluster) {
		List<Propeptide> propeptides = new ArrayList<Propeptide>();
		String leader = OrfAnalyzer.getParentOrf(domain, cluster).sequence();

		MotifList motifs = domain.motifs();
		Motif motif = motifs.getBestMotif();
		
		if (motifs.size() == 0) {
			System.out.println("[LinearAzoleContainingPeptidesPrecursor] "
					+ "No motifs identified in precursor "
					+ "with sequence " + leader);
		} else {
			int start = -1; 
			if (motif.getType() == RibosomalPrecursorMotifs.LAP_1) {
				start = motif.getEnd() + 3; 
			} else if (motif.getType() == RibosomalPrecursorMotifs.LAP_2) {
				start = motif.getEnd() + 5;
			} else if (motif.getType() == RibosomalPrecursorMotifs.LAP_3) {
				start = motif.getEnd() + 1;
			} else if (motif.getType() == RibosomalPrecursorMotifs.LAP_4) {
				start = motif.getEnd() - 5;
			}
			
			if (start < 0)
				start = 0;
			int length = leader.length();
			if (length > start) {
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
		return RibosomalClusterTypes.LINEAR_AZOLE_CONTAINING_PEPTIDE;
	}

}