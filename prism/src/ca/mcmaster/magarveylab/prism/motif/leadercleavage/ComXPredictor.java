package ca.mcmaster.magarveylab.prism.motif.leadercleavage;

import java.util.List;
import java.util.ArrayList;

import ca.mcmaster.magarveylab.enums.clusters.RibosomalClusterTypes;
import ca.mcmaster.magarveylab.prism.cluster.analysis.OrfAnalyzer;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.data.Propeptide;
import ca.mcmaster.magarveylab.prism.motif.*;

/**
 * Predicts cleavage sites for ComX peptides.
 * 
 * @author skinnider
 */
public class ComXPredictor extends AbstractLeaderPredictor implements
		LeaderPredictor {

	public List<Propeptide> cleaveLeader(Domain domain, Cluster cluster) {
		List<Propeptide> propeptides = new ArrayList<Propeptide>();
		String leader = OrfAnalyzer.getParentOrf(domain, cluster).sequence();
		
		MotifList motifs = domain.motifs();
		Motif motif = motifs.getBestMotif();

		if (motif == null) {
			System.out.println("[ComXPredictor] "
					+ "No motifs identified in precursor with sequence "
					+ leader);
		} else {
			int end = motif.getEnd() + 11;
			int length = leader.length();
			
			// cleavage is also known to occur at +6 or +9
			// if +11 is impossible/implausible, pick one of these 
			if (end > length || length - end < 6) 
				end = motif.getEnd() + 9;
			if (end > length || length - end < 6)
				end = motif.getEnd() + 6;
			
			if (length > end) {
				// remove the N-terminus leader peptide from the precursor
				String sequence = leader.substring(end);
				Propeptide propeptide = new Propeptide(sequence, end, length);
				propeptide.addMotif(motif);
				propeptides.add(propeptide);

				// save information for html output
				Orf orf = OrfAnalyzer.getParentOrf(domain, cluster);
				orf.addPropeptides(propeptides);
			}
		}

		return propeptides;

	}

	public RibosomalClusterTypes type() {
		return RibosomalClusterTypes.COMX;
	}

}
