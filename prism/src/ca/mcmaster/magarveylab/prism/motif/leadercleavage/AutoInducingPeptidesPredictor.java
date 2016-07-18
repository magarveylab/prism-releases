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
 * Predicts cleavage sites for auto-inducing peptides.
 * 
 * @author skinnider
 */

public class AutoInducingPeptidesPredictor extends AbstractLeaderPredictor
		implements LeaderPredictor {

	public List<Propeptide> cleaveLeader(Domain domain, Cluster cluster) {
		List<Propeptide> propeptides = new ArrayList<Propeptide>();
		String leader = OrfAnalyzer.getParentOrf(domain, cluster).sequence();

		MotifList motifs = domain.motifs();
		Motif motif = motifs.getBestMotif();

		if (motif == null) {
			System.out.println("[AutoInducingPeptidesPredictor] "
					+ "No motifs identified in precursor with sequence "
					+ leader);
		} else {
			int length = leader.length();
			int end = (length > motif.getStart() + 4) ? motif.getStart() + 4
					: length;
			int[] starts = new int[] { motif.getStart() - 3,
					motif.getStart() - 4, motif.getStart() - 5 };
			
			for (int start : starts) {
				if (start < 0)
					continue;
				String sequence = leader.substring(start, end);
				Propeptide propeptide = new Propeptide(sequence, start, end);
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
		return RibosomalClusterTypes.AUTO_INDUCING_PEPTIDE;
	}
	
}
