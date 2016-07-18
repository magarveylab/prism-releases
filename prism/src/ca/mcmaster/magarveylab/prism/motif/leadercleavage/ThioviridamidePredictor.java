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
 * Predicts cleavage sites for thioviridamide.
 * 
 * @author Robyn Edgar, skinnider
 */
public class ThioviridamidePredictor extends AbstractLeaderPredictor implements
		LeaderPredictor {

	public List<Propeptide> cleaveLeader(Domain domain, Cluster cluster) {
		List<Propeptide> propeptides = new ArrayList<Propeptide>();
		String leader = OrfAnalyzer.getParentOrf(domain, cluster).sequence();

		MotifList motifs = domain.motifs();
		Motif motif = motifs.getBestMotif();

		if (motif == null) {
			System.out.println("[ThioviridamidePredictor] "
					+ "No motifs identified in precursor with sequence "
					+ leader);
		} else {
			// remove the N-terminus leader peptide from the precursor
			int start = motif.getEnd();
			if (leader.length() > start) {
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
		return RibosomalClusterTypes.THIOVIRIDAMIDE;
	}

}
