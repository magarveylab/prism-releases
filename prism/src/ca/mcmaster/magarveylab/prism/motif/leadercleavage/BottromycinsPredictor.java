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
 * Predicts cleavage sites for bottromycins.
 * 
 * @author Robyn Edgar, skinnider
 */
public class BottromycinsPredictor extends AbstractLeaderPredictor implements
		LeaderPredictor {

	public List<Propeptide> cleaveLeader(Domain domain, Cluster cluster) {
		List<Propeptide> propeptides = new ArrayList<Propeptide>();
		String leader = OrfAnalyzer.getParentOrf(domain, cluster).sequence();

		MotifList motifs = domain.motifs();
		Motif motif = motifs.getBestMotif();

		if (motif == null) {
			System.out.println("[BottromycinsPredictor] "
					+ "No motifs identified in precursor with sequence "
					+ leader);
		} else {
			// always cleave the N-terminal 'M' for bottromycins, and
			// remove the C-terminus leader peptide from the precursor
			if (leader.length() >= motif.getStart()) {
				String sequence = leader.substring(1, motif.getStart() - 1);
				Propeptide propeptide = new Propeptide(sequence, 1,
						motif.getStart() - 1);
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
		return RibosomalClusterTypes.BOTTROMYCIN;
	}

}
