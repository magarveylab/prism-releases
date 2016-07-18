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
 * Predicts cleavage sites for thiopeptides.
 * 
 * @author skinnider
 */
public class ThiopeptidesPredictor extends AbstractLeaderPredictor implements
		LeaderPredictor {

	public List<Propeptide> cleaveLeader(Domain domain, Cluster cluster) {
		List<Propeptide> propeptides = new ArrayList<Propeptide>();
		String leader = OrfAnalyzer.getParentOrf(domain, cluster).sequence();

		MotifList motifs = domain.motifs();
		Motif motif = motifs.getBestMotif();

		if (motif == null) {
			System.out.println("[ThiopeptidesPredictor] "
					+ "No motifs identified in precursor with sequence "
					+ leader);
		} else {
			int start = -1;
			int end = -1;
			if (motif.getType() == RibosomalPrecursorMotifs.Thiopeptide_A) {
				System.out.println("[ThiopeptidesPredictor] "
						+ "Using motif Thiopeptide_A");
				// thiostrepton/siomycin-type thiopeptides
				start = motif.getStart() + 32;
			} else if (motif.getType() == RibosomalPrecursorMotifs.Thiopeptide_B) {
				System.out.println("[ThiopeptidesPredictor] "
						+ "Using motif Thiopeptide_B");
				// nosiheptide/nocathiacin-type thiopeptides
				start = motif.getStart() + 33;
			} else if (motif.getType() == RibosomalPrecursorMotifs.Thiopeptide_C) {
				System.out.println("[ThiopeptidesPredictor] "
						+ "Using motif Thiopeptide_C");
				// large thiopeptide clade #1 
				start = motif.getEnd() - 3;
			} else if (motif.getType() == RibosomalPrecursorMotifs.Thiopeptide_D) {
				System.out.println("[ThiopeptidesPredictor] "
						+ "Using motif Thiopeptide_D");
				// thiopeptides w/ C-terminal cleavage (GE2270/thiomuracin-type)
				start = motif.getStart() + 29;
				end = motif.getEnd() + 3;
			} else if (motif.getType() == RibosomalPrecursorMotifs.Thiopeptide_E) {
				System.out.println("[ThiopeptidesPredictor] "
						+ "Using motif Thiopeptide_E");
				// large thiopeptide clade #2 
				start = motif.getStart() + 4;
			}

			// remove the N-terminus leader peptide from the precursor
			if (start > -1 && leader.length() > start) {
				if (end > -1 && leader.length() >= end) {
					String sequence = leader.substring(start, end);
					Propeptide propeptide = new Propeptide(sequence, start,
							end);
					propeptide.addMotif(motif);
					propeptides.add(propeptide);
				} else {
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
		}

		return propeptides;

	}

	public RibosomalClusterTypes type() {
		return RibosomalClusterTypes.THIOPEPTIDE;
	}

}
