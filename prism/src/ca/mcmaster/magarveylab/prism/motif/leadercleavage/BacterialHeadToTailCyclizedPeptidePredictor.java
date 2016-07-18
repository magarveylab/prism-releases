package ca.mcmaster.magarveylab.prism.motif.leadercleavage;

import ca.mcmaster.magarveylab.enums.RibosomalPrecursorMotifs;
import ca.mcmaster.magarveylab.enums.clusters.RibosomalClusterTypes;
import ca.mcmaster.magarveylab.prism.cluster.analysis.OrfAnalyzer;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.data.Propeptide;
import ca.mcmaster.magarveylab.prism.motif.*;

import java.util.List;
import java.util.ArrayList;

/**
 * Predicts cleavage sites for bacterial head-to-tail cyclized peptides.
 * 
 * @author skinnider
 */
public class BacterialHeadToTailCyclizedPeptidePredictor extends
		AbstractLeaderPredictor implements LeaderPredictor {

	public List<Propeptide> cleaveLeader(Domain domain, Cluster cluster) {
		List<Propeptide> propeptides = new ArrayList<Propeptide>();
		String leader = OrfAnalyzer.getParentOrf(domain, cluster).sequence();

		MotifList motifs = domain.motifs();
		Motif motif = motifs.getBestMotif();

		if (motif == null) {
			System.out.println("[BacterialHeadToTailCyclizedPeptidePredictor] "
					+ "No motifs identified in precursor with sequence "
					+ leader);
		} else { 
			int start = -1;
			if (motif.getType() == RibosomalPrecursorMotifs.Bacterial_HTT_A) {
				System.out.println("[BacterialHeadToTailCyclizedPeptidePredictor] " + 
						"Cleaving clade A bacterial HTT peptide");
				start = motif.getStart() - 5;
			} else {
				if (motifs.contains(RibosomalPrecursorMotifs.Bacterial_HTT_B)
						&& motifs.contains(RibosomalPrecursorMotifs.Bacterial_HTT_C)) {
					System.out.println("[BacterialHeadToTailCyclizedPeptidePredictor] " + 
							"Cleaving clade C1 bacterial HTT peptide");
					motif = motifs.getBestMotif(RibosomalPrecursorMotifs.Bacterial_HTT_C);
					start = motif.getStart() - 3;
				} else if (motif.getType() == RibosomalPrecursorMotifs.Bacterial_HTT_C) {
					System.out.println("[BacterialHeadToTailCyclizedPeptidePredictor] " + 
							"Cleaving clade C2 bacterial HTT peptide");
					start = motif.getStart() - 6;
				}
			}
			
			if (start > 0) {
				String sequence = leader.substring(start);
				Propeptide propeptide = new Propeptide(sequence, start,
						leader.length());
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
		return RibosomalClusterTypes.BACTERIAL_HEAD_TO_TAIL_CYCLIZED;
	}

}