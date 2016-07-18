package ca.mcmaster.magarveylab.prism.motif;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.enums.RibosomalPrecursorMotifs;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.util.exception.FimoSearchException;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Contains core functionality used by leader peptide cleavage site predictors
 * irrespective of ribosomal peptide class.
 * 
 * @author skinnider
 *
 */
public abstract class AbstractLeaderPredictor implements LeaderPredictor {

	public void runFimo(Domain domain, Cluster cluster, Session session)
			throws FimoSearchException, InterruptedException {
		for (RibosomalPrecursorMotifs motif : motifs()) {
			try {
				MotifSearch fimo = new MotifSearch(domain, motif, cluster,
						session);
				fimo.run();
			} catch (IOException e) {
				throw new FimoSearchException("Error executing FIMO search: "
						+ e.getMessage());
			}
		}
	}

	public List<RibosomalPrecursorMotifs> motifs() {
		List<RibosomalPrecursorMotifs> motifs = new ArrayList<RibosomalPrecursorMotifs>();
		for (RibosomalPrecursorMotifs motif : RibosomalPrecursorMotifs.values()) {
			if (motif.getClusterType().equals(type())) {
				motifs.add(motif);
			}
		}
		return motifs;
	}
	
}
