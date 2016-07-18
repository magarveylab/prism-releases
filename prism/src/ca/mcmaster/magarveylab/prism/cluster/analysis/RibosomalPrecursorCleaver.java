package ca.mcmaster.magarveylab.prism.cluster.analysis;

import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Propeptide;
import ca.mcmaster.magarveylab.prism.enums.LeaderPredictors;
import ca.mcmaster.magarveylab.prism.motif.LeaderPredictor;
import ca.mcmaster.magarveylab.prism.util.exception.FimoSearchException;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Cleaves putative RiPP precursor peptides.
 * 
 * @author skinnider
 *
 */
public class RibosomalPrecursorCleaver {

	/**
	 * Predict cleavage sites for all ribosomally synthesized and
	 * post-translationally modified peptides within a cluster, generating
	 * cleaved propeptides.
	 * 
	 * @param cluster
	 *            cluster being analyzed
	 * @param session
	 *            the current session
	 * @throws FimoSearchException
	 * @throws InterruptedException
	 */
	public static void cleave(Cluster cluster, Session session)
			throws FimoSearchException, InterruptedException {
		// cleave leader sequence from precursor for a final propeptide
		System.out.println(
				"[RibosomalPrecursorCleaver] Getting propeptide module permutations");
		for (Domain domain : RibosomalPrecursorAnalyzer
				.getPrecursors(cluster)) {
			List<Propeptide> propeptides = new ArrayList<Propeptide>();
			LeaderPredictor predictor = null;
			LeaderPredictors[] predictors = LeaderPredictors.values();
			for (LeaderPredictors p : predictors)
				if (p.domain() == domain.type())
					predictor = p.predictor();
			if (predictor == null) {
				System.out.println("[RibosomalPrecursorCleaver] Error: "
						+ "could not get leader "
						+ "cleavage site(s) predictor for domain "
						+ domain.type());
				continue;
			} else {
				System.out.println("Got " + predictor.toString()
						+ " leader peptide cleavage site(s) predictor");
			}

			// run FIMO
			predictor.runFimo(domain, cluster, session);

			// create cleaved propeptides
			propeptides.addAll(predictor.cleaveLeader(domain, cluster));

			// print propeptides to command line
			if (propeptides.isEmpty()) {
				System.out.println("[RibosomalPrecursorCleaver] Error: "
						+ "could not identify any propeptides");
			} else {
				StringBuffer sb = new StringBuffer();
				for (Propeptide p : propeptides)
					sb.append(p.getSequence() + " ");
				System.out.println("[RibosomalPrecursorCleaver] "
						+ "Identified " + propeptides.size() + " propeptides: "
						+ sb.toString());
			}
		}
	}

}