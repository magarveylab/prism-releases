package ca.mcmaster.magarveylab.prism.motif;

import java.io.IOException;
import java.util.List;

import ca.mcmaster.magarveylab.enums.RibosomalPrecursorMotifs;
import ca.mcmaster.magarveylab.enums.clusters.RibosomalClusterTypes;
import ca.mcmaster.magarveylab.prism.util.exception.FimoSearchException;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Propeptide;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Determines the cleavage point of the leader peptide in the precursor sequence
 * of ribosomal peptides.
 * 
 * @author Robyn Edgar, skinnider 
 */

public interface LeaderPredictor {

	/**
	 * Runs a FIMO search using the sequence from the domain that has been
	 * identified as a ribosomal precursor.
	 * 
	 * @param domain
	 *            a ribosomal precursor domain
	 * @param cluster
	 *            the parent cluster
	 * @param session
	 *            the current session
	 * @throws InterruptedException 
	 * @throws IOException
	 */
	public void runFimo(Domain domain, Cluster cluster, Session session)
			throws FimoSearchException, InterruptedException;

	/**
	 * Uses the FIMO results from #runFimo to cleave the leader peptide into the
	 * final aa sequence.
	 * 
	 * @param domain
	 *            a ribosomal precursor domain
	 * @param cluster
	 *            the parent cluster
	 * @return a list of cleaved propeptides
	 */
	public List<Propeptide> cleaveLeader(Domain domain, Cluster cluster);
	
	/**
	 * Get all of the ribosomal peptide precursor cleavage site motifs used by
	 * this class of leader peptide cleavage predictor.
	 * 
	 * @return all of the cleavage site motifs used by this class of leader
	 *         cleavage predictor
	 */
	public List<RibosomalPrecursorMotifs> motifs();

	/**
	 * Get the type of the ribosomal cluster (i.e. biosynthetic family of RiPPs)
	 * associated with this leader predictor.
	 * 
	 * @return	the ribosomal cluster type associated with this leader predictor
	 */
	public RibosomalClusterTypes type(); 
	
}