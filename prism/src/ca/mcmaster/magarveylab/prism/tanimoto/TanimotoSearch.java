package ca.mcmaster.magarveylab.prism.tanimoto;

import java.io.IOException;
import java.sql.Connection;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.similarity.Tanimoto;

import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.database.AccessDatabase;
import ca.mcmaster.magarveylab.prism.database.data.SmallMolecule;
import ca.mcmaster.magarveylab.prism.tanimoto.data.TanimotoScore;
import ca.mcmaster.magarveylab.prism.util.Sorter;
import ca.mcmaster.magarveylab.prism.util.exception.BadSmilesToFingerprinterException;
import ca.mcmaster.magarveylab.prism.util.exception.DatabaseConnectException;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Search the BIRD database for chemically similar small molecules using Tanimoto scores.
 * @author skinnider
 *
 */
public class TanimotoSearch {
	
	/**
	 * Generate Tanimoto scores for a cluster scaffold library using the fcfp6 fingerprinting algorithm.
	 * @param cluster	cluster in question
	 * @param database	filepath of database
	 * @param session	current PRISM session
	 * @return			list of fcfp Tanimoto scores 
	 * @throws InterruptedException 
	 * @throws BadSmilesToFingerprinterException 
	 * @throws CDKException 
	 * @throws IOException 
	 * @throws DatabaseConnectException 
	 * @throws Exception 
	 */
	public static List<TanimotoScore> scoreFingerprints(Cluster cluster, Session session) 
			throws IOException, CDKException, BadSmilesToFingerprinterException, InterruptedException, DatabaseConnectException {
		List<TanimotoScore> scores = new ArrayList<TanimotoScore>();

		List<SmallMolecule> queries = Ecfp6Fcfp6Fingerprinter.getClusterFingerprints(cluster, session);
		List<SmallMolecule> targets = new ArrayList<SmallMolecule>();
		
		// get database fingerprints
		try {
			Connection conn = AccessDatabase.getConnection("jdbc:neo4j://localhost:7474/db/data/cypher");
			targets = AccessDatabase.getSmallMolecules(conn);
		} catch (Exception e) {
			throw new DatabaseConnectException("Error computing Tanimoto coefficients: could not connect to database", e);
		}
		System.out.println("[TanimotoSearch] Got " + targets.size() + " database fingerprints");
		
		// generate Tanimoto scores 
		for (SmallMolecule query : queries) {
			for (SmallMolecule target : targets) {
				TanimotoScore score = score(query, target);
				scores.add(score);
			}
		}
		System.out.println("[TanimotoSearch] Generated " + scores.size() + " Tanimoto coefficients");
		
		// sort Tanimoto coefficients and retain only top 100 (need >10 for duplicates)
		Sorter.sortTanimotoScoresByEcfp6Coefficient(scores);
		if (scores.size() > 100) 
			scores.subList(100, scores.size() - 1).clear();

		return scores;
	}

	/**
	 * Score a predicted cluster scaffold (query) molecule against a known database (target) molecule.
	 * @param query		predicted molecule
	 * @param target	database molecule
	 * @return			a list of Tanimoto coefficients containing different fingerprinting algorithm scores
	 * @throws CDKException
	 */
	public static TanimotoScore score(SmallMolecule query, SmallMolecule target) throws CDKException {
		BitSet queryEcfp6 = query.ecfp6Fingerprint();
		BitSet targetEcfp6 = target.ecfp6Fingerprint();
		BitSet queryFcfp6 = query.fcfp6Fingerprint();
		BitSet targetFcfp6 = target.fcfp6Fingerprint();
		float fcfp6Score = Tanimoto.calculate(queryFcfp6, targetFcfp6);
		float ecfp6Score = Tanimoto.calculate(queryEcfp6, targetEcfp6);
		TanimotoScore score = new TanimotoScore(query, target);
		score.addScore("ecfp6", ecfp6Score);
		score.addScore("fcfp6", fcfp6Score);
		return score;
	}

}
