package ca.mcmaster.magarveylab.prism.tanimoto;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.exception.CDKException;

import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.database.data.SmallMolecule;
import ca.mcmaster.magarveylab.prism.tanimoto.data.Fingerprint;
import ca.mcmaster.magarveylab.prism.util.exception.BadSmilesToFingerprinterException;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Generate fingerprints using the ecfp6 and fcfp6 algorithms. Requires Python
 * and RDKit.
 * 
 * @author skinnider, cDejong
 * 
 */
public class Ecfp6Fcfp6Fingerprinter {

	private final static String script = "ecfp6_fcfp6.py";

	/**
	 * Get fingerprints for a cluster scaffold library.
	 * 
	 * @param cluster
	 *            cluster in question
	 * @param session
	 *            current session 
	 * @return fingerprints for all scaffolds in the cluster library
	 * @throws IOException
	 * @throws CDKException
	 * @throws BadSmilesToFingerprinterException
	 * @throws InterruptedException
	 */
	public static List<SmallMolecule> getClusterFingerprints(Cluster cluster, Session session) 
			throws IOException, CDKException, BadSmilesToFingerprinterException, InterruptedException {
		List<SmallMolecule> molecules = new ArrayList<SmallMolecule>();
		Map<String, String> library = FingerprintUtil.getLibraryForFingerprinting(cluster);

		// set permissions of fcfp6 script
		chmodFingerprinter(session);

		for (Map.Entry<String, String> entry : library.entrySet()) {
			SmallMolecule molecule = new SmallMolecule();

			String name = entry.getKey();
			String smiles = entry.getValue();

			molecule.setSmiles(smiles);
			molecule.setName(name);

			String dir = session.subDir("tanimoto");
			List<Fingerprint> fingerprints = getCompoundFingerprints(smiles,
					dir);
			for (Fingerprint fingerprint : fingerprints) {
				if (fingerprint.type().equals("fcfp6")) {
					molecule.setFcfp6Fingerprint(fingerprint.bitset());
				} else if (fingerprint.type().equals("ecfp6")) {
					molecule.setEcfp6Fingerprint(fingerprint.bitset());
				}
			}

			// ignore molecules with null fingerprints
			if (fingerprints.size() < 2)
				continue;

			molecules.add(molecule);
		}

		return molecules;
	}

	/**
	 * Get the fingerprints of a single cluster scaffold library compound.
	 * 
	 * @param smiles
	 *            SMILES of the compound in question
	 * @param session
	 *            current PRISM session
	 * @return compound fingerprints
	 * @throws IOException
	 * @throws BadSmilesToFingerprinterException
	 * @throws InterruptedException
	 */
	public static List<Fingerprint> getCompoundFingerprints(String smiles,
			String dir) throws IOException, BadSmilesToFingerprinterException,
			InterruptedException {
		// run python script
		ProcessBuilder pb = new ProcessBuilder("python", script, smiles);
		Map<String, String> environment = pb.environment();
		environment.put("PYTHONPATH", "/opt/rdkit/");
		environment.put("LD_LIBRARY_PATH", "/opt/rdkit/lib/");
		pb.directory(new File(dir));
		Process p = pb.start();
		p.waitFor();

		List<Fingerprint> fingerprints = FingerprintUtil
				.outputToFingerprints(p);
		return fingerprints;
	}

	/**
	 * Set permissions of ecfp6/fcfp6 python script.
	 * 
	 * @param session
	 *            current PRISM session
	 * @throws IOException
	 */
	private static void chmodFingerprinter(Session session) throws IOException {
		// chmod fcfp6 executable
		String dir = session.subDir("tanimoto");
		String executable = dir + script;
		Runtime.getRuntime().exec("chmod +x " + executable);
	}

}
