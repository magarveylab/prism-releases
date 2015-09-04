package ca.mcmaster.magarveylab.prism.tanimoto;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Library;
import ca.mcmaster.magarveylab.prism.tanimoto.data.Fingerprint;

/**
 * Utilities class for operations involving chemical fingerprints.
 * 
 * @author cDejong, skinnider
 *
 */
public class FingerprintUtil {

	/**
	 * Read a string output by the BitSet.toString() method back in to a BitSet.
	 * 
	 * @param in
	 *            input string
	 * @param size
	 *            bitset size
	 * @return a BitSet corresponding to the input string
	 */
	public static BitSet bitset(String in, int size) {
		String array = in.replace("{", "").replace("}", "");
		String[] bits = array.split(",");
		BitSet bitset = new BitSet(size);
		for (String set : bits) {
			String bit = set.trim();
			int index = Integer.parseInt(bit);
			if (index > size)
				throw new IllegalArgumentException("Error: fingerprint bitset greater than " + size);
			bitset.set(index);
		}
		return bitset;
	}

	/**
	 * Convert the output from Chris's RDKit python script to a set of
	 * fingerprint objects, including the fingerprinting algorithm used, and the
	 * fingerprint bitset itself.
	 * 
	 * @param p
	 *            python process
	 * @return fingerprints
	 * @throws IOException
	 */
	public static List<Fingerprint> outputToFingerprints(Process p)
			throws IOException {
		List<Fingerprint> fingerprints = new ArrayList<Fingerprint>();

		BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
		String line = null;
		while ((line = in.readLine()) != null) {
			String parts[] = line.split(" ");
			String type = parts[0];
			BitSet bitset = bitSetFromString(parts[1]);
			Fingerprint fingerprint = new Fingerprint(type, bitset);
			fingerprints.add(fingerprint);
		}
		in.close();

		return fingerprints;
	}

	public static BitSet bitSetFromString(String s) {
		BitSet t = new BitSet(s.length());
		int lastBitIndex = s.length() - 1;
		int i = lastBitIndex;
		while (i >= 0) {
			if (s.charAt(i) == '1') {
				t.set(lastBitIndex - i);
				i--;
			} else
				i--;
		}
		return t;
	}

	/**
	 * Convert a cluster scaffold library object into a map of name-SMILES pairs
	 * for the generation of chemical fingerprints.
	 * 
	 * @param cluster
	 *            cluster to convert
	 * @return
	 */
	public static Map<String, String> getLibraryForFingerprinting(Cluster cluster) {
		HashMap<String, String> input = new HashMap<String, String>();
		Library library = cluster.library();
		List<String> scaffolds = library.scaffolds();

		// convert scaffold library into name-smiles pairs
		for (int i = 0; i < scaffolds.size(); i++) {
			String scaffold = scaffolds.get(i);
			String name = "Scaffold_" + (i + 1);
			input.put(name, scaffold);
		}
		System.out.println("[FingerprintUtil] Generated " + input.keySet().size() 
				+ " scaffold-smiles " + "pre-fingerprint pairs");
		return input;
	}

}
