package ca.mcmaster.magarveylab.prism.motif;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.prism.util.Numbers;

/**
 * Reads output from FIMO motif searches.
 * 
 * @author skinnider
 *
 */
public class FimoSearchReader {

	/**
	 * Motif p-value cutoff, set by default to p = 0.1.
	 */
	private double minimumPValue = 0.1;

	private BufferedReader br;

	/**
	 * Initialize a new FIMO search output output reader.
	 * 
	 * @param br
	 *            buffered reader for the FIMO output file
	 */
	public FimoSearchReader(BufferedReader br) {
		this.br = br;
	}

	/**
	 * Parse FIMO search output and return a list of FIMO search results.
	 * 
	 * @throws NumberFormatException
	 * @throws IOException
	 */
	public List<FimoSearchResult> read() throws NumberFormatException,
			IOException {
		List<FimoSearchResult> results = new ArrayList<FimoSearchResult>();
		
		String line = null;
		boolean nextStart = false;
		while ((line = br.readLine()) != null) {
			if (nextStart) {
				String[] splitLine = line.split("\t");
				if (splitLine.length < 7)
					continue;
				
				String match = splitLine[1];
				int start = Integer.parseInt(splitLine[2]);
				int end = Integer.parseInt(splitLine[3]);
				// String frame = splitLine[4];
				double score = Double.parseDouble(splitLine[5]);
				double eValue = Numbers.newDoubleFromString(splitLine[6]);

				if (eValue < minimumPValue) {
					FimoSearchResult result = new FimoSearchResult(match,
							start, end, score, eValue);
					results.add(result);
				} else {
					System.out.println("[FimoSearchReader] "
							+ "Skipped motif with p-value " + eValue);
				}
			}
			if (line.startsWith("#pattern"))
				nextStart = true;
		}

		br.close();
		return results;
	}

}
