package ca.mcmaster.magarveylab.prism.blast;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.prism.util.Numbers;

/**
 * Reads output of a BLAST search.
 * 
 * @author skinnider
 *
 */
public class BlastSearchReader { 

	private BufferedReader br;

	/**
	 * Instantiate a new BLAST search reader.
	 * 
	 * @param input
	 *            location of the input file to read, i.e. the output file from
	 *            the BLAST search
	 */
	public BlastSearchReader(BufferedReader br) {
		this.br = br;
	}

	/**
	 * Parse the output from a BLAST search into Java data packages.
	 * 
	 * @param readLength
	 *            maximum number of subject hits per query to read. Set to -1 to
	 *            read all subject hits.
	 * @return a list of BLAST search results, each representing the
	 *         intersection of a query and subject
	 * @throws IOException
	 */
	public List<BlastSearchResult> read(int readLength) throws IOException {
		List<BlastSearchResult> results = new ArrayList<BlastSearchResult>();

		String line = null;
		StringBuffer frame = new StringBuffer();
		StringBuffer query = new StringBuffer();
		StringBuffer subject = new StringBuffer();
		double score = 0, eValue = 0, identity = 0, positive = 0;
		int j = 0, length = 0, coverage = 0;
		boolean readingDetails = false;

		while ((line = br.readLine()) != null) {
			if (line.startsWith("Query=")) {
				// reset variables
				query.delete(0, query.length());
				j = 0;

				String queryString = line.split("\\s+")[1];	
				query.append(queryString);
			}
			if (j < readLength || readLength <= 0) {
				if (line.startsWith(">")) {
					subject.append(line.split("\\s+")[1].replace("...", ""));
					readingDetails = true;

					// determine frame
					String[] split = line.split("\\|");
					if (split.length >= 2)
						if (split[2].contains("+")) {
							frame.append("+");
						} else if (split[2].contains("-")) {
							frame.append("-");
						}
				}
				if (readingDetails) {
					if (line.contains("Length="))
						length = Integer.parseInt(line.split("=")[1]);
					if (line.contains("Score = ")) {
						// parse bitscore
						String trimmed = line.trim();
						String[] split = trimmed.split(",");
						String rawScore = split[0].split("\\(")[0].trim().replace("Score = ", "")
								.replace(" bits", "");
						score = Double.parseDouble(rawScore);

						// parse E-value
						String rawValue = split[1].trim().replace("Expect = ", "");
						eValue = Numbers.newDoubleFromString(rawValue);
					}
					if (line.contains("Identities = ")) {
						String[] split = line.trim().split("\\s+");
						String identityRaw = split[2].split("/")[0];
						String positiveRaw = split[6].split("/")[0];
						String coverageRaw = split[2].split("/")[1];

						coverage = Integer.parseInt(coverageRaw);
						identity = Double.parseDouble(identityRaw) / coverage;
						positive = Double.parseDouble(positiveRaw) / coverage;

						readingDetails = false;

						BlastSearchResult result = new BlastSearchResult(query.toString(), subject.toString(), length, 
								score, eValue, identity, positive, coverage);
						results.add(result);
						j++;

						// reset variables
						length = 0;
						subject.delete(0, subject.length());	
						frame.delete(0, frame.length());
					}
				}
			}
		}
		br.close();

		return results;
	}

}
