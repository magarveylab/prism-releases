package ca.mcmaster.magarveylab.prism.genome;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.prism.genome.data.HmmSearchResult;
import ca.mcmaster.magarveylab.prism.genome.data.HmmSearchResultAnnotation;
import ca.mcmaster.magarveylab.prism.util.Numbers;

/**
 * Read output from a hidden Markov model search.
 * @author skinnider
 *
 */
public class HmmSearchReader {

	private BufferedReader br;

	/**
	 * Initialize a new hmmsearch output reader.
	 * 
	 * @param session
	 *            the current Session
	 */
	public HmmSearchReader(BufferedReader br) {
		this.br = br;
	}

	/**
	 * Parse hmmsearch output and return a list of {@code HmmSearchResult}s with domains annotated.
	 * @return	all HmmSearchResults
	 * @throws IOException
	 */
	public List<HmmSearchResult> read() throws IOException {
		List<HmmSearchResult> results = new ArrayList<HmmSearchResult>();

		boolean readingScoreTable = false;
		String line = null;
		String name = null;
		Integer counter1 = null, counter2 = null; // counter1 detects table title, counter2 detects score details line
		while ((line = br.readLine()) != null) {
			// Read scores table
			if (counter1 != null) counter1++;
			if (line.indexOf("Scores for complete sequences") != -1) {
				// the title of the table has been detected -- read 3 lines then start parsing
				counter1 = 0;
				readingScoreTable = true;
			}
			if (counter1 != null && counter1 > 3) {
				line = line.trim();
				if (line.indexOf("------ inclusion threshold ------") != -1 || line.equals("") || line == null) 
					readingScoreTable = false;
				if (readingScoreTable) {
					HmmSearchResult row = newHmmSearchResultFromRow(line);
					results.add(row);	
				}
			}

			// Match annotations to rows
			if (counter2 != null) counter2++;
			if (line.indexOf(">>") != -1) {
				name = line.split("\\s+")[1];
				counter2 = 0;
				// start of domain annotations -- read 2 line header then parse table
			}
			if (counter2 != null && counter2 > 2) {
				if (line.indexOf('?') != -1)
					continue;
				if (line.indexOf('!') != -1 && line.indexOf("!!") == -1) {
					//	double independentEValue = convert(line.split("\\s+")[5]);
					double score = Numbers.newDoubleFromString(line.split("\\s+")[2]);
					int envfrom = Integer.parseInt(line.split("\\s+")[12]);
					int envto = Integer.parseInt(line.split("\\s+")[13]);
					for (HmmSearchResult row : results) {
						if (row.name.equals(name)) {
							HmmSearchResultAnnotation annotation = new HmmSearchResultAnnotation(envfrom, envto, score, row);
							row.addAnnotation(annotation);
						}
					}
				}
			}
		}
		br.close();

		return results;
	}

	/**
	 * Parses a single row of a hmmsearch reading score table to instantiate a new <code>HmmRow</code>.
	 * @param columns	the reading score table line, represented 
	 * @return
	 */
	public HmmSearchResult newHmmSearchResultFromRow(String line) {
		String[] columns = line.split("\\s+");

		String fullSequenceEValue = null;
		String fullSequenceScore = null;
		String fullSequenceBias = null;
		String bestDomainEValue = null;
		String bestDomainScore = null;
		String bestDomainBias = null;
		String sequence = null;
		String description = null;

		if (columns.length > 0 && columns[0] != null)
			fullSequenceEValue = columns[0];
		if (columns.length > 1 && columns[1] != null) {
			fullSequenceScore = columns[1];
		}
		if (columns.length > 2 && columns[2] != null) {
			fullSequenceBias = columns[2];
		}
		if (columns.length > 3 && columns[3] != null) {
			bestDomainEValue = columns[3];
		}
		if (columns.length > 4 && columns[4] != null) {
			bestDomainScore = columns[4];
		}
		if (columns.length > 5 && columns[5] != null) {
			bestDomainBias = columns[5];
		}
		// 6 & 7 unimportant
		if (columns.length > 8 && columns[8] != null) {
			sequence = columns[8];
		}
		if (columns.length > 9 && columns[9] != null) {
			description = columns[9];
		}

		HmmSearchResult result = new HmmSearchResult(fullSequenceEValue, fullSequenceScore, fullSequenceBias,
				bestDomainEValue, bestDomainScore, bestDomainBias, sequence, description);
		return result;
	}

}
