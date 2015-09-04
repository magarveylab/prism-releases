package ca.mcmaster.magarveylab.prism.database;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.codehaus.jackson.JsonGenerationException;
import org.codehaus.jackson.map.JsonMappingException;
import org.codehaus.jackson.map.ObjectMapper;

import ca.mcmaster.magarveylab.prism.Prism;
import ca.mcmaster.magarveylab.prism.data.Contig;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.wasp.session.Session;
import ca.mcmaster.magarveylab.prism.util.JsonUtils;

/**
 * Generates a JSON report representing the PRISM analysis of a sequence file.
 * 
 * @author prees
 */
public class JsonOutput {

	private Prism prism;
	private Map<String, Object> results;

	/**
	 * Instantiate a new JSON report.
	 * 
	 * @param session
	 *            A session representing the prism analysis
	 */
	public JsonOutput(Session session) {
		this.prism = (Prism) session.webapp();
	}

	/**
	 * Write the JSON report to a file.
	 * 
	 * @param filepath
	 *            A filepath to be created to output the json object to.
	 */
	public void writeToFile(String filepath) throws JsonGenerationException,
	JsonMappingException, IOException {
		ObjectMapper mapper = new ObjectMapper();
		mapper.writeValue(new File(filepath + "/" + prism.genome().file().getName()
				+ ".json"), results);
	}

	/**
	 * Generate the JSON object.
	 */
	public void generateJson() throws JsonGenerationException,
	JsonMappingException, IOException {
		Map<String, Object> prismResult = new HashMap<String, Object>();
		Map<String, Object> prismRun = new HashMap<String, Object>();
		
		Map<String, Object> input = JsonUtils.generateInputJson(prism);
		prismRun.put("input", input);

		Map<String, Object> config = JsonUtils.generateConfigurationJson(prism);
		prismRun.put("configuration", config);

		List<Object> all16s = JsonUtils.generate16sJson(prism);
		if (all16s.size() > 0) {
			prismRun.put("16s_sequences", all16s);
		} else {
			prismRun.put("16s_sequences", null);
		}

		List<Object> clusters = new ArrayList<Object>();
		List<Object> allGlobalDomains = new ArrayList<Object>();
		for (Contig contig : prism.genome().contigs()) {
			List<Object> globalDomains = JsonUtils.generateRegulatoryJson(contig);
			allGlobalDomains.add(globalDomains);

			for (Cluster cluster : contig.clusters()) {
				Map<String, Object> clusterData = JsonUtils.generateClusterJson(prism, cluster);
				clusters.add(clusterData);
			}
		}
		
		prismRun.put("regulatory_genes", allGlobalDomains);
		prismRun.put("clusters", clusters);
		prismResult.put("prism_results", prismRun);

		results = prismResult;
	}
	
}