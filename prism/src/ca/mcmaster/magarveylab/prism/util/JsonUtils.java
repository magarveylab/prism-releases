package ca.mcmaster.magarveylab.prism.util;

import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.codehaus.jackson.JsonGenerationException;
import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.map.JsonMappingException;
import org.codehaus.jackson.map.ObjectMapper;
import org.codehaus.jackson.type.TypeReference;

import ca.mcmaster.magarveylab.enums.ClusterFamilies;
import ca.mcmaster.magarveylab.enums.DomainFamilies;
import ca.mcmaster.magarveylab.enums.clusters.ClusterType;
import ca.mcmaster.magarveylab.prism.Prism;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Contig;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.data.Organism;
import ca.mcmaster.magarveylab.prism.data.RnaSequence;
import ca.mcmaster.magarveylab.prism.data.sugar.Sugar;
import ca.mcmaster.magarveylab.prism.database.Utils;
import ca.mcmaster.magarveylab.prism.web.PrismConfig;

/**
 * @author prees
 *
 */
public class JsonUtils {

	/**
	 * Given a Map object and file path will convert the map to a JSON object
	 * and write it to the file.
	 * 
	 * @param results
	 *            A map object contain relevant data.
	 * @param filepath
	 *            A path to a file that will be created with the data written in
	 *            JSON format.
	 * @throws JsonGenerationException
	 * @throws JsonMappingException
	 * @throws IOException
	 */
	public static void writeToFile(Map<String, Object> results, String filepath) 
			throws JsonGenerationException, JsonMappingException, IOException {
		ObjectMapper mapper = new ObjectMapper();
		mapper.writeValue(new File(filepath), results);
	}
	
	/**
	 * Read a JSON file into a Map object.
	 * 
	 * @param filepath
	 *            filepath of the JSON file to read
	 * @throws JsonParseException
	 * @throws JsonMappingException
	 * @throws IOException
	 */
	public static Map<String, Object> readFile(String filepath) 
			throws JsonParseException, JsonMappingException, IOException {
		String json = Files.readFile(filepath);
		ObjectMapper mapper = new ObjectMapper();
		try {
			Map<String, Object> results = mapper.readValue(json, new TypeReference<HashMap<String,Object>>(){});
			return results;
		} catch (JsonParseException e) {
			throw new IOException("Could not parse input file: input file is not in JSON format!");
		} catch (JsonMappingException e) {
			throw new IOException("Could not parse input file: JSON format is invalid!");
		}
	}
	
	public static Map<String, Object> generateConfigurationJson(Prism prism) {
		DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		PrismConfig pc = prism.config();
		Map<String, Object> config = new HashMap<String, Object>();
		config.put("max_substrates", pc.display);
		config.put("homology_cutoff", pc.homologyCutoff);
		config.put("scaffold_limit", pc.scaffoldLimit);
		config.put("score", pc.score);
		config.put("tanimoto_cutoff", pc.tanimotoCutoff);
		config.put("cluster_window", pc.window);
		if (pc.date != null)
			config.put("date", dateFormat.format(pc.date));
		config.put("version", pc.version);
		config.put("aminoglycoside", pc.aminoglycoside);
		config.put("ribosomal", pc.ribosomal);
		config.put("resistance", pc.resistance);
		config.put("regulation", pc.regulation);
		return config;
	}
	
	public static Map<String, Object> generateInputJson(Prism prism) {
		Map<String, Object> input = new HashMap<String, Object>();
		
		Organism organism = prism.genome().organism();
		input.put("description", organism.rawfastaheader());
		input.put("genus", organism.genus());
		input.put("species", organism.species());
		input.put("strain", organism.strain());
		input.put("source", organism.source());
		input.put("accession", organism.accession());
		input.put("classification", organism.classification());
		input.put("filename", prism.genome().filename());
		
		StringBuilder seq = new StringBuilder();
		for (Contig f : prism.genome().contigs()) {
			seq.append(f.sequence());
		}
		
		input.put("length", seq.length());
		return input;
	}
	
	public static List<Object> generate16sJson(Prism prism) {
		List<Object> all16s = new ArrayList<Object>();
		for (RnaSequence r : prism.genome().ribosomalSequences()) {
			Map<String, String> sixteenS = new HashMap<String, String>();
			sixteenS = r.getData();
			sixteenS.put("end", String.valueOf(r.getEnd()));
			sixteenS.put("start", String.valueOf(r.getStart()));
			sixteenS.put("frame", String.valueOf(r.getFrame()));
			all16s.add(sixteenS);
		}
		return all16s;
	}
	
	public static List<Object> generateRegulatoryJson(Contig contig) {
	List<Object> globalDomains = new ArrayList<Object>();
		for (Orf orf : contig.orfs()) {
			for (Domain domain : orf.domains(DomainFamilies.REGULATOR)) {
				String[] domainData = Utils.getDomainData(domain);
				Map<String, Object> dom = new HashMap<String, Object>();
				dom.put("name", domainData[0]);
				dom.put("start", Integer.valueOf(domainData[1]));
				dom.put("stop", Integer.valueOf(domainData[2]));
				dom.put("score", Double.valueOf(domainData[3]));
				dom.put("family", domain.family());
				globalDomains.add(dom);
			}
		}
	return globalDomains;
	}
	
	/**
	 * For a given cluster will generate a Map object representing that cluster. This map can be easily outputed as JSON.
	 * @param prism A Prism object containing the cluster to be used.
	 * @param cluster The cluster to generate the object for.
	 *
	 */
	public static Map<String, Object> generateClusterJson(Prism prism, Cluster cluster) {
		Map<String, Object> clusterData = new HashMap<String, Object>();
		clusterData.put("frame", cluster.frame());
		List<String> molecules = cluster.library().scaffolds();
		if (molecules.size() == 0) 
			molecules = null;
		clusterData.put("molecules", molecules);
		clusterData.put("start", cluster.start());
		clusterData.put("end", cluster.end());
		
		List<String> families = new ArrayList<String>();
		for (ClusterFamilies family : cluster.families())
			families.add(family.toString());
		clusterData.put("family", families);
		
		List<String> types = new ArrayList<String>();
		for (ClusterType type : cluster.types())
			types.add(type.toString());
		clusterData.put("type", types);

		List<Object> sugars = new ArrayList<Object>();
		for (List<Sugar> sugar : cluster.sugars()) {
			List<Object> combo = new ArrayList<Object>();
			if (sugar != null) {
				for (Sugar s : sugar) {
					Map<String, String> subarData = new HashMap<String, String>();
					subarData.put("name", s.name());
					subarData.put("type", s.family().toString());
					subarData.put("smiles", s.smiles());
					combo.add(subarData);
				}
				sugars.add(combo);
			}
		}

		if (sugars.size() > 0) {
			clusterData.put("sugars_combinations", sugars);
		} else {
			clusterData.put("sugars_combinations", null);
		}

		List<Object> orfs = new ArrayList<Object>();

		for (Orf orf : cluster.orfs()) {
			
			//only add orf if it has domain data
			if (orf.domains().size() > 0) {
				Map<String, Object> orfData = new HashMap<String, Object>();
				orfData.put("start", orf.start());
				orfData.put("stop", orf.end());
				orfData.put("frame", orf.frame());
				orfData.put("name", orf.name());
				if (orf.type() != null)
					orfData.put("type", orf.type().toString());
				if (prism.config().saveSequences)
					orfData.put("sequence", orf.sequence());

				List<Object> domains = new ArrayList<Object>();
				// int d = 1;

				for (Domain dm : orf.domains()) {
					String[] domainData = Utils.getDomainData(dm);
					Map<String, Object> domain = new HashMap<String, Object>();
					domain.put("name", domainData[0]);
					domain.put("start", Integer.valueOf(domainData[1]));
					domain.put("stop", Integer.valueOf(domainData[2]));
					domain.put("score", Double.valueOf(domainData[3]));
					domain.put("full_name", domainData[5]);
					domain.put("ks_starter:", Boolean.valueOf(domainData[4]));
					domain.put("modifying", dm.family() == DomainFamilies.TAILORING);
					domain.put("family", dm.family());
					
					boolean onAModule = false;
					int m = 1;
					for (Module md : orf.modules()) {
						if (md.contains(dm)) {
							onAModule = true;
							domain.put("on_module", true);
							domain.put("module_type", md.type().name());
							domain.put("module_number", m);
							domain.put("active", md.isActive());
						}
						m++;
					}
					if (onAModule == false) {
						domain.put("on_module", false);
						if (dm.family() != DomainFamilies.TAILORING) {
							domain.put("orphaned", true);
						}
					}

					List<Object> substrates = new ArrayList<Object>();
					for (ArrayList<String[]> data : Utils.getSubstrateData(prism, dm)) {
						Map<String, Object> substrateData = new HashMap<String, Object>();
						for (String[] substrate : data) {
							substrateData.put("name", substrate[0]);
							substrateData.put("score",
									Double.valueOf(substrate[3]));
							substrateData.put("start",
									Integer.valueOf(substrate[1]));
							substrateData.put("stop",
									Integer.valueOf(substrate[2]));
						}

						substrates.add(substrateData);
					}
					if (Utils.getSubstrateData(prism, dm).size() > 0) {
						domain.put("substrates", substrates);
					} else {
						domain.put("substrates", null);
					}

					List<Object> sources = new ArrayList<Object>();
					for (List<String[]> data : Utils.getSourceData(prism, dm)) {
						Map<String, Object> sourceData = new HashMap<String, Object>();
						for (String[] source : data) {
							sourceData.put("name", source[0]);
							sourceData.put("score",
									Double.valueOf(source[1]));
						}
						sources.add(sourceData);
					}

					if (Utils.getSourceData(prism, dm).size() > 0) {
						domain.put("sources", sources);
					} else {
						domain.put("sources", null);
					}
					domains.add(domain);
					// d++;
				}
				
				orfData.put("domains", domains);
				orfs.add(orfData);
				
			}
		}
		clusterData.put("orfs", orfs);
		return clusterData;
	}
}
