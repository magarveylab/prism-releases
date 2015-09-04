package ca.mcmaster.magarveylab.prism.database;

import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.Map;

import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.map.JsonMappingException;

import ca.mcmaster.magarveylab.enums.ClusterFamilies;
import ca.mcmaster.magarveylab.enums.DeoxySugars;
import ca.mcmaster.magarveylab.enums.DomainFamilies;
import ca.mcmaster.magarveylab.enums.Frames;
import ca.mcmaster.magarveylab.enums.HexoseSugars;
import ca.mcmaster.magarveylab.enums.OrfTypes;
import ca.mcmaster.magarveylab.enums.clusters.ClusterType;
import ca.mcmaster.magarveylab.enums.domains.AminoglycosideDomains;
import ca.mcmaster.magarveylab.enums.domains.BetaLactamDomains;
import ca.mcmaster.magarveylab.enums.domains.DeoxySugarDomains;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.PrerequisiteDomains;
import ca.mcmaster.magarveylab.enums.domains.RegulatorDomains;
import ca.mcmaster.magarveylab.enums.domains.ResistanceDomains;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.enums.domains.TypeIIPolyketideDomains;
import ca.mcmaster.magarveylab.enums.substrates.AcylAdenylatingSubstrates;
import ca.mcmaster.magarveylab.enums.substrates.AcyltransferaseSubstrates;
import ca.mcmaster.magarveylab.enums.substrates.AdenylationSubstrates;
import ca.mcmaster.magarveylab.enums.substrates.SubstrateType;
import ca.mcmaster.magarveylab.prism.Prism;
import ca.mcmaster.magarveylab.prism.blast.BlastSearchResult;
import ca.mcmaster.magarveylab.prism.cluster.analysis.type.ClusterTypeAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.LibraryGenerator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Contig;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Genome;
import ca.mcmaster.magarveylab.prism.data.Library;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.data.Organism;
import ca.mcmaster.magarveylab.prism.data.Substrate;
import ca.mcmaster.magarveylab.prism.data.sugar.DeoxySugar;
import ca.mcmaster.magarveylab.prism.data.sugar.HexoseSugar;
import ca.mcmaster.magarveylab.prism.data.sugar.Sugar;
import ca.mcmaster.magarveylab.prism.util.JsonUtils;
import ca.mcmaster.magarveylab.prism.web.PrismConfig;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Read saved output, in JSON format, back into Prism objects and data
 * structures.
 * 
 * @author skinnider
 *
 */
public class JsonInput {
	
	/**
	 * Read a JSON file into Prism objects.<br>
	 * <br>
	 * This method parses the configuration, the genome, the organism, and the
	 * clusters present in the original PRISM search. The session should be
	 * instantiated based on a HTTP request. 
	 * 
	 * @param filepath
	 *            filepath of the JSON file
	 * @param session
	 *            current session
	 * @return a Prism search object based on the saved JSON results
	 * @throws JsonParseException
	 * @throws JsonMappingException
	 * @throws IOException
	 * @throws ParseException
	 */
	@SuppressWarnings("unchecked")
	public static Prism read(String filepath, Session session) 
			throws JsonParseException, JsonMappingException, IOException, ParseException {
		Map<String,Object> json = JsonUtils.readFile(filepath);
		if (json.get("prism_results") == null)
			throw new IOException("Could not parse JSON: JSON does not contain any PRISM results!");
		json = (Map<String, Object>) json.get("prism_results");
		
		PrismConfig config = parseConfig(json);
		Prism prism = new Prism(config, session);
		
		Genome genome = parseGenome(json, filepath);
		prism.setGenome(genome);
		
		Organism organism = parseOrganism(json);
		genome.setOrganism(organism);

		List<Cluster> contigClusters = genome.contigs().get(0).clusters();
		if (json.get("clusters") != null) {
			System.out.println("[JsonInput] Reading clusters from genomic JSON input");
			List<Cluster> clusters = parseAllClusters(json, config);
			contigClusters.addAll(clusters);
		} else {
			System.out.println("[JsonInput] No clusters found");
		}
		
		for (Cluster cluster : genome.clusters())
			LibraryGenerator.write(cluster.library(), cluster, session);
	
		System.out.println("[JsonInput] Read " + genome.clusters().size() + " clusters from JSON input");
		return prism;
	}

	/**
	 * Parse the configuration of a saved PRISM search result.<br>
	 * <br>
	 * This method attempts to read the display limit, homology and Tanimoto
	 * cutoffs, scaffold limit, clustering window, PRISM version, the date of
	 * the original search, and whether or not bio/cheminformatic dereplication
	 * was enabled.
	 * 
	 * @param json
	 *            map read in from the saved JSON file
	 * @return the configuration of the saved result, if available
	 * @throws ParseException
	 */
	@SuppressWarnings("unchecked")
	public static PrismConfig parseConfig(Map<String,Object> json) throws ParseException {
		PrismConfig config = new PrismConfig();
		config.web = true;
		config.saveSequences = true;
		
		if (json.get("configuration") != null) {
			json = (Map<String, Object>) json.get("configuration");
		} else {
			return config;
		}
		
		DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");

		Integer display = (Integer) json.get("max_substrates");
		if (display != null)
			config.display = display;
		Double homologyCutoff = (Double) json.get("homology_cutoff");
		if (homologyCutoff != null)
			config.homologyCutoff = homologyCutoff;
		Double tanimotoCutoff = (Double) json.get("tanimoto_cutoff");
		if (tanimotoCutoff != null)
			config.tanimotoCutoff = tanimotoCutoff;
		Integer scaffoldLimit = (Integer) json.get("scaffold_limit");
		if (scaffoldLimit != null)
			config.scaffoldLimit = scaffoldLimit;
		Boolean score = (Boolean) json.get("score");
		if (score != null)
			config.score = score;
		Integer window = (Integer) json.get("window");
		if (window != null)
			config.window = window;
		String version = (String) json.get("version");
		if (version != null) {
			config.version = version;
		} else {
			config.version = "";
		}
		String date = (String) json.get("date");
		if (date != null) {
			config.date = dateFormat.parse(date);
		} else {
			Date d = new Date();
			config.date = d;
		}
		Boolean aminoglycoside = (Boolean) json.get("aminoglycoside");
		if (aminoglycoside != null)
			config.aminoglycoside = aminoglycoside; 
		Boolean ribosomal = (Boolean) json.get("ribosomal");
		if (ribosomal != null)
			config.ribosomal = ribosomal; 
		Boolean resistance = (Boolean) json.get("resistance");
		if (resistance != null)
			config.resistance = resistance; 
		Boolean regulation = (Boolean) json.get("regulation");
		if (regulation != null)
			config.regulation = regulation; 
		
		return config;
	}
	
	/**
	 * Parse organism data from a saved PRISM search result. <br>
	 * <br>
	 * This method attempts to read the raw FASTA header, genus, species,
	 * strain, accession number, and classification(s) of the organism, if
	 * available.
	 * 
	 * @param json
	 *            map read in from the saved JSON file
	 * @return the organism data package
	 */
	@SuppressWarnings("unchecked")
	public static Organism parseOrganism(Map<String,Object> json) {
		Organism organism = new Organism();
		if (json.get("organism") != null) {
			json = (Map<String, Object>) json.get("organism");
		} else {
			return organism;
		}

		String description = (String) json.get("description");
		if (description != null)
			organism.setRawFastaHeader(description);
		String genus = (String) json.get("genus");
		if (genus != null)
			organism.setGenus(genus);
		String species = (String) json.get("species");
		if (species != null)
			organism.setSpecies(species);
		String strain = (String) json.get("strain");
		if (strain != null)
			organism.setStrain(strain);
		String source = (String) json.get("source");
		if (source != null)
			organism.setSource(source);
		String accession = (String) json.get("accession");
		if (accession != null)
			organism.setAccession(accession);
		List<String> classification = (List<String>) json.get("classification");
		if (classification != null)
			organism.setClassification(classification);
		
		return organism;
	}
	
	/**
	 * Parse genome data from a saved PRISM search result.<br>
	 * <br>
	 * This method attempts to read the filename of the original sequence.
	 * 
	 * @param json
	 *            map read in from the saved JSON file
	 * @return the genome object
	 */
	public static Genome parseGenome(Map<String,Object> json, String filepath) {
		Genome genome = null;
		String filename = (String) json.get("filename");
		if (filename != null) {
			File file = new File(filename);
			genome = new Genome(file);
		}
		if (genome == null) {
			File file = new File(filepath);
			genome = new Genome(file);
		}
		Contig artificial = new Contig("Artificial contig", "");
		genome.contigs().add(artificial);
		return genome;
	}
	
	/**
	 * Parse a list of clusters from a genome-wide saved JSON result (i.e., the
	 * complete JSON package for an entire PRISM search, as opposed to just the
	 * JSON for a single biosynthetic gene cluster).
	 * 
	 * @param json
	 *            map read in from the saved JSON file
	 * @return all clusters identified in the original PRISM search
	 * @throws ParseException 
	 */
	@SuppressWarnings("unchecked")
	public static List<Cluster> parseAllClusters(Map<String,Object> json, PrismConfig config) 
			throws ParseException {
		List<Cluster> clusters = new ArrayList<Cluster>();
		int i = 1;
		
		if (json.get("clusters") != null) {
			List<Map<String,Object>> maps = (List<Map<String, Object>>) json.get("clusters");
			for (Map<String,Object> map : maps) {
				Cluster cluster = parseCluster(map, config);
				clusters.add(cluster);
				cluster.setIndex(i);
				i++;
			}
		}
		return clusters;
	}
	
	/**
	 * Parse a single cluster, either from a map within a genome-wide saved JSON
	 * result, or from a JSON file representing only a single cluster.
	 * 
	 * Currently set to parse family value as a string for version 1.2.4 and parses it as a list of family values for new versions.
	 * @param json
	 *            map read in from the saved JSON file
	 * @return a single cluster
	 * @throws ParseException
	 */
	@SuppressWarnings("unchecked")
	public static Cluster parseCluster(Map<String,Object> json, PrismConfig config) 
			throws ParseException {
		Cluster cluster = new Cluster();
		
		// parse family
		List<String> families = new ArrayList<String>();
		if (config.version.startsWith("1.2.4")) {
			families.add((String) json.get("family"));
		}else {
			families = (List<String>) json.get("family");
		}
		if (families != null)
			for (String s : families)
				for (ClusterFamilies family : ClusterFamilies.values())
					if (family.toString().equals(s))
						cluster.addFamily(family);
		
		// parse type
		List<String> types = (List<String>) json.get("type");
		if (types != null)
			for (String t : types)
				for (ClusterType type : ClusterTypeAnalyzer.getAllClusterTypes())
					if (type.toString().equals(t))
						cluster.addType(type);
		
		// parse frame 
		String frame = (String) json.get("frame");
		if (frame != null)
			for (Frames clusterFrame : Frames.values())
				if (clusterFrame.toString().equals(frame))
					cluster.setFrame(clusterFrame);
		
		// parse scaffold library 
		if (json.get("molecules") != null) {
			List<String> scaffolds = (List<String>) json.get("molecules");
			Library library = new Library();
			library.scaffolds().addAll(scaffolds);
			cluster.setLibrary(library);
		}
		
		// parse sugars
		List<List<Sugar>> sugars = parseSugars(json);
		cluster.setSugars(sugars);
		
		// parse orfs
		List<Orf> orfs = parseOrfs(json);
		cluster.orfs().addAll(orfs);
		
		return cluster;
	}
	
	/**
	 * Parse a list of sugar combinations for a single cluster from a saved JSON
	 * result.
	 * 
	 * @param json
	 *            map read in from the saved JSON file
	 * @return a cluster's sugar combinations saved in the JSON file
	 */
	@SuppressWarnings("unchecked")
	public static List<List<Sugar>> parseSugars(Map<String,Object> json) {
		List<List<Sugar>> combinations = new ArrayList<List<Sugar>>();
		if (json.get("sugar") != null) {
			List<List<Map<String,String>>> sugarCombinations = (List<List<Map<String, String>>>) json.get("sugar");
			for (List<Map<String,String>> sugarCombination : sugarCombinations) {
				List<Sugar> sugars = new ArrayList<Sugar>();
				for (Map<String,String> map : sugarCombination) {
					Sugar sugar = null;
					String family = (String) map.get("type");
					String name = (String) map.get("name");
					
					if (family.equals("DEOXY")) {
						DeoxySugars type = null;
						for (DeoxySugars ds : DeoxySugars.values())
							if (ds.toString().equals(name))
								type = ds;
						sugar = new DeoxySugar(type);
					} else if (family.equals("HEXOSE")) {
						HexoseSugars type = null;
						for (HexoseSugars hs : HexoseSugars.values())
							if (hs.toString().equals(name))
								type = hs;
						sugar = new HexoseSugar(type);
					}
					
					if (sugar != null)
						sugars.add(sugar);
				}
				combinations.add(sugars);
			}
		}
		System.out.println("[JsonInput] Read " + combinations.size() + " sugar combinations from JSON input");
		return combinations;
	}

	/**
	 * Parse the open reading frames in a cluster from a saved JSON file.
	 * 
	 * @param json
	 *            map read in from the saved JSON file
	 * @return the cluster open reading frames in the saved JSON file.
	 */
	@SuppressWarnings("unchecked")
	public static List<Orf> parseOrfs(Map<String,Object> json) {
		List<Orf> orfs = new ArrayList<Orf>();
		List<Map<String,Object>> orfsMap = (List<Map<String, Object>>) json.get("orfs");
		for (Map<String,Object> orfMap : orfsMap) {
			Integer start = (Integer) orfMap.get("start");
			Integer end = (Integer) orfMap.get("stop");
			String name = (String) orfMap.get("name");
			String frame = (String) orfMap.get("frame");
			
			String sequence = null;
			if (orfMap.get("sequence") != null) {
				sequence = (String) orfMap.get("sequence");
			} else {
				sequence = "";
			}
			
			Orf orf = new Orf(name, sequence); 
			orf.setEnd(end);
			orf.setStart(start);
			orf.setFrame(frame);
			
			// parse type
			String type = (String) orfMap.get("type");
			for (OrfTypes t : OrfTypes.values())
				if (t.toString().equals(type))
					orf.setType(t);
			
			List<Map<String,Object>> domainMaps = (List<Map<String, Object>>) orfMap.get("domains");
			for (Map<String,Object> domainMap : domainMaps) {
				Domain domain = parseDomain(domainMap);
				orf.add(domain);
			}
			orfs.add(orf);
		}
		System.out.println("[JsonInput] Read " + orfs.size() + " orfs from JSON input");
		return orfs;
	}

	/**
	 * Read a biosynthetic domain from a saved JSON file.
	 * 
	 * @param json
	 *            map read in from the saved JSON file
	 * @return the biosynthetic domain in the saved JSON file
	 */
	@SuppressWarnings("unchecked")
	public static Domain parseDomain(Map<String,Object> json) {
		String domainName = (String) json.get("name");
		Integer domainStart = (Integer) json.get("start");
		Integer domainStop = (Integer) json.get("stop");
		Double domainScore = (Double) json.get("score");

		Domain domain = new Domain(domainStart, domainStop, domainScore, domainName);
		
		// get family
		String domainFamily = (String) json.get("family");
		DomainFamilies family = parseDomainFamily(domainFamily);
		domain.setFamily(family);
		
		// get domain type 
		String domainType = (String) json.get("full_name");
		DomainType type = parseDomainType(domainType, domain.family());
		domain.setType(type);
		
		// get substrates
		if (json.get("substrates") != null) {
			List<Map<String,Object>> map = (List<Map<String, Object>>) json.get("substrates");
			List<Substrate> substrates = parseSubstrates(map);
			domain.setSubstrates(substrates);
		}
		
		// get blast 
		if (json.get("sources") != null) {
			List<Map<String,Object>> map = (List<Map<String, Object>>) json.get("sources");
			List<BlastSearchResult> results = parseBlastResults(map);
			domain.blastResults().addAll(results);
		}
		
		return domain;
	}
	
	/**
	 * Parse the family of a biosynthetic domain from a saved JSON result.
	 * 
	 * @param domainFamily
	 *            the raw string associated with the 'family' key in the map
	 *            read in from the saved JSON file
	 * @return the domain family
	 */
	public static DomainFamilies parseDomainFamily(String domainFamily) {
		DomainFamilies f = null;
		for (DomainFamilies family : DomainFamilies.values())
			if (family.toString().equals(domainFamily))
				f = family;
		if (f == null)
			f = DomainFamilies.NULL;
		System.out.println("Parsed domain from family " + f);
		return f;
	}

	/**
	 * Parse the type of a biosynthetic domain from a saved JSON result.
	 * 
	 * @param domainType
	 *            the raw string associated with the 'full_name' key in the map
	 *            read in from the saved JSON file
	 * @return the domain type
	 */
	public static DomainType parseDomainType(String domainType, DomainFamilies family) {
		List<DomainType> types = new ArrayList<DomainType>();
		if (family == DomainFamilies.AMINOGLYCOSIDE) {
			types.addAll(Arrays.asList(AminoglycosideDomains.values()));
		} else if (family == DomainFamilies.BETA_LACTAM) {
			types.addAll(Arrays.asList(BetaLactamDomains.values()));
		} else if (family == DomainFamilies.PREREQUISITE) {
			types.addAll(Arrays.asList(PrerequisiteDomains.values()));
		} else if (family == DomainFamilies.REGULATOR) {
			types.addAll(Arrays.asList(RegulatorDomains.values()));
		} else if (family == DomainFamilies.RESISTANCE) {
			types.addAll(Arrays.asList(ResistanceDomains.values()));
		} else if (family == DomainFamilies.RIBOSOMAL) {
			types.addAll(Arrays.asList(RibosomalDomains.values()));
		} else if (family == DomainFamilies.SUGAR) {
			types.addAll(Arrays.asList(DeoxySugarDomains.values()));
		} else if (family == DomainFamilies.TAILORING) {
			types.addAll(Arrays.asList(TailoringDomains.values()));
		} else if (family == DomainFamilies.THIOTEMPLATED) {
			types.addAll(Arrays.asList(ThiotemplatedDomains.values()));
		} else if (family == DomainFamilies.TYPE_II_POLYKETIDE) {
			types.addAll(Arrays.asList(TypeIIPolyketideDomains.values()));
		}
		
		DomainType t = null;
		for (DomainType type : types)
			if (type.toString().equals(domainType))
				t = type;
		if (t == null)
			t = ThiotemplatedDomains.NULL;
		System.out.println("Read " + t + " domain from JSON input");
		return t;
	}

	/**
	 * Parse the substrates associated with a biosynthetic domain from a saved
	 * JSON result.
	 * 
	 * @param json
	 *            map read in from the saved JSON file
	 * @return list of substrates associated with the domain
	 */
	public static List<Substrate> parseSubstrates(List<Map<String,Object>> json) {
		List<Substrate> substrates = new ArrayList<Substrate>();
		for (Map<String,Object> map : json) {
			String name = (String) map.get("name");
			Double score = (Double) map.get("score");
			Integer start = (Integer) map.get("start");
			Integer end = (Integer) map.get("stop");
			
			SubstrateType type = parseSubstrateType(name);
			Substrate substrate = new Substrate(start, end, score, type);
			substrates.add(substrate);
		}
		System.out.println("[JsonInput] Read " + substrates.size() + " substrates from JSON input");
		return substrates;
	}

	/**
	 * Parse the substrate type of a substrate from a saved JSON result
	 * 
	 * @param name
	 *            the raw string associated with the 'name' key in map read in
	 *            from the saved JSON file
	 * @return the substrate type
	 */
	public static SubstrateType parseSubstrateType(String name) {
		List<SubstrateType> types = new ArrayList<SubstrateType>();
		types.addAll(Arrays.asList(AdenylationSubstrates.values()));
		types.addAll(Arrays.asList(AcylAdenylatingSubstrates.values()));
		types.addAll(Arrays.asList(AcyltransferaseSubstrates.values()));
		SubstrateType t = null;
		for (SubstrateType type : types)
			if (name.equals(type.abbreviation()))
				t = type;
		System.out.println("[JsonInput] Read " + t.toString() + " substrate from JSON input");
		return t;
	}
	
	/**
	 * Parse all saved BLAST search results associated with a biosynthetic
	 * domain from a saved JSON result.<br>
	 * <br>
	 * Because at this point only the name and score of the BLAST search result
	 * are saved in the JSON, all other values in the BLAST search result data
	 * package will be initialized to either zero or an empty string.
	 * 
	 * @param json
	 *            map read in from the saved JSON file
	 * @return all BLAST search results for a given domain
	 */
	public static List<BlastSearchResult> parseBlastResults(List<Map<String,Object>> json) {
		List<BlastSearchResult> results = new ArrayList<BlastSearchResult>();
		for (Map<String,Object> map : json) {
			String name = (String) map.get("name");
			Double score = (Double) map.get("score");
			BlastSearchResult result = new BlastSearchResult("", name, 0, score, 0, 0, 0, 0);
			results.add(result);
		}
		return results;
	}
		
}
