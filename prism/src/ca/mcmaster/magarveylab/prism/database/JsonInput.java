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
import ca.mcmaster.magarveylab.enums.domains.DeoxySugarDomains;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.PrerequisiteDomains;
import ca.mcmaster.magarveylab.enums.domains.RegulatorDomains;
import ca.mcmaster.magarveylab.enums.domains.ResistanceDomains;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.enums.domains.TypeIIPolyketideDomains;
import ca.mcmaster.magarveylab.enums.interfaces.SubstrateType;
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
import ca.mcmaster.magarveylab.prism.database.data.SmallMolecule;
import ca.mcmaster.magarveylab.prism.enums.hmms.AcylAdenylatingHmms;
import ca.mcmaster.magarveylab.prism.enums.hmms.AcyltransferaseHmms;
import ca.mcmaster.magarveylab.prism.enums.hmms.AdenylationHmms;
import ca.mcmaster.magarveylab.prism.homology.data.HomologousCluster;
import ca.mcmaster.magarveylab.prism.tanimoto.data.TanimotoScore;
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

		if (json.get("clusters") != null) {
			getAllContigCluster(genome, json, config);
			/*
			 * List<Cluster> clusters = parseAllClusters(json, config);
			 * contigClusters.addAll(clusters);
			 */
		} else {
			// List<Cluster> contigClusters = genome.contigs().get(0).clusters();
		}
		
		for (Cluster cluster : genome.clusters())
			LibraryGenerator.write(cluster.library(), cluster, session);
	
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
		Boolean thiotemplated = (Boolean) json.get("thiotemplated");
		if (thiotemplated != null)
			config.thiotemplated = thiotemplated; 
		Boolean sugar = (Boolean) json.get("sugar");
		if (sugar != null)
			config.sugar = sugar; 
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
	@SuppressWarnings("unchecked")
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
		
		List<Integer> contigs = (List<Integer>) json.get("contigs");
		for( int i = 0; i<contigs.size(); i++){
			Contig art = new Contig("Artificial contig","");
			art.setLength(contigs.get(i));
			genome.contigs().add(art);
		}

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

	@SuppressWarnings("unchecked")
	public static void getAllContigCluster(Genome genome, Map<String, Object> clusterJson, PrismConfig config) throws ParseException{
		if(clusterJson.get("clusters") != null) {
			List<Map<String,Object>> maps = (List<Map<String, Object>>) clusterJson.get("clusters");
			int i = 1;
			for(Map<String, Object> map : maps) {
				Cluster cluster = parseCluster(map, config);
				cluster.setIndex(i);
				Integer clusterIndex = getClusterIndex(map);
				genome.contigs().get(clusterIndex).clusters().add(cluster);
				i++;
			}
		}
	}
	
	public static Integer getClusterIndex(Map<String,Object> json) throws ParseException{
		return (Integer) json.get("contig");
	}
	
	
	
	/**
	 * Parse a single cluster, either from a map within a genome-wide saved JSON
	 * result, or from a JSON file representing only a single cluster.
	 * 
	 * Parses biosynthetic family as a single string (for version 1.2.4) and
	 * as a list of family values for newer versions.
	 * 
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
		// parse type
		List<String> families = new ArrayList<String>();
		List<String> types = new ArrayList<String>();
		if (config.version.startsWith("1.2.4")) {
			families.add((String) json.get("family"));
			types.add((String) json.get("type"));
		}else {
			
			try{
				families = (List<String>) json.get("family");
				types = (List<String>) json.get("type");
			}
			catch(ClassCastException e){
				String family = (String) json.get("family");
				families.add(family);
				
				String type = (String) json.get("type");
				types.add(type);
			}
			
			
		}
		
		if (families != null)
			for (String s : families)
				for (ClusterFamilies family : ClusterFamilies.values())
					if (family.toString().equals(s))
						cluster.addFamily(family);
		
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
		
		//parse homologous clusters
		cluster.setHomologs(parseHomologousClusters(json));
		
		//parse homologous molecules
		cluster.setScores(parseHomologousMolecules(json));

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
		return combinations;
	}

	/**
	 * Parse the homologous cluster data for a cluster
	 * 
	 * @param json
	 *            map read in from the saved JSON file
	 * @return the list of homologous cluster data
	 */
	@SuppressWarnings("unchecked")
	public static List<HomologousCluster> parseHomologousClusters(Map<String, Object> json) {
		List<HomologousCluster> homologsReturn = new ArrayList<HomologousCluster>();
		if (json.get("homolog_clusters") != null) {
			List<Map<String,Object>> homologs = (List<Map<String, Object>>) json.get("homolog_clusters");
			for (Map<String,Object> homolog : homologs) {
				HomologousCluster h = new HomologousCluster((String) homolog.get("name"));
				h.setCoverage((double) homolog.get("coverage"));
				h.setDomainScore((double) homolog.get("domain_score"));
				h.setIdentityScore((double) homolog.get("identity_score"));
				//TODO need to add weighted domain score to object
				//homolog.get("weighted_domain_score");

				homologsReturn.add(h);
			}
		}
		return homologsReturn;
	}

	/**
	 * Parse the homologous molecule data for a cluster
	 * 
	 * @param json
	 *            map read in from the saved JSON file
	 * @return the list of homologous molecules for the cluster
	 */
	@SuppressWarnings("unchecked")
	public static List<TanimotoScore> parseHomologousMolecules(Map<String,Object> json) {
		List<TanimotoScore> scores = new ArrayList<TanimotoScore>();

		if (json.get("homolog_molecules") != null) {
			List<Map<String,Object>> molecules = (List<Map<String, Object>>) json.get("homolog_molecules");
			
			for (Map<String,Object> mol : molecules) {				
				SmallMolecule query = new SmallMolecule();
				SmallMolecule target = new SmallMolecule();

				query.setName((String) mol.get("query_name"));
				target.setSmiles((String) mol.get("smiles"));
				target.setName((String) mol.get("name"));

				TanimotoScore ts = new TanimotoScore(query, target);
				ts.addScore("ecfp6", ((Number)mol.get("ecfp6")).floatValue());
				ts.addScore("fcfp6", ((Number)mol.get("fcfp6")).floatValue());
				
				scores.add(ts);
			}
						
		}
		return scores;
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
		
		
		Double domainScore = null;
		try{
			domainScore = Double.valueOf(json.get("score").toString());
		}
		catch(Exception e){
			System.out.println(json.get("score").toString());
			System.out.println(json.getClass().getName());
			e.printStackTrace();
			System.exit(0);
		}
			
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
		if (family == DomainFamilies.PREREQUISITE) {
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
			Double score = Double.valueOf(map.get("score").toString());
			Integer start = (Integer) map.get("start");
			Integer end = (Integer) map.get("stop");
			
			SubstrateType type = parseSubstrateType(name);
			Substrate substrate = new Substrate(start, end, score, type);
			substrates.add(substrate);
		}
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
		types.addAll(Arrays.asList(AdenylationHmms.values()));
		types.addAll(Arrays.asList(AcylAdenylatingHmms.values()));
		types.addAll(Arrays.asList(AcyltransferaseHmms.values()));
		SubstrateType t = null;
		for (SubstrateType type : types)
			if (name.equals(type.abbreviation()))
				t = type;
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
			Double score = Double.valueOf((map.get("score").toString()));
			BlastSearchResult result = new BlastSearchResult("", name, 0, score, 0, 0, 0, 0);
			results.add(result);
		}
		return results;
	}
		
}