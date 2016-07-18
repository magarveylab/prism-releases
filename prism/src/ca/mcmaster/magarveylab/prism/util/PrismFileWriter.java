package ca.mcmaster.magarveylab.prism.util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Properties;

import org.codehaus.jackson.JsonGenerationException;
import org.codehaus.jackson.map.JsonMappingException;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.io.listener.PropertiesListener;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.smiles.FixBondOrdersTool;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ca.mcmaster.magarveylab.enums.DomainFamilies;
import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.prism.Prism;
import ca.mcmaster.magarveylab.prism.cluster.analysis.DomainAnalyzer;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Contig;
import ca.mcmaster.magarveylab.prism.data.Genome;
import ca.mcmaster.magarveylab.prism.data.Library;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.data.sugar.Sugar;
import ca.mcmaster.magarveylab.prism.database.JsonOutput;
import ca.mcmaster.magarveylab.prism.fasta.FastaWriter;
import ca.mcmaster.magarveylab.prism.web.PrismConfig;
import ca.mcmaster.magarveylab.prism.web.html.graph.CircularContigGraph;
import ca.mcmaster.magarveylab.prism.web.html.graph.CircularGenomeGraph;
import ca.mcmaster.magarveylab.wasp.session.Session;
import ca.mcmaster.magarveylab.prism.orfs.GenePredictionModes;
import ca.mcmaster.magarveylab.prism.util.SmilesIO;

/**
 * Writes files from PRISM searches.
 * 
 * @author skinnider
 *
 */
public class PrismFileWriter {

	/**
	 * Write all files for a PRISM search, including config file, genome graphs,
	 * FASTA files, GRAPE file, and sugar file.
	 * 
	 * @param genome
	 *            genome object from the PRISM search
	 * @param session
	 *            current session
	 * @throws IOException
	 */
	public static void writeAllFiles(Genome genome, Session session) throws IOException {
		Prism prism = (Prism) session.webapp();
		PrismConfig config = prism.config();
		
		if (config.web) {
			// write genome graphs
			writeGenomeGraphs(genome, session);

			for (Contig contig : genome.contigs()) {
				for (Cluster cluster : contig.clusters()) {
					// write FASTA files
					writeBiosyntheticClusterOrfs(cluster, contig, session);
					writeAllClusterOrfs(cluster, contig, session);
					writeGenomicClusterSequence(cluster, contig, session);

					// write GRAPE file
					writeGrapeFile(cluster, session);

					// write sugar file
					writeSugarFile(cluster, session);
					
					//write json file
					//writeClusterJson(cluster, session);
				}
			}
		}
	
		// write config file
		writeConfigFile(session);

		// write JSON file 
		//writeJsonFile(session);
	}
	
	public static void writeAllJson(Genome genome, Session session) throws IOException{
		for (Contig contig : genome.contigs()) {
			for (Cluster cluster : contig.clusters()) {
				//write json file
				writeClusterJson(cluster, session, contig.index());
			}
		}
		// write JSON file
		writeJsonFile(session);
	}
	
	/**
	 * Write the JSON 
	 * @param session	current session
	 * @throws JsonGenerationException
	 * @throws JsonMappingException
	 * @throws IOException
	 */
	public static void writeJsonFile(Session session) throws JsonGenerationException, JsonMappingException, IOException {
		Prism prism = (Prism) session.webapp();
		PrismConfig config = prism.config();
		
		JsonOutput json = new JsonOutput(session);
		json.generateJson();
		
		String folder = config.output != null ? config.output : session.dir();
		json.writeToFile(folder);
		
		System.out.println("[PrismFileWriter] Generated JSON file");
	}
	
	/**
	 * Write a single json object for a single cluster.
	 *
	 * @param cluster
	 *            cluster in question
	 * @param session
	 *            current session
	 * @param contigIndex
	 *            Contig that the cluster is within
	 * @throws IOException
	 */
	public static void writeClusterJson(Cluster cluster, Session session,
			Integer contigIndex) throws IOException {
		String clusterJsonPath = session.dir() + "cluster_" + cluster.index()
				+ ".json";
		Map<String, Object> data = JsonUtils.generateClusterJsonFile(
				(Prism) session.webapp(), cluster, contigIndex);
		JsonUtils.writeToFile(data, clusterJsonPath);
		cluster.addFile("clusterJson", clusterJsonPath);
	}
	
	/**
	 * Create and write circular genome graphs and write as SVG images, if the
	 * input genome has a single FASTA item.
	 * 
	 * @param genome
	 *            genome object from the PRISM search
	 * @param session
	 *            current session
	 */
	public static void writeGenomeGraphs(Genome genome, Session session) {
		// circular genome graph 
		try {
			for (Contig contig : genome.contigs()) {
				for (Cluster cluster : contig.clusters()) {
					CircularContigGraph contigGraph = new CircularContigGraph(
							contig, cluster, 400);
					contigGraph.write(session.dir() + "cluster_contigs_"
							+ cluster.index() + ".svg");
					cluster.setContigGraph(contigGraph);
				}
			}
			
			if (genome.contigs().size() == 1) {
				Contig contig = genome.contigs().get(0);
				if (contig.sequence() == null
						|| contig.sequence().isEmpty()
						|| contig.sequence().equals(""))
					return; // for JSON
				CircularGenomeGraph genomeGraph = new CircularGenomeGraph(genome, 400);
				genomeGraph.write(session.dir() + "genome.svg");
				genome.setGraph(genomeGraph); 

				for (Cluster cluster : genome.clusters()) {
					// create and embed circular genome graphs for each cluster
					CircularGenomeGraph clusterGraph = new CircularGenomeGraph(genome, cluster, 200);
					clusterGraph.write(session.dir() + "cluster_" + cluster.index() + ".svg");
					cluster.setGraph(clusterGraph);
				}
			}
			
		} catch (Exception e) {
			System.out.println("Error: could not write circular genome graphs");
		}
	}
	
	/**
	 * Write all cluster open reading frames which contain a biosynthetic domain
	 * to a multi-FASTA file.
	 * 
	 * @param cluster
	 *            cluster in question
	 * @param contig
	 *            parent contig
	 * @param session
	 *            current session
	 * @throws IOException
	 */
	public static void writeBiosyntheticClusterOrfs(Cluster cluster, Contig contig, Session session) throws IOException {
		String scaffoldOrfsPath = session.dir() + "cluster_" + cluster.index() + ".fasta";
		File fastaFile = new File(scaffoldOrfsPath);
		if (!fastaFile.exists())
			fastaFile.createNewFile();

		List<Orf> orfs = new ArrayList<Orf>();
		for (Orf orf : cluster.orfs())
			if (orf.domains().size() > 0)
				orfs.add(orf);
		
		FastaWriter.printOrfsToFasta(orfs, scaffoldOrfsPath);		
		cluster.addFile("scaffoldOrfs", scaffoldOrfsPath);
	}
	
	/**
	 * Write all of the open reading frames within the boundary of a cluster to a file, regardless of whether
	 * they contain biosynthetic domains.
	 * @param cluster	cluster in question
	 * @param contig		parent contig
	 * @param session	current session 
	 * @throws IOException 
	 */
	public static void writeAllClusterOrfs(Cluster cluster, Contig contig, Session session) throws IOException {
		Prism prism = (Prism) session.webapp();
		PrismConfig config = prism.config();
		int window = config.window;

		List<Orf> orfs = null;
		if (contig.sequence() == null
				|| contig.sequence().isEmpty()
				|| contig.sequence().equals("")) {
			// for JSON: just write all orfs
			orfs = cluster.orfs();
		} else {
			orfs = contig.getAllOrfs(cluster, window);
		}

		String clusterOrfsPath = session.dir() + "cluster_" + cluster.index() + "_full.fasta";
		FastaWriter.printOrfsToFasta(orfs, clusterOrfsPath);
		cluster.addFile("clusterOrfs", clusterOrfsPath);
	}
	
	/**
	 * Set the genomic sequence of the cluster, and print it to a FASTA file.
	 * @param cluster	cluster in question
	 * @param contig	parent contig
	 * @param session	current session
	 * @throws IOException 
	 */
	public static void writeGenomicClusterSequence(Cluster cluster, Contig contig, Session session) throws IOException {
		if (contig.sequence() == null
				|| contig.sequence().isEmpty()
				|| contig.sequence().equals(""))
			return; // for JSON
		
		Prism prism = (Prism) session.webapp();
		PrismConfig config = prism.config();

		Orf first = cluster.first();
		Orf last = cluster.last();
		int start = first.start() - config.window;
		if (start < 0)
			start = 0;
		int end = last.end() + config.window;
		if (end > contig.length())
			end = contig.length();

		String sequence = contig.sequence().substring(start, end);

		String header = "FASTA_item_" + contig.index() + "_Cluster_" + cluster.index();
		String path = session.dir() + "cluster_" + cluster.index() + "_genomic.fasta";
		FastaWriter.printToFasta(header, sequence, path);

		cluster.addFile("sequence", path);
	}

	/**
	 * Generate mol files from a list of clusters.
	 * @param clusters	clusters to write mol files for
	 * @return			the absolute path of each mol file
	 * @throws CDKException 
	 * @throws IOException 
	 */
	public static List<String> writeMolFile(List<Cluster> clusters, Session session) throws CDKException, IOException {
		// force MDLV2000Writer to write 2D coordinates, even when 3D available
		Properties mdlSettings = new Properties();
		mdlSettings.setProperty("ForceWriteAs2DCoordinates","true");
		PropertiesListener propertyListener = new PropertiesListener(mdlSettings);
		MDLV2000Writer writer;
		StructureDiagramGenerator sdg = new StructureDiagramGenerator();
		
		List<String> results = new ArrayList<String>();
		for (Cluster cluster : clusters) {
			Library library = cluster.library();
			String scaffold = library.scaffolds().get(0);
			if (scaffold != null) {
				IAtomContainer molecule = SmilesIO.molecule(scaffold);

				File peptideFile = new File(session.dir() + "cluster_" + cluster.index() + ".mol");
				writer = new MDLV2000Writer(new FileWriter(peptideFile));
				results.add(peptideFile.getAbsolutePath());
				writer.addChemObjectIOListener(propertyListener);
				
				// fix aromaticity
				AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
				FixBondOrdersTool fbot = new FixBondOrdersTool();
				molecule = fbot.kekuliseAromaticRings((IMolecule) molecule);
				
				// use StructureDiagramGenerator to generate coordinates for IMolecule of parsed smiles 
				sdg.setMolecule((IMolecule) molecule);
				// necessary for appropriate rendering, default bond length is too compressed
				sdg.setBondLength( 3 );
				sdg.generateCoordinates();
				molecule = sdg.getMolecule();
				writer.write(molecule);
				writer.close();
			}
		}
		return results;
	}

	/**
	 * Write a scaffold library to a tab-delimited text file.
	 * @param library	library to write
	 * @param path		location of the file
	 * @throws IOException
	 */
	public static void writeLibraryFile(Library library, String path) throws IOException {
		// check file exists
		File file = new File(path);
		if (!file.exists())
			file.createNewFile();
		FileWriter fw = new FileWriter(file);
		BufferedWriter bw = new BufferedWriter(fw);
	
		List<String> scaffolds = library.scaffolds();		
		for (int i = 0; i < scaffolds.size(); i++) {
			String scaffold = scaffolds.get(i);
			if (scaffold != null) {
				String name = "Scaffold_" + (i + 1);
				bw.append(name + "\t" + scaffold + "\n");
			} else {
				System.out.println("[LibraryUtil] Error: Scaffold " + i + " is null!");
			}
		}
		
		bw.close();
	}

	/**
	 * Write a scaffold library to a file formatted for use with iSNAP as a user-input library.
	 * @param library	library to write
	 * @param path		location of the file
	 * @throws IOException
	 */
	public static void writeGNPFile(Cluster cluster, String path) throws IOException {
		// check file exists
		File file = new File(path);
		if (!file.exists())
			file.createNewFile();
		FileWriter fw = new FileWriter(file);
		BufferedWriter bw = new BufferedWriter(fw);
	 
		bw.append("; Cluster " + cluster.index() + "\n");
		
		Library library = cluster.library();
		List<String> scaffolds = library.scaffolds();
		for (int i = 0; i < scaffolds.size(); i++) {
			String scaffold = scaffolds.get(i);
			if (scaffold != null) {
				String name = "Scaffold_" + (i + 1);
				bw.append("#define $" + name + " = " + scaffold + "\n");
			} else {
				System.out.println("[LibraryUtil] Error: Scaffold " + i + " is null!");
			}
		}
		
		bw.close();
	}
	
	/**
	 * Write the settings used in this PRISM search to a configuration file. 
	 * @param session		current session
	 * @throws IOException
	 */
	public static void writeConfigFile(Session session) throws IOException {
		Prism prism = (Prism) session.webapp();
		PrismConfig config = prism.config();
		
		String path = session.dir() + "config.txt";
		File file = new File(path);
		if (!file.exists())
			file.createNewFile();
		FileWriter fw = new FileWriter(file);
		BufferedWriter bw = new BufferedWriter(fw);
		
		bw.append("PRISM configuration" + "\n");
		
		if (config.version != null)
			bw.append("Version: " + config.version + "\n");
		bw.append("Search ID: " + session.id() + "\n");
		if (config.date != null)
			bw.append("Search date: " + config.date.toString() + "\n");
		bw.append("Query sequence: " + prism.genome().filename() + "\n");
		
		bw.append("\n");
		bw.append("Global settings" + "\n");
		
		// get orf prediction
		StringBuffer sb = new StringBuffer();
		for (GenePredictionModes mode : config.genePredictionModes)
			sb.append(mode.toString() + " ");
		bw.append("Open reading frame prediction: " + sb.toString() + "\n");
		
		bw.append("Cluster window: " + config.window + "\n");
		bw.append("Generated structure limit: " + config.scaffoldLimit + "\n");
		bw.append("Result display size: " + config.display + "\n");
		bw.append("Thiotemplated gene detection: " + config.thiotemplated + "\n");
		bw.append("Deoxysugar gene detection: " + config.sugar + "\n");
		bw.append("RiPP gene detection: " + config.ribosomal + "\n");
		bw.append("Resistance gene detection: " + config.resistance + "\n");
		bw.append("Regulatory gene detection: " + config.regulation + "\n");
		
		bw.append("\n");
		bw.append("Similarity search" + "\n");
		bw.append("Detecting knowns: ");
		if (config.score) {
			bw.append("YES" + "\n");
			bw.append("Tanimoto score cutoff: " + config.tanimotoCutoff + "\n");
			bw.append("Homologous cluster cutoff: " + config.homologyCutoff + "\n");
		} else {
			bw.append("NO" + "\n");
		}
		bw.close();

	}

	/**
	 * Write a GRAPE (Genetic Reverse Assembly Prediction Engine) cluster file. 
	 * @param cluster	cluster to write
	 * @param file		location of the file
	 * @throws IOException
	 */
	public static void writeGrapeFile(Cluster cluster, Session session) throws IOException {
		String grapeFilePath = session.dir() + "cluster_" + cluster.index() + "_grape.txt";
		File file = new File(grapeFilePath);
		FileWriter fw = new FileWriter(file);
		BufferedWriter bw = new BufferedWriter(fw);
		
		StringBuffer sb = new StringBuffer();
		PrismStringBuffer psb = new PrismStringBuffer(sb);
	 
		for (Orf orf : cluster.orfs()) {
			if (orf.hasActiveModule()) {
				psb.appendLine(orf.name());
				for (Module module : orf.modules()) {
					psb.appendLine("\t" + "Module: " + module.type());
					for (Domain domain : module.domains()) {
						if (DomainAnalyzer.isCyclization(domain)) {
							psb.appendLine("\t\tCyclization");
						} else {
							psb.appendLine("\t\t" + domain.type());
						}
					}
					if (module.isAdenylationModule()) {
						Domain scaffold = module.scaffold();
						if (scaffold != null && scaffold.topSubstrate() != null) {
							String abbreviation = scaffold.topSubstrate().type().abbreviation();
							psb.appendLine("\t\t\t" + abbreviation);
						}
					} else if (module.type() == ModuleTypes.C_STARTER) {
						for (Domain domain : module.domains())
							psb.appendLine("\t\t\t" + domain.blastResults().get(0).subject());
					}
				}
				psb.append("\n");
			}
		}

		// also, all tailoring domains
		psb.appendLine("\nTailoring");
		for (Domain domain : cluster.domains()) {
			if (domain.family() == DomainFamilies.TAILORING) {
				psb.appendLine(domain.type() + "");
			}
		}
		
		bw.append(psb.toString());
		bw.flush();
		bw.close();
		
		cluster.addFile("grape", grapeFilePath);
	}
	
	public static void writeSugarFile(Cluster cluster, Session session) throws IOException {
		String path = session.dir() + "cluster_" + cluster.index() + "_sugars.txt";
		List<List<Sugar>> sugars = cluster.sugars();
		
		if (sugars.size() == 0)
			return;
		
		File file = new File(path);
		if (!file.exists())
			file.createNewFile();
		FileWriter fw = new FileWriter(file);
		BufferedWriter bw = new BufferedWriter(fw);
		
		for (List<Sugar> combination : sugars) {
			if (combination != null && combination.size() > 0) {
				StringBuffer sb = new StringBuffer();
				for (int i = 0; i < combination.size() - 1; i++)
					sb.append(combination.get(i).toString() + ", ");
				sb.append(combination.get(combination.size() - 1).toString());
				sb.append("\t");
				for (int i = 0; i < combination.size() - 1; i++)
					sb.append(combination.get(i).smiles() + ",");
				sb.append(combination.get(combination.size() - 1).smiles());
				bw.append(sb.toString());
				bw.append("\n");
			}
		}
		
		bw.close();
	}
	
}
