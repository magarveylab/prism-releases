package ca.mcmaster.magarveylab.prism;

import java.io.File;
import java.io.UnsupportedEncodingException;
import java.net.URLDecoder;
import java.text.SimpleDateFormat;
import java.util.Date;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import ca.mcmaster.magarveylab.prism.web.PrismConfig;
import ca.mcmaster.magarveylab.prism.web.html.PrismReport;
import ca.mcmaster.magarveylab.wasp.exception.ExceptionHandler;
import ca.mcmaster.magarveylab.wasp.session.BasicSession;
import ca.mcmaster.magarveylab.wasp.session.Session;
import ca.mcmaster.magarveylab.wasp.session.SessionListener;
import ca.mcmaster.magarveylab.wasp.session.SessionManager;
import ca.mcmaster.magarveylab.wasp.util.TimeUtil;

/**
 * Execute a PRISM search from the command line. 
 * @author skinnider 
 *
 */
public class PrismDesktop {
	
	private Prism prism;
	
	/**
	 * Create a new command-line PRISM search.
	 * @param genome	the nucleotide sequence file to run on 
	 */
	public PrismDesktop(String genome) {
		try {
			PrismConfig config = new PrismConfig();
			Session session = createSession(config);			
			config.input = genome;
			
			prism = new Prism(config, session);
			session.setWebapp(prism);
		} catch (UnsupportedEncodingException e) {
			handleException(e, "Error: unsupported encoding!");
		}
	}
	
	/**
	 * Get the PRISM search created from the command line.
	 * @return	PRISM search 
	 */
	public Prism search() {
		return prism;
	}
	
	/**
	 * Execute a Prism search from a jar file.
	 * Usage: java -jar prism.jar -file sequence.fa [options]
	 * @param args	command line arguments
	 */
	public static void main(String[] args) {
		try {
			Options help = createHelp();
			Options version = createVersion();
			Options options = createOptions();
			
			CommandLineParser parser = new GnuParser();
			
			// parse help first
			CommandLine line = parser.parse(help, args, true);
			if (line.hasOption("help")) {
				HelpFormatter formatter = new HelpFormatter();
				 String header = "PRISM: PRediction Informatics for Secondary Metabolomes \n" +
						 "Search a genome for novel nonribosomal peptide and polyketide molecules.\n" +
						 "Requires blastp, hmmsearch, and BioPerl as dependencies.\n\n";
				 String footer = "\nGUI available at http://magarveylab.com/prism";
				 formatter.printHelp("prism", header, options, footer, true);
			} else { 
				line = parser.parse(version, args, true);
				if (line.hasOption("version")) {
					System.out.println("PRISM: version " + new PrismConfig().version);
				} else {
					line = parser.parse(options, args);
					PrismConfig config = parseCommandLine(line);
					Session session = createSession(config);

					Prism prism = new Prism(config, session);
					session.setWebapp(prism);
					prism.run();
				}
			} 
		} catch (UnsupportedEncodingException e) {
			handleException(e, "Error: unsupported encoding!");
		} catch (ParseException e) {
			handleException(e, "Error parsing command line arguments!");
		} catch (IllegalArgumentException e) {
			handleException(e, "Error: must specify input sequence file!");
		} catch (Exception e) {
			handleException(e, "Error!");
		} finally {
			System.exit(0);
		}
	}

	/**
	 * Create the help option.
	 * @return	newly created help option 
	 */
	public static Options createHelp() {
		Options options = new Options();
		Option help = new Option("help", "help", false, "Print this message");
		options.addOption(help);
		return options;
	}


	/**
	 * Create the version option.
	 * @return	newly created version option 
	 */
	public static Options createVersion() {
		Options options = new Options();
		Option version = new Option("version", "version", false, "Print the current version and exit");
		options.addOption(version);
		return options;
	}

	/**
	 * Create all other command line options. 
	 * @return	newly created options
	 */
	@SuppressWarnings("static-access")
	public static Options createOptions() {
		Options options = new Options();
		
		// construct options with values
		Option file = OptionBuilder.withLongOpt("file").hasArg().withArgName("FILE")
				.withDescription("Sequence file").create("f");
		Option window = OptionBuilder.withLongOpt("window").withArgName("WINDOW")
				.hasArg().withDescription("Biosynthetic cluster window").create("w");
		Option structureLimit = OptionBuilder.withLongOpt("limit").withArgName("LIMIT")
				.hasArg().withDescription("Scaffold library limit").create("l");
		Option tanimotoCutoff = OptionBuilder.withLongOpt("tanimoto").withArgName("TANIMOTO").hasArg().
				withDescription("Tanimoto cutoff").create("c");
		Option homologyCutoff = OptionBuilder.withLongOpt("homology").withArgName("HOMOLOGY").hasArg().
				withDescription("Homology score cutoff").create("h");
		Option root = OptionBuilder.withLongOpt("root").withArgName("ROOT").hasArg().
				withDescription("Root directory for prism").create("r");
		Option setOutput = OptionBuilder.withLongOpt("outputfolder").withArgName("OUTPUT FOLDER").hasArg().
				withDescription("Set an output folder to output various files to. Ex. JSON").create("out");
		Option geneSearchMode = OptionBuilder.withLongOpt("genesearch").withArgName("GENE PREDICTION MODE")
				.hasArg().withDescription("Method to use for gene prediction (prodigal, genemark, or allorfs)").create("g");
		
		// construct boolean options
		Option organism = new Option("o", "Organism", false, "Parse organism information from GenBank header");
		Option score = new Option("s", "Score", false, "Enable known natural product scoring");
		Option resistance = new Option("resistance", "Resistance", false, "Enable resistance domain search");
		Option regulator = new Option("regulator", "Regulator", false, "Enable regulator domain search");
		Option aminoglycoside = new Option("aminoglycoside", "Aminoglycoside", false, "Enable aminoglycoside domain search");
		Option ribosomal = new Option("ribosomal", "Ribosomal", false, "Enable RiPP domain search");
		Option find16s = new Option("16S", "16S", false, "Detect 16S sequences");
		Option web = new Option("web", "web", false, "Turn on HTML output and graphical output generation");
		Option saveSequences = new Option("saveseq", "saveSequences", false, 
				"Write open reading frame sequences in JSON output");
		Option help = new Option("help", "help", false, "Print this message");
		Option version = new Option("version", "version", false, "Print the current version and exit");
		
		options.addOption(file);
		options.addOption(window);
		options.addOption(structureLimit);
		options.addOption(tanimotoCutoff);
		options.addOption(homologyCutoff);
		options.addOption(score);
		options.addOption(root);
		options.addOption(geneSearchMode);
		options.addOption(organism);
		options.addOption(resistance);
		options.addOption(regulator);
		options.addOption(aminoglycoside);
		options.addOption(ribosomal);
		options.addOption(find16s);
		options.addOption(setOutput);
		options.addOption(web);
		options.addOption(saveSequences);
		options.addOption(help);
		options.addOption(version);

		return options;
	}
	
	/**
	 * Parse the command line input.
	 * 
	 * @param line
	 *            line to parse
	 * @return a PRISM configuration package
	 * @throws ParseException
	 */
	public static PrismConfig parseCommandLine(CommandLine line) throws ParseException {
		PrismConfig config = new PrismConfig();

		// set HMM cutoffs 
		if (line.hasOption("f")) {
			String file = line.getOptionValue("f");
			config.input = file;
			System.out.println("[Prism] Set file to " + config.input);		
		}

		// set other options
		if (line.hasOption("s")) 
			config.score = true;
		if (line.hasOption("w")) {
			String value = line.getOptionValue("w");
			config.window = Integer.parseInt(value);
			System.out.println("[Prism] Set clustering window to " + config.window);
		}
		if (line.hasOption("l")) {
			String value = line.getOptionValue("l");
			config.scaffoldLimit = Integer.parseInt(value);
			System.out.println("[Prism] Set scaffold limit to " + config.scaffoldLimit);
		}
		if (line.hasOption("c")) {
			String value = line.getOptionValue("c");
			config.tanimotoCutoff = Double.parseDouble(value);
			System.out.println("[Prism] Set tanimoto cutoff to " + config.tanimotoCutoff);
		}
		if (line.hasOption("h")) {
			String value = line.getOptionValue("h");
			config.homologyCutoff = Double.parseDouble(value);
			System.out.println("[Prism] Set homology cutoff to " + config.homologyCutoff);
		}
		if (line.hasOption("r")) {
			String value = line.getOptionValue("r");
			config.root = value;
			System.out.println("[Prism] Set the root folder to " + config.root);
		}
		if (line.hasOption("o")) {
			config.parseorganism = true;
			System.out.println("[Prism] Parsing organism data from GenBank-style header");
		}
		if (line.hasOption("16S")) {
			config.find16s = true;
			System.out.println("[Prism] Finding 16S sequences");
		}
		if (line.hasOption("web")) {
			config.web = true;
			System.out.println("[Prism] Set to output html results.");
		}
		if (line.hasOption("saveseq")) {
			config.saveSequences = true;
			System.out.println("[Prism] Set to write open reading frame sequences to JSON results.");
		}
		if (line.hasOption("out")) {
			String value = line.getOptionValue("out");
			config.output = value;
			System.out.println("[Prism] Set output folder to " + value);
		}
		if (line.hasOption("ribosomal")) {
			config.ribosomal = true;
			System.out.println("[Prism] Executing ribosomal domain search");
		}
		if (line.hasOption("resistance")) {
			config.resistance = true;
			System.out.println("[Prism] Executing resistance domain search");
		}
		if (line.hasOption("regulator")) {
			config.regulation = true;
			System.out.println("[Prism] Executing regulator domain search");
		}
		if (line.hasOption("aminoglycoside")) {
			config.regulation = true;
			System.out.println("[Prism] Executing aminoglycoside domain search");
		}
		return config;
	}
	
	/**
	 * Create a new session to use in a command line PRISM search. 
	 * @param config	current PRISM configuration
	 * @return			new session 
	 * @throws UnsupportedEncodingException
	 */
	public static Session createSession(PrismConfig config) throws UnsupportedEncodingException {
		SessionManager sessionManager = SessionManager.getSessionManager();

		// create a new session with random ID
		Date date = config.date;
		SimpleDateFormat dateFormat = new SimpleDateFormat("yyyyMMdd-kkmm");
		String prefix = dateFormat.format(date);
		
		String sessionID = prefix + "-" + (int) Math.floor(Math.random() * 1000000000);
		Session old = sessionManager.getSession(sessionID); 
		if (old != null)
			sessionManager.removeSession(sessionID);
		Session session = new BasicSession();
		session.setID(sessionID);
		
		if (config.root == null) {
			// get jar location to set root, dir
			String path = PrismDesktop.class.getProtectionDomain().getCodeSource().getLocation().getPath();
			String decodedPath = URLDecoder.decode(path, "UTF-8");
			String root = decodedPath.substring(0, decodedPath.lastIndexOf(File.separator)) + File.separator;
			System.out.println("[PrismDesktop] Root was not set; setting root to " + root);
			config.root = root;
		}

		session.setRoot(config.root);
		String dir = config.root + "prism" + File.separator + sessionID + File.separator;
		session.setDir(dir);
		
		// make sure session directory exists
		File sessionDir = new File(dir);
		if (!sessionDir.isDirectory()) 
			sessionDir.mkdirs();
		 
		// set heart beat
		String heartbeat = TimeUtil.getTimeTag();
		session.setLastHeartBeat(heartbeat);

		// set report
		PrismReport report = new PrismReport(session);
		session.setReport(report);
		
		// set listener
		SessionListener listener = new SessionListener();
		session.setListener(listener);

		// set exception handler
		ExceptionHandler handler = new ExceptionHandler(listener);
		session.setExceptionHandler(handler);
		
		// register session
		sessionManager.addSession(sessionID, session);
		
		return session;
	}
	
	/**
	 * Handle an exception by printing the stack trace to the command line, and exiting.  
	 * @param e			exception
	 * @param message	message associated with the exception cause
	 */
	public static void handleException(Exception e, String message) {
		System.out.println("[Prism] " + message);
		e.printStackTrace();
		System.exit(0);
	}

}
