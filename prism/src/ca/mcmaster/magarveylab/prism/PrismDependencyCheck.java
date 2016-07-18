package ca.mcmaster.magarveylab.prism;

import java.io.File;
import java.io.IOException;
import java.util.List;

import ca.mcmaster.magarveylab.enums.version.Version;
import ca.mcmaster.magarveylab.prism.blast.BlastSearchResult;
import ca.mcmaster.magarveylab.prism.blast.BlastpSearch;
import ca.mcmaster.magarveylab.prism.data.Genome;
import ca.mcmaster.magarveylab.prism.fasta.FastaReader;
import ca.mcmaster.magarveylab.prism.genome.HmmSearch;
import ca.mcmaster.magarveylab.prism.genome.data.HmmSearchResult;
import ca.mcmaster.magarveylab.prism.motif.FimoSearch;
import ca.mcmaster.magarveylab.prism.motif.FimoSearchResult;
import ca.mcmaster.magarveylab.prism.orfs.GenePredictionModes;
import ca.mcmaster.magarveylab.prism.orfs.ProdigalSearch;
import ca.mcmaster.magarveylab.prism.tanimoto.Ecfp6Fcfp6Fingerprinter;
import ca.mcmaster.magarveylab.prism.tanimoto.data.Fingerprint;
import ca.mcmaster.magarveylab.prism.util.exception.BadSmilesToFingerprinterException;
import ca.mcmaster.magarveylab.prism.util.exception.DependencyBinaryException;
import ca.mcmaster.magarveylab.prism.util.exception.DependencyException;
import ca.mcmaster.magarveylab.prism.util.exception.ProdigalSearchException;
import ca.mcmaster.magarveylab.prism.util.exception.VersionCompatibilityException;
import ca.mcmaster.magarveylab.prism.util.exception.VersionFormatException;
import ca.mcmaster.magarveylab.prism.web.PrismConfig;
import ca.mcmaster.magarveylab.wasp.session.Session;

/**
 * Ensure dependencies are present before running PRISM.
 * 
 * @author skinnider
 *
 */
public class PrismDependencyCheck {

	private PrismConfig config;
	private Session session;

	public PrismDependencyCheck(PrismConfig config, Session session) {
		this.config = config;
		this.session = session;
	}

	/**
	 * Ensure all PRISM dependencies are present.
	 * 
	 * @throws DependencyException
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws BadSmilesToFingerprinterException
	 */
	public void run() throws DependencyException, IOException,
			InterruptedException, BadSmilesToFingerprinterException {
		checkEnumsVersion();
		checkScripts();
		checkBinaryDependencies();
	}

	/**
	 * Ensure that the version of the enums package matches that of PRISM to
	 * avoid unexplained null pointer exceptions.
	 * 
	 * @throws DependencyException
	 *             if the enums package version number is not compatible with
	 *             the PRISM version number
	 */
	public void checkEnumsVersion() throws DependencyException {
		String enumsVersionRaw = Version.getValue();
		String prismVersionRaw = config.version;

		// get major.minor.patch
		String enumsVersion = enumsVersionRaw.split("-")[0]; // remove SNAPSHOT
		String[] enums = enumsVersion.split("\\.");
		if (enums.length < 3)
			throw new VersionFormatException("Error: could not parse version "
					+ "in format MAJOR.MINOR.PATCH for version string: "
					+ enumsVersionRaw);
		String prismVersion = prismVersionRaw.split("-")[0]; // remove SNAPSHOT
		String[] prism = prismVersion.split("\\.");
		if (prism.length < 3)
			throw new VersionFormatException("Error: could not parse version "
					+ "in format MAJOR.MINOR.PATCH for version string: "
					+ prismVersionRaw);

		int enumsMajor = Integer.parseInt(enums[0]);
		int prismMajor = Integer.parseInt(prism[0]);
		if (enumsMajor != prismMajor)
			throw new VersionCompatibilityException(
					"Error: enums package major version " + enumsMajor
							+ "(enums " + enumsVersionRaw
							+ ") is not compatible with PRISM major version "
							+ prismMajor + "(PRISM " + prismVersionRaw + ")");

		int enumsMinor = Integer.parseInt(enums[1]);
		int prismMinor = Integer.parseInt(prism[1]);
		if (enumsMinor != prismMinor)
			throw new VersionCompatibilityException(
					"Error: enums package minor version " + enumsMinor
							+ "(enums " + enumsVersionRaw
							+ ") is not compatible with PRISM minor version "
							+ prismMinor + "(PRISM " + prismVersionRaw + ")");

		int enumsPatch = Integer.parseInt(enums[2]);
		int prismPatch = Integer.parseInt(prism[2]);
		if (enumsPatch != prismPatch)
			throw new VersionCompatibilityException(
					"Error: enums package patch version " + enumsPatch
							+ "(enums " + enumsVersionRaw
							+ ") is not compatible with PRISM patch version "
							+ prismPatch + "(PRISM " + prismVersionRaw + ")");
	}

	/**
	 * To be able to execute installed binaries (hmmsearch, blastp, fimo,
	 * prodigal) from Java, it must be possible to locate bash scripts that run
	 * each program. If these cannot be located, PRISM will terminate.<br>
	 * <br>
	 * This method checks the prodigal.sh script only if prodigal gene search is
	 * enabled, the fimo.sh script only if ribosomal domain search is enabled,
	 * and the ecfp6_fcfp6.py script only if chemical similarity search is
	 * enabled.
	 * 
	 * @throws DependencyException
	 *             if a script does not exist
	 */
	public void checkScripts() throws DependencyException {
		// check blastp
		String blastFilepath = session.subDir("blast") + "blastp.sh";
		checkScriptFile(blastFilepath);

		// check hmmsearch
		String hmmFilepath = session.subDir("hmm") + "hmmsearch.sh";
		checkScriptFile(hmmFilepath);

		// check fimo
		if (requireFimo()) {
			String fimoFilepath = session.subDir("motifs") + "fimo.sh";
			checkScriptFile(fimoFilepath);
		}

		// check prodigal
		if (requireProdigal()) {
			String prodigalFilepath = session.subDir("orfs") + "prodigal.sh";
			checkScriptFile(prodigalFilepath);
		}

		// check fingerprinter
		if (requireFingerprinter()) {
			String fingerprintFilepath = session.subDir("tanimoto")
					+ "ecfp6_fcfp6.py";
			checkScriptFile(fingerprintFilepath);
		}
	}

	/**
	 * Check that a script file exists.
	 * 
	 * @param filepath
	 *            the filepath of the script to check
	 * @throws DependencyException
	 *             if no such file exists
	 */
	private void checkScriptFile(String filepath) throws DependencyException {
		File file = new File(filepath);
		String filename = file.getName();
		if (!file.exists())
			throw new DependencyException("Error: couldn't locate script "
					+ filename + " at filepath " + file.getAbsolutePath());
	}

	/**
	 * Check that all binary dependencies (blastp, hmmsearch, and optionally
	 * fimo, prodigal, and RDKit) are installed.
	 * 
	 * @throws DependencyException
	 * @throws InterruptedException
	 * @throws IOException
	 * @throws BadSmilesToFingerprinterException
	 */
	public void checkBinaryDependencies()
			throws IOException, InterruptedException, DependencyException,
			BadSmilesToFingerprinterException {
		checkTestsFolder();
		checkBlast();
		checkHmmsearch();
		if (requireFimo())
			checkFimo();
		if (requireProdigal())
			checkProdigal();
		if (requireFingerprinter())
			checkRDKit();
	}

	/**
	 * Make sure that the folder which contains files used to test the
	 * installation of binary dependencies exists.
	 * 
	 * @throws DependencyException
	 *             if the tests folder doesn't exist
	 */
	private void checkTestsFolder() throws DependencyException {
		String filepath = session.subDir("tests");
		File file = new File(filepath);
		if (!file.exists())
			throw new DependencyException("Error: couldn't locate "
					+ "tests folder at filepath " + file.getAbsolutePath());
	}

	/**
	 * Test for BLASTP installation by running a simple alignment of one
	 * sequence against a database also containing a single sequence. If no
	 * alignments can be retrieved then BLASTP is not installed.
	 * 
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DependencyException
	 */
	private void checkBlast()
			throws IOException, InterruptedException, DependencyException {
		session.listener().updateLastDetail("Ensuring BLASTP is installed...");

		String database = session.subDir("tests") + "blast" + File.separator
				+ "database";
		String query = session.subDir("tests") + "blast" + File.separator
				+ "query.fasta";

		BlastpSearch blast = new BlastpSearch(database, query, session);
		List<BlastSearchResult> results = blast.run(-1);
		if (results.size() == 0)
			throw new DependencyBinaryException(
					"Error: BLASTP is not installed!");
	}

	/**
	 * Test for hmmsearch installation by running a hidden Markov model against
	 * a small sequence database. If no hmm search results are produced then
	 * hmmsearch is not installed.
	 * 
	 * @throws DependencyBinaryException
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void checkHmmsearch() throws DependencyBinaryException, IOException,
			InterruptedException {
		session.listener()
				.updateLastDetail("Ensuring hmmsearch is installed...");

		String model = session.subDir("tests") + "hmm" + File.separator
				+ "database.hmm";
		String query = session.subDir("tests") + "hmm" + File.separator
				+ "query.fasta";
		HmmSearch hmmsearch = new HmmSearch(model, query, session);
		hmmsearch.run();
		List<HmmSearchResult> results = hmmsearch.reader().read();
		if (results.size() == 0)
			throw new DependencyBinaryException(
					"Error: hmmsearch is not installed!");
	}

	/**
	 * Test for prodigal installation by finding orfs in a small (~3 kb)
	 * cluster. If no orfs are identified then prodigal is not installed.
	 * 
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DependencyBinaryException
	 */
	private void checkProdigal() throws IOException, InterruptedException,
			DependencyBinaryException {
		session.listener()
				.updateLastDetail("Ensuring Prodigal is installed...");

		String filepath = session.subDir("tests") + "prodigal" + File.separator
				+ "streptide.fasta";
		Genome genome = new Genome(new File(filepath));
		FastaReader.readFastaFile(genome);
		try {
			ProdigalSearch prodigal = new ProdigalSearch(genome, session);
			prodigal.run();
		} catch (ProdigalSearchException e) {
			throw new DependencyBinaryException(
					"Error: Prodigal is not installed!");
		}
	}

	/**
	 * Test for FIMO installation by finding a motif in a small precursor
	 * peptide. If no motifs are identified, then FIMO is not installed.
	 * 
	 * @throws DependencyBinaryException
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void checkFimo() throws DependencyBinaryException, IOException,
			InterruptedException {
		session.listener().updateLastDetail("Ensuring FIMO is installed...");

		String query = session.subDir("tests") + "fimo" + File.separator
				+ "query.fasta";
		String motif = session.subDir("tests") + "fimo" + File.separator
				+ "motif.txt";

		FimoSearch search = new FimoSearch(query, motif, session);
		List<FimoSearchResult> results = search.run();
		if (results.size() == 0)
			throw new DependencyBinaryException(
					"Error: FIMO is not installed!");
	}

	/**
	 * Determine whether RDKit is installed by generating ECFP6 and FCFP6
	 * fingerprints for an amino acid. If no fingerprints are generated, then
	 * RDKit is not installed.
	 * 
	 * @throws DependencyBinaryException
	 * @throws IOException
	 * @throws BadSmilesToFingerprinterException
	 * @throws InterruptedException
	 */
	private void checkRDKit() throws DependencyBinaryException, IOException,
			BadSmilesToFingerprinterException, InterruptedException {
		session.listener().updateLastDetail("Ensuring RDKit is installed...");

		String directory = session.subDir("tanimoto");
		String smiles = "O=C(O)CN";
		List<Fingerprint> fingerprints = Ecfp6Fcfp6Fingerprinter
				.getCompoundFingerprints(smiles, directory);
		if (fingerprints.size() == 0)
			throw new DependencyBinaryException(
					"Error: RDKit is not installed!");
	}

	/**
	 * Determine whether FIMO is required as a dependency.
	 * 
	 * @return true if FIMO dependency checks must be performed
	 */
	private boolean requireFimo() {
		return config.ribosomal && config.scaffoldLimit > 0;
	}

	/**
	 * Determine whether Prodigal is required as a dependency.
	 * 
	 * @return true if Prodigal dependency checks must be performed
	 */
	private boolean requireProdigal() {
		return config.genePredictionModes
				.contains(GenePredictionModes.PRODIGAL);
	}

	/**
	 * Determine whether the RDKit and a Python script to execute it are
	 * required.
	 * 
	 * @return true if RDKit dependency checks must be performed
	 */
	private boolean requireFingerprinter() {
		return config.score && config.scaffoldLimit > 0;
	}

}
