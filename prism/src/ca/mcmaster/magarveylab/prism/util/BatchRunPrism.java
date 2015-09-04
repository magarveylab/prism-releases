package ca.mcmaster.magarveylab.prism.util;

import java.io.File;

import ca.mcmaster.magarveylab.prism.Prism;
import ca.mcmaster.magarveylab.prism.PrismDesktop;
import ca.mcmaster.magarveylab.prism.web.PrismConfig;
import ca.mcmaster.magarveylab.wasp.session.Session;

public class BatchRunPrism {
	
	/**
	 * Location of the cluster database directory. 
	 */
	private static String dir = "/Users/michaelskinnider/Desktop/accuracy_analysis/output_2/glycopeptide/";

	/**
	 * Extension of the FASTA files to run on. 
	 */
	private static String extension = ".fasta";
	
	public static void main(String[] args) throws Exception {
		// optionally set directory variable
		if (args.length > 0)
			dir = args[0];
		if (args.length > 1)
			extension = args[1];
			
		System.out.println("Running PRISM on all " + extension + " files in directory " + dir);
		
		// iterate over all .fasta files in directory 
		File[] files = Files.getDirectoryFiles(dir, extension);
		for (File file : files) {
			String filepath = file.getAbsolutePath();
			System.out.println("Starting new PRISM search on file " + filepath);
			
			PrismDesktop pd = new PrismDesktop(filepath);
			Prism prism = pd.search();

			Session session = pd.search().session();
			Runtime rt = Runtime.getRuntime();
			rt.exec("cp " + filepath + " " + session.dir() + Strings.name(filepath));
			
			// set config window 
			PrismConfig config = prism.config();
			config.window = 100_000;

			// execute search
			prism.run();
		}
		
		System.exit(0);
	}


}
