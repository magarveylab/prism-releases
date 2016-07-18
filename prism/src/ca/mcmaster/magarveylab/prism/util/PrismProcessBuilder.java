package ca.mcmaster.magarveylab.prism.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 * Launches an external application using Java's ProcessBuilder class and
 * returns output in a BufferedReader.
 * 
 * @author skinnider
 *
 */
public class PrismProcessBuilder {

	private String[] cmd;

	/**
	 * Instantiate a new custom process builder.
	 * 
	 * @param cmd
	 *            command to run in command line
	 */
	public PrismProcessBuilder(String[] cmd) {
		this.cmd = cmd;
	}

	/**
	 * Launch an external application using Java's ProcessBuilder.
	 * 
	 * @return a reader for the command line output
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public BufferedReader run() throws IOException, InterruptedException {
		ProcessBuilder pb = new ProcessBuilder(cmd);
		pb.redirectErrorStream(true);
		Process shell = pb.start();
		InputStream shellIn = shell.getInputStream();
		BufferedReader br = new BufferedReader(new InputStreamReader(shellIn));
		// shell.waitFor();
		return br;
	}

}
