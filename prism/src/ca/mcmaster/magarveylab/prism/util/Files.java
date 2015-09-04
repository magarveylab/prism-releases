package ca.mcmaster.magarveylab.prism.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;

/**
 * Utilities class for file operations.
 * 
 * @author skinnider
 *
 */
public class Files { 
	
	/**
	 * Get all files in a directory with the given extension.
	 * 
	 * @param directory
	 *            directory to search
	 * @param extension
	 *            extension to search for
	 * @return all files in query directory with that extension
	 */
	public static File[] getDirectoryFiles(final String directory,
			final String extension) {
		FilenameFilter filter = getFilter(extension);
		File[] files = new File(directory).listFiles(filter);
		return files;
	}

	/**
	 * Get all files in a directory, regardless of extension.
	 * 
	 * @param directory
	 *            directory to search
	 * @return all files in that directory
	 */
	public static File[] getAllFiles(final String directory) {
		File dir = new File(directory);
		File[] files = dir.listFiles();
		return files;
	}

	/**
	 * Get the FilenameFilter interface object used to consider files with a
	 * given extension.
	 * 
	 * @param extension
	 *            the file extension to filter by
	 * @return appropriate FilenameFilter
	 */
	public static FilenameFilter getFilter(final String extension) {
		FilenameFilter textFilter = new FilenameFilter() {
			public boolean accept(File dir, String name) {
				String lowercaseName = name.toLowerCase();
				if (lowercaseName.endsWith(extension)) {
					return true;
				} else {
					return false;
				}
			}
		};
		return textFilter;
	}
	
	/**
	 * Read the contents of a file into a string.
	 * 
	 * @param filepath
	 *            file to read
	 * @return entire file contents as string
	 * @throws IOException
	 */
	public static String readFile(String filepath) throws IOException {
		StringBuffer sb = new StringBuffer();
		BufferedReader br = new BufferedReader(new FileReader(filepath));
		String line;
		while ((line = br.readLine()) != null)
			sb.append(line + "\n");
		br.close();
		return sb.toString();
	}
	
	/**
	 * Write a string to a file.
	 * 
	 * @param file
	 *            string to write
	 * @param filepath
	 *            file to write to
	 * @throws IOException
	 */
	public static void writeFile(String file, String filepath)
			throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(filepath));
		bw.append(file);
		bw.close();
	}
	
	/**
	 * Get the names (not filepaths!) of all subdirectories within a parent
	 * directory.
	 * 
	 * @param directory
	 *            parent directory
	 * @return names of all directories within parent directory
	 */
	public static String[] getAllDirectories(String directory) {
		File parent = new File(directory);
		String[] directories = parent.list(new FilenameFilter() {
			@Override
			public boolean accept(File current, String name) {
				return new File(current, name).isDirectory();
			}
		});
		return directories;
	}

}
