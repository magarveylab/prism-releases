package ca.mcmaster.magarveylab.prism.util;

import java.io.File;
import java.util.List;

/**
 * Utilities class for string operations
 * @author skinnider
 *
 */
public class Strings {
	
	/**
	 * Get this string in capitalized format, with the first letter in uppercase and the rest in lowercase.
	 * @param input	string to capitalize
	 * @return		capitalized string
	 */
	public static String capitalize(String input) {
		return input.substring(0,1).toUpperCase() + input.substring(1).toLowerCase();
	}
	
	/**
	 * Get this string with the first letter capitalized.
	 * @param input	string to modify
	 * @return		modified string
	 */
	public static String capitalizeFirstLetter(String input) {
		return input.substring(0,1).toUpperCase() + input.substring(1);
	}

	/**
	 * Get the name and extension of a file, given its absolute path, i.e. removing the path to the file.
	 * @param file	filepath
	 * @return		filename
	 */
	public static String name(String file) {
		String[] split = file.split(File.separator);
		return split[split.length - 1];
	}
	
	/**
	 * Replace all the underscores in a string with spaces.
	 * @param input	input string
	 * @return		underscore-replaced string
	 */
	public static String removeUnderscores(String input) {
		return input.replace('_', ' ');
	}
	
	/**
	 * Test whether a list of strings contains a given query string.
	 * @param query		string to query
	 * @param strings	list of strings
	 * @return
	 */
	public static boolean contains(String query, List<String> strings) {
		boolean flag = false;
		for (String string : strings)
			if (string.equals(query))
				flag = true;
		return flag;
	}
	
	/**
	 * Format a list of strings as a comma-delineated list within a single string.
	 * This function transforms a list { "a", "b", "c" } -> "a, b, c"
	 * @param input		list of strings
	 * @return			list as a single, comma-delineated string 
	 */
	public static String arrayAsCommaDelinatedList(List<String> input) {
		StringBuffer out = new StringBuffer();
		int last = input.size() - 1;
		for (int i = 0; i < last; i++)
			out.append(input.get(i) + ", ");
		if (input.size() > 0)
			out.append(input.get(last));
		return out.toString();
	}

	/**
	 * Calculate an appropriate font size based on the length of an assembly domain's abbreviation. 
	 * @param abbreviation	the domain's abbreviation
	 */
	public static int calculateFontSize(String abbreviation) {
		int size = 10; // default
		String replaced = abbreviation.replace("<sub>", "").replace("</sub>", "").replace("&alpha;", "a")
				.replace("&Alpha;", "a").replace("&Beta;", "b").replace("&beta;", "b"); // strip HTML tags
		int length = replaced.length();
		if (length <= 2) {
			size = 16;
		} else if (length == 3) {
			size = 12;
		} else if (length == 4) {
			size = 10;
		} else if (length == 5) {
			size = 9;
		} else if (length > 5) {
			size = 8;
		}
		return size;
	}
	
}
