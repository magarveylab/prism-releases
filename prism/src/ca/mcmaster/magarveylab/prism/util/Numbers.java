package ca.mcmaster.magarveylab.prism.util;

import java.util.List;

/**
 * Utilities class for number operations.
 * @author skinnider
 *
 */
public class Numbers {

	/**
	 * Convert a number output by hmmsearch or blast as a String to a double.
	 * @param score		a single cell of numerical output
	 * @return			the number as a double
	 */
	public static double newDoubleFromString(String score) {
		double result;
		if (score.indexOf("e") != -1) {
			String[] temp = score.split("e");
			Double base = Double.parseDouble(temp[0]);
			Double exponent = Double.parseDouble(temp[1]);
			result = base * Math.pow(10,exponent);
		} else {
			result = Double.parseDouble(score);
		}
		return result;
	}
	
	/**
	 * Swap the values of two doubles.
	 * @param d1	first double
	 * @param d2	second double
	 */
	public static void swap(double d1, double d2) {
		double temp = d1;
		d1 = d2;
		d2 = temp;
	}
	
	public static int maximum(List<Integer> integers) {
		int max = 0;
		if (integers.size() > 0) {
			max = integers.get(0);
			for (int i : integers)
				if (i > max)
					max = i;
		}
		return max;
	}
	
	public static int minimum(List<Integer> integers) {
		int min = 0;
		if (integers.size() > 0) {
			min = integers.get(0);
			for (int i : integers)
				if (i < min)
					min = i;
		}
		return min;
	}
	
	public static String humanLength(Integer bytes){
	    int unit = 1000;
	    if (bytes < unit) return bytes + " B";
	    int exp = (int) (Math.log(bytes) / Math.log(unit));
	    String pre = "kMGTPE".charAt(exp-1) + "";
	    return String.format("%.1f %sbp", bytes / Math.pow(unit, exp), pre);
	}
	

}
