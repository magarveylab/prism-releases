package ca.mcmaster.magarveylab.prism.fasta;

/**
 * Remove unsupported nucleotide codes from DNA sequences.  	
 * @author skinnider
 *
 */
public class SequenceCleaner {
	
	private final static char[] valid = { 
		'k', 'r', 'y', 'm', 's', 'w', 'b', 'd', 'v', 'h', 'a', 'c', 'g', 't', 'n',
		'K', 'R', 'Y', 'M', 'S', 'W', 'B', 'D', 'V', 'H', 'A', 'C', 'G', 'T', 'N'
	};

	/**
	 * Remove all unsupported nucleotide codes (RYKMSWBDHV) from a sequence.
	 * @param seq	the sequence to clean
	 * @return		the cleaned sequence
	 */
	public static String clean(String seq) {
		StringBuilder sequence = new StringBuilder(seq);
		
		for (int i = 0; i < seq.length(); i++) {
			char c = sequence.charAt(i);
			if (!isValid(c))
				replaceN(i, sequence);
			if (c == 'k' || c == 'K')
				replaceK(i, sequence);
			if (c == 'r' || c == 'R')
				replaceR(i, sequence);
			if (c == 'y' || c == 'Y')
				replaceY(i, sequence);
			if (c == 'm' || c == 'M')
				replaceM(i, sequence);
			if (c == 's' || c == 'S')
				replaceS(i, sequence);
			if (c == 'w' || c == 'W')
				replaceW(i, sequence);
			if (c == 'b' || c == 'B')
				replaceB(i, sequence);
			if (c == 'd' || c == 'D')
				replaceD(i, sequence);
			if (c == 'v' || c == 'V')
				replaceV(i, sequence);
			if (c == 'h' || c == 'H')
				replaceH(i, sequence);
		}
		return sequence.toString();
	}
	
	/**
	 * Test whether the array of valid characters contains a given FASTA nucleotide character.
	 * @param valid		array of valid FASTA input characters
	 * @param c			character to test
	 * @return			false if this is invalid
	 */
	public static boolean isValid(char c) {
		boolean flag = false;
		for (char ch : valid)
			if (ch == c)
				flag = true;
		return flag;
	}
	
	/**
	 * Replace an invalid character with N. 
	 * @param i			index of the invalid character within the sequence
	 * @param sequence	the sequence
	 */
	public static void replaceN(int i, StringBuilder sequence) {
		sequence.setCharAt(i, 'N');
	}
	
	/**
	 * Replace a K (ketone nucleotide) with either G or T. 
	 * @param i			index of the K within the sequence
	 * @param sequence	the sequence
	 */
	public static void replaceK(int i, StringBuilder sequence) {
		int random = (int) Math.floor(Math.random() * 2);
		if (random == 1) {
			sequence.setCharAt(i, 'G');
		} else {
			sequence.setCharAt(i, 'T');
		}
	}
	
	/**
	 * Replace a R (purine nucleotide) with either A or G.
	 * @param i			index of the R within the sequence
	 * @param sequence	the sequence
	 */
	public static void replaceR(int i, StringBuilder sequence) {
		int random = (int) Math.floor(Math.random() * 2);
		if (random == 1) {
			sequence.setCharAt(i, 'A');
		} else {
			sequence.setCharAt(i, 'G');
		}
	}

	/**
	 * Replace a Y (pyrimidine nucleotide) with either T or C.
	 * @param i			index of the R within the sequence
	 * @param sequence	the sequence
	 */
	public static void replaceY(int i, StringBuilder sequence) {
		int random = (int) Math.floor(Math.random() * 2);
		if (random == 1) {
			sequence.setCharAt(i, 'T');
		} else {
			sequence.setCharAt(i, 'C');
		}
	}
	
	/**
	 * Replace a M (amino nucleotide) with either A or C.
	 * @param i			index of the R within the sequence
	 * @param sequence	the sequence
	 */
	public static void replaceM(int i, StringBuilder sequence) {
		int random = (int) Math.floor(Math.random() * 2);
		if (random == 1) {
			sequence.setCharAt(i, 'A');
		} else {
			sequence.setCharAt(i, 'C');
		}
	}

	/**
	 * Replace a S (strong interaction nucleotide) with either C or G.
	 * @param i			index of the R within the sequence
	 * @param sequence	the sequence
	 */
	public static void replaceS(int i, StringBuilder sequence) {
		int random = (int) Math.floor(Math.random() * 2);
		if (random == 1) {
			sequence.setCharAt(i, 'G');
		} else {
			sequence.setCharAt(i, 'C');
		}
	}

	/**
	 * Replace a W (weak interaction nucleotide) with either A or T.
	 * @param i			index of the R within the sequence
	 * @param sequence	the sequence
	 */
	public static void replaceW(int i, StringBuilder sequence) {
		int random = (int) Math.floor(Math.random() * 2);
		if (random == 1) {
			sequence.setCharAt(i, 'T');
		} else {
			sequence.setCharAt(i, 'A');
		}
	}
	
	/**
	 * Replace a B (not-A nucleotide) with either C, T, or G.
	 * @param i			index of the R within the sequence
	 * @param sequence	the sequence
	 */
	public static void replaceB(int i, StringBuilder sequence) {
		int random = (int) Math.floor(Math.random() * 3);
		if (random == 1) {
			sequence.setCharAt(i, 'T');
		} else if (random == 2) {
			sequence.setCharAt(i, 'G');
		} else {
			sequence.setCharAt(i, 'C');
		}
	}

	/**
	 * Replace a D (not-C nucleotide) with either A, T, or G.
	 * @param i			index of the R within the sequence
	 * @param sequence	the sequence
	 */
	public static void replaceD(int i, StringBuilder sequence) {
		int random = (int) Math.floor(Math.random() * 3);
		if (random == 1) {
			sequence.setCharAt(i, 'T');
		} else if (random == 2) {
			sequence.setCharAt(i, 'G');
		} else {
			sequence.setCharAt(i, 'A');
		}
	}
	
	/**
	 * Replace a H (not-G nucleotide) with either C, T, or A.
	 * @param i			index of the R within the sequence
	 * @param sequence	the sequence
	 */
	public static void replaceH(int i, StringBuilder sequence) {
		int random = (int) Math.floor(Math.random() * 3);
		if (random == 1) {
			sequence.setCharAt(i, 'T');
		} else if (random == 2) {
			sequence.setCharAt(i, 'A');
		} else {
			sequence.setCharAt(i, 'C');
		}
	}
	
	/**
	 * Replace a V (not-T nucleotide) with either C, A, or G.
	 * @param i			index of the R within the sequence
	 * @param sequence	the sequence
	 */
	public static void replaceV(int i, StringBuilder sequence) {
		int random = (int) Math.floor(Math.random() * 3);
		if (random == 1) {
			sequence.setCharAt(i, 'A');
		} else if (random == 2) {
			sequence.setCharAt(i, 'G');
		} else {
			sequence.setCharAt(i, 'C');
		}
	}

}
