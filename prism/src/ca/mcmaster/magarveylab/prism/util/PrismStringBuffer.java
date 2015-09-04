package ca.mcmaster.magarveylab.prism.util;

/**
 * Extends the functionality of Java's default string buffer.
 * @author skinnider
 *
 */
public class PrismStringBuffer {

	private StringBuffer sb;
	
	public PrismStringBuffer(StringBuffer sb) {
		this.sb = sb;
	}
	
	public void appendLine(String line) {
		sb.append(line + "\n");
	}
	
	public void append(String line) {
		sb.append(line);
	}
	
	public String toString() {
		return sb.toString();
	}
	
}