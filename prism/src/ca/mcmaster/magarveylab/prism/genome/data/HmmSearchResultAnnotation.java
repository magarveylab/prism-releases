package ca.mcmaster.magarveylab.prism.genome.data;

/**
 * Data package for a hmmsearch domain annotation. 
 * @author skinnider
 *
 */
public class HmmSearchResultAnnotation {

	private int start;
	private int end;
	private double score;
	private String sequence;
	
	private HmmSearchResult row;
	
	public HmmSearchResultAnnotation(int start, int end, double score, HmmSearchResult row) {
		this.start = start;
		this.end = end;
		this.score = score;
		this.row = row;
	}
	
	/**
	 * Get the name of this annotation, expressed as the name of the parent <code>HmmRow</code> with 
	 * "<code>|start|end</code>" appended.
	 * @return	the name of this annotation
	 */
	public String name() {
		return row.name + "|" + start + "|" + end;
	}
		
	public int start() {
		return start;
	}
	
	public int end() {
		return end;
	}
	
	public double score() {
		return score;
	}
	
	public String sequence() {
		return sequence;
	}
	
	public HmmSearchResult row() {
		return row;
	}
}
