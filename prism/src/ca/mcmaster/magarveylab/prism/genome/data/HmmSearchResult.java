package ca.mcmaster.magarveylab.prism.genome.data;

import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.prism.util.Numbers;

/**
 * Data package for a single row in a hmmsearch score table. The results of a hmmsearch for a single sequence. 
 * @author skinnider
 *
 */
public class HmmSearchResult {

	private List<HmmSearchResultAnnotation> annotations = new ArrayList<HmmSearchResultAnnotation>();
	
	public double fullSequenceEValue, fullSequenceScore, fullSequenceBias,
		bestDomainEValue, bestDomainScore, bestDomainBias;
	public String description, name;
	
	/**
	 * Initialize a new HmmRow. 
	 * @param fullSequenceEValue	the full sequence E value
	 * @param fullSequenceScore		the full sequence score
	 * @param fullSequenceBias		the full sequence bias
	 * @param bestDomainEValue		the best domain E value
	 * @param bestDomainScore		the best domain score
	 * @param bestDomainBias		the best domain bias
	 * @param description			the 
	 * @param name
	 */
	public HmmSearchResult(String fullSequenceEValue, String fullSequenceScore,
			String fullSequenceBias, String bestDomainEValue, String bestDomainScore, 
			String bestDomainBias, String name, String description) {
		this.fullSequenceEValue = Numbers.newDoubleFromString(fullSequenceEValue);
		this.fullSequenceScore = Numbers.newDoubleFromString(fullSequenceScore);
		this.fullSequenceBias = Numbers.newDoubleFromString(fullSequenceBias);
		this.bestDomainEValue = Numbers.newDoubleFromString(bestDomainEValue);
		this.bestDomainScore = Numbers.newDoubleFromString(bestDomainScore);
		this.bestDomainBias = Numbers.newDoubleFromString(bestDomainBias);
		this.description = description;
		this.name = name;
	}
	
	/**
	 * Get all the hmmsearch annotations associated with this result. 
	 * @return	all hmmsearch annotations associated with this result
	 */
	public List<HmmSearchResultAnnotation> annotations() {
		return annotations;
	}
	
	/**
	 * Associate a new hmmsearch annotation to this result. 
	 * @param annotation	the hmmsearch annotation to associate with this result
	 */
	public void addAnnotation(HmmSearchResultAnnotation annotation) {
		annotations.add(annotation);
	}
	
}
