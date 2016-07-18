package ca.mcmaster.magarveylab.prism.motif;

import java.util.Comparator;

/**
 * Motif comparison by start of motif.
 * 
 * @author Robyn Edgar
 *
 */
public class CompareMotifs implements Comparator<Motif> {

	@Override
	public int compare(Motif motif1, Motif motif2){
		return motif1.getStart() - motif2.getStart();
	}

	
}