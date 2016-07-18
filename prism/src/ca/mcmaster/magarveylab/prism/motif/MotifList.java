package ca.mcmaster.magarveylab.prism.motif;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;

import ca.mcmaster.magarveylab.enums.RibosomalPrecursorMotifs;
import ca.mcmaster.magarveylab.prism.motif.CompareMotifs;

/**
 * A list of motifs identified within a sequence using FIMO.
 * 
 * @author Robyn Edgar, skinnider 
 *
 */
public class MotifList implements Iterable<Motif> {

	private List<Motif> motifs;

	/**
	 * Instantiate a new motif list.
	 * 
	 */
	public MotifList() {
		this.motifs = new ArrayList<Motif>();
	}

	public void addMotif(Motif motif) {
		motifs.add(motif);
	}

	public void removeMotif(int index) {
		motifs.remove(index);
	}

	public List<Motif> getMotifs() {
		return motifs;
	}

	/**
	 * Get the motif at a given index.
	 * 
	 * @param index
	 *            index of the motif to get
	 * @return the motif at that index
	 */
	public Motif getMotifAtIndex(int index) {
		return motifs.get(index);
	}
	
	/**
	 * Get a sublist containing only motifs of the given type.
	 * 
	 * @param type
	 *            type of motifs to get
	 * @return a sublist of this motif list, containing only motifs of the given
	 *         type
	 */
	public MotifList getMotifs(RibosomalPrecursorMotifs type) {
		MotifList motifs = new MotifList();
		for (Motif motif : this.motifs)
			if (motif.getType() == type)
				motifs.addMotif(motif);
		return motifs;
	}

	public int size() {
		return motifs.size();
	}

	public void sortMotifsByStart() {
		Collections.sort(motifs, new CompareMotifs());
	}

	public void removeOverlappingMotifs() {
		if (motifs.isEmpty()) {
			System.out.println("[MotifList] Could not remove "
					+ "overlapping motifs: list is empty");
		} else {
			Iterator<Motif> itr1 = motifs.iterator();
			Iterator<Motif> itr2 = motifs.iterator();
			while (itr1.hasNext()) {
				Motif motif1 = itr1.next();
				while (itr2.hasNext()) {
					Motif motif2 = itr2.next();
					if (motif1.equals(motif1))
						continue;
					if (motif1.overlaps(motif2)) 
						if (motif1.getPValue() < motif2.getPValue()) {
							itr2.remove();
						} else {
							itr1.remove();
						}
				}
			}
		}
	}

	/**
	 * Get the motif with the lowest p-value from a list of motifs.
	 * 
	 * @return the motif with the lowest p-value in the list, or
	 *         null if the list is empty 
	 */
	public Motif getBestMotif() {
		double pValue = 1;
		Motif motif = null;
		for (Motif m : motifs)
			if (m.getPValue() < pValue) {
				motif = m;
				pValue = m.getPValue();
			}
		return motif;
	}

	/**
	 * Get the first motif from a list of motifs.
	 * 
	 * @return the first motif in the list, or null if the list is empty
	 */
	public Motif getFirstMotif() {
		int start = -1;
		Motif motif = null;
		for (Motif m : motifs)
			if (m.getStart() < start || start == -1) {
				motif = m;
				start = m.getStart();
			}
		return motif;
	}

	/**
	 * Get the motif with the lowest p-value from a list of motifs, of a given
	 * type.
	 * 
	 * @param type
	 *            type of motif to search for
	 * @return the motif of the given type and with the lowest p-value in the
	 *         list, or null if no such motif exists
	 */
	public Motif getBestMotif(RibosomalPrecursorMotifs type) {
		double pValue = 1;
		Motif motif = null;
		for (Motif m : motifs)
			if (m.getType() == type && m.getPValue() < pValue) {
				motif = m;
				pValue = m.getPValue();
			}
		return motif;
	}
	
	/**
	 * Get the first occurrence of a motif of a given type from a list of motifs.
	 * 
	 * @param type
	 *            type of motif to search for
	 * @return the first motif of the given type in the
	 *         list, or null if no such motif exists
	 */
	public Motif getFirstMotif(RibosomalPrecursorMotifs type) {
		int start = -1;
		Motif motif = null;
		for (Motif m : motifs)
			if (m.getType() == type && (m.getStart() < start || start == -1)) {
				motif = m;
				start = m.getStart();
			}
		return motif;
	}
	
	public Iterator<Motif> iterator() {
		return motifs.iterator();
	}
	
	public boolean isEmpty() {
		return motifs.isEmpty();
	}
	
	public boolean contains(Motif motif) {
		return motifs.contains(motif);
	}
	
	/**
	 * Check whether this list contains a motif of a given type.
	 * 
	 * @param type
	 *            motif type to check for
	 * @return true if the list contains a motif of that type
	 */
	public boolean contains(RibosomalPrecursorMotifs type) {
		boolean flag = false;
		for (Motif m : motifs)
			if (m.getType() == type)
				flag = true;
		return flag;
	}
	
}