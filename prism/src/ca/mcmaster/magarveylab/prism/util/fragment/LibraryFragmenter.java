package ca.mcmaster.magarveylab.prism.util.fragment;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import ca.mcmaster.magarveylab.prism.data.Library;

public class LibraryFragmenter {

	private Library library;
	public LibraryFragmenter(Library library){
		this.library = library;
	}

	public void fragment(){

		library.scaffolds();
		if (library != null && library.scaffolds().size() > 0){
			List<String> fragMFString = new ArrayList<String>();
			for(String s : library.scaffolds()){		
				try{
					SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
					IMolecule mol = sp.parseSmiles(s);
					NRPFragmenter nrpf = new NRPFragmenter(mol);
					ArrayList<IMolecule> frags = nrpf.getAllFragments();

					for(IMolecule frag: frags){
						IMolecularFormula imf = MolecularFormulaManipulator.getMolecularFormula(frag);
						double mass = MolecularFormulaManipulator.getNaturalExactMass(imf);
						mass = mass*10000;
						mass = (double)((int)mass);
						mass = mass / 10000;
						fragMFString.add(new Double(mass).toString());
					}
				}
				catch(Exception e){
					System.out.println("Invalid smiles in library: " + s);
					continue;
				}
			}
			// Getting the unique molecular formulas
			Set<String> unique_mass = new HashSet<String>();
			unique_mass.addAll(fragMFString);
			Map<String, Integer> counts = new LinkedHashMap<String, Integer>();
			for(String mass: unique_mass){
				int occurrences = Collections.frequency(fragMFString, mass);
				counts.put(mass, occurrences);
			}
			// Sorting the mass according to frequency
			List<Entry<String, Integer>> counts2 = new ArrayList<Entry<String, Integer>>( counts.entrySet());
			Collections.sort(counts2, new Comparator<Entry<String, Integer>>(){
				@Override
				public int compare(Entry<String, Integer> o1, Entry<String, Integer> o2) {
					Integer val0 = o1.getValue();
					Integer val1 = o2.getValue();
					if(val0.equals(val1)){
						return 0;
					}
					else if(val1 > val0){
						return 1;
					}
					else if(val1 < val0){
						return -1;
					}
					else{
						return 0;
					}
				}
			});

			// The max value for i is the total number of fragments we want
			// Change in here if we want to change the number of fragments we want
			int fragmentCutoff = 50;
			if(counts.size() < 50){
				fragmentCutoff = counts.size();
			}
			List<Double> sorted_mass = new ArrayList<Double>();
			for(int i = 0; i < fragmentCutoff; i++){
				sorted_mass.add( Double.parseDouble(counts2.get(i).getKey()) );
			}
			library.setFragments(sorted_mass);
		}
	}
}
