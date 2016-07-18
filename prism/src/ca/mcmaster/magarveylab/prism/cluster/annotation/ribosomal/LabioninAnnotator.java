package ca.mcmaster.magarveylab.prism.cluster.annotation.ribosomal;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.analysis.RibosomalClusterAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.cluster.reactions.RibosomalUtil;
import ca.mcmaster.magarveylab.prism.combinatorialization.Permutations;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;

/**
 * Find potential substrates for the lantibiotic LanKC cyclase/dehydratase,
 * which catalyzes labionin bond formation (i.e., lantithione bonds with the
 * enolate attacking a second Dha/Dhb residue).
 * 
 * @author skinnider
 *
 */
public class LabioninAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.LanKC };
	}

	public List<SubstrateSet> findSubstrates(Domain domain,
			List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		// get all serine/threonine
		List<Module> st = new ArrayList<Module>();
		for (Module module : permutation)
			if (module.scaffold() != null
					&& module.scaffold().topSubstrate() != null
					&& module.scaffold().topSubstrate().type() != null
					&& (module.scaffold().topSubstrate().type() == ProteinogenicAminoAcids.SERINE || module
							.scaffold().topSubstrate().type() == ProteinogenicAminoAcids.THREONINE))
				st.add(module);
//		st.addAll(RibosomalClusterAnalyzer.getModules(
//				ProteinogenicAminoAcids.SERINE, permutation));
//		st.addAll(RibosomalClusterAnalyzer.getModules(
//				ProteinogenicAminoAcids.THREONINE, permutation));

		// get all cysteine
		List<Module> cys = RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.CYSTEINE, permutation);

		// first, get all possible combinations of labionins, assuming that all
		// possible labionins get formed
		List<SubstrateSet> labionins = getLabionins(st, cys, permutation);

		// next, for each combination of labionins, get potential subsets of
		// (methyl)lanthionines
		List<SubstrateSet> labioninsWithLanthionines = addLanthionines(
				labionins, domain, permutation, cluster);

		// finally, for each labionin/(methyl)lanthionine combination, get
		// potential subsets of remaining residues
		List<SubstrateSet> substrates = addDehydrations(
				labioninsWithLanthionines, permutation);

		// print this combination 
		for (SubstrateSet substrate : substrates) {
			StringBuffer sb = new StringBuffer();
			for (Module module : substrate.getAllModules())
				if (module == null) {
					sb.append(" -- ");
				} else {
					sb.append(module.scaffold().topSubstrate().type()
							.abbreviation()
							+ permutation.indexOf(module) + " ");
				}
			System.out.println(sb.toString());
		}
		
		return substrates;
	}
	
	public static List<SubstrateSet> addLanthionines(List<SubstrateSet> labionins,
			Domain domain, List<Module> permutation, Cluster cluster) throws InvalidSmilesException, IOException {
		List<SubstrateSet> lanthionines = new ArrayList<SubstrateSet>();
		if (labionins.size() == 0) {
			SubstrateSet empty = new SubstrateSet();
			lanthionines = addLanthionines(empty, domain, permutation, cluster);
		} else {
			for (SubstrateSet labionin : labionins) {
				List<SubstrateSet> substrates = addLanthionines(labionin, domain, permutation, cluster);
				lanthionines.addAll(substrates);
			}
		}
		return lanthionines;
	}
	
	private static List<SubstrateSet> addLanthionines(SubstrateSet labionin,
			Domain domain, List<Module> permutation, Cluster cluster) throws InvalidSmilesException, IOException {
		// get all unused modules 
		List<Module> unused = new ArrayList<Module>();
		for (Module module : permutation)
			if (!labionin.getAllModules().contains(module))
				unused.add(module);
			
		LanthionineAnnotator annotator = new LanthionineAnnotator();
		List<SubstrateSet> lanthionines = annotator.findSubstrates(domain, unused, cluster);
		
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		if (lanthionines.size() == 0) {
			SubstrateSet substrate = new SubstrateSet(labionin);
			substrate.add(null); // spacer for Reaction implementation 
			substrates.add(substrate);
		} else {
			for (SubstrateSet lanthionine : lanthionines) {
				SubstrateSet substrate = new SubstrateSet(labionin);
				substrate.add(null); // spacer for Reaction implementation 
				substrate.addAll(lanthionine.getAllModules());
				substrates.add(substrate);
			}
		}
		return substrates;
	}

	public static List<SubstrateSet> addDehydrations(
			List<SubstrateSet> labLans, List<Module> permutation)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		if (labLans.size() == 0) {
			SubstrateSet empty = new SubstrateSet();
			substrates = addDehydrations(empty, permutation);
		} else {
			for (SubstrateSet labLan : labLans) {
				List<SubstrateSet> s = addDehydrations(labLan, permutation);
				substrates.addAll(s);
			}
		}
		return substrates;
	}

	public static List<SubstrateSet> addDehydrations(SubstrateSet labLan,
			List<Module> permutation) throws InvalidSmilesException,
			IOException {
		// get all unused modules
		List<Module> unused = new ArrayList<Module>();
		for (Module module : permutation)
			if (!labLan.getAllModules().contains(module))
				unused.add(module);
		
		// get ser/thr
		List<Module> st = new ArrayList<Module>();
		st.addAll(RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.SERINE, unused));
		st.addAll(RibosomalClusterAnalyzer.getModules(ProteinogenicAminoAcids.THREONINE, unused));
		
		// get all subsets 
		List<List<Module>> subsets = RibosomalUtil.getSubsetsWithMinSize(st, 0);

		// convert to substrates
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		if (subsets.size() == 0) {
			SubstrateSet substrate = new SubstrateSet(labLan);
			substrate.add(null); // spacer for Reaction implementation 
			substrates.add(substrate);
		} else {
			for (List<Module> dehydration : subsets) {
				SubstrateSet substrate = new SubstrateSet(labLan);
				substrate.add(null);
				substrate.addAll(dehydration);
				substrates.add(substrate);
			}
		}
		return substrates;
	}

	public static List<SubstrateSet> getLabionins(List<Module> serines, List<Module> cysteines,
			List<Module> permutation) {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();
		
		int c = cysteines.size() * 2;
		int s = serines.size();
		if (s >= c) {
			// get all permutations of serine of size [# cysteines] * 2
			List<int[]> serThrPermutations = Permutations.permutations(s, c, 2_000);
			permutationLoop: for (int[] serThrPermutation : serThrPermutations) {
				// make sure serines are in order
				if (!ordered(serThrPermutation))
					continue permutationLoop;
				
				// convert int array to SubstrateSet 
				SubstrateSet substrate = new SubstrateSet(); 
				for (int i = 0; i < serThrPermutation.length - 1; i += 2) {
					Module serThr1 = serines.get(serThrPermutation[i] - 1);
					Module serThr2 = serines.get(serThrPermutation[i+1] - 1);
					Module cys = cysteines.get(i/2);

					int serThr1Idx = permutation.indexOf(serThr1);
					int serThr2Idx = permutation.indexOf(serThr2);
					int cysIdx = permutation.indexOf(cys);

					// ring 1 must have 2 aa between
					if (Math.abs(serThr2Idx - serThr1Idx) != 3)
						continue permutationLoop;
					
					// ring 2 must have 2-5 aa between
					if (Math.abs(cysIdx - serThr2Idx) < 3
							|| Math.abs(cysIdx - serThr2Idx) > 6) 
						continue permutationLoop;
					
					// make sure cysteines are in order
					if (!(cysIdx > serThr1Idx && cysIdx > serThr2Idx))
						continue permutationLoop;
					
					// make sure all previous cysteines are in order, too
					for (int j = 0; j < i; j += 2) {
						Module previousCys = cysteines.get(j/2);
						if (permutation.indexOf(previousCys) > serThr1Idx
								|| permutation.indexOf(previousCys) > serThr2Idx)
							continue permutationLoop;
					}
					
					substrate.add(serThr2);
					substrate.add(cys);
					substrate.add(serThr1);
				}
				
				if (substrate.size() > 0)
					substrates.add(substrate);
				
				System.out.println("Substrates size = " + substrates.size());
			}
		} else {
			// if there are fewer serines than cysteines, only use as many as possible
			System.out.println("[LabioninAnnotator] Permuting " + c + " cysteines and " + s + " serines");
			
			List<int[]> serThrPermutations = Permutations.permutations(s, s, 500);
			List<int[]> cysPermutations = Permutations.permutations(cysteines.size(), s/2, 500);
			
			permutationLoop: for (int[] serThrPermutation : serThrPermutations) {
				if (!ordered(serThrPermutation))
					continue permutationLoop;
				
				cysPermutationLoop: for (int[] cysPermutation : cysPermutations) {
					if (!ordered(cysPermutation)) 
						continue cysPermutationLoop;
					
					// convert int array to SubstrateSet 
					SubstrateSet substrate = new SubstrateSet(); 
					for (int i = 0; i < serThrPermutation.length - 1; i += 2) {
						Module serThr1 = serines.get(serThrPermutation[i] - 1);
						Module serThr2 = serines.get(serThrPermutation[i+1] - 1);
						Module cys = cysteines.get(cysPermutation[i/2] - 1); 
						
						int serThr1Idx = permutation.indexOf(serThr1);
						int serThr2Idx = permutation.indexOf(serThr2);
						int cysIdx = permutation.indexOf(cys);
						
						// ring 1 must have 2 aa between
						if (Math.abs(serThr2Idx - serThr1Idx) != 3)
							continue permutationLoop;
						
						// ring 2 must have 2-5 aa between
						if (Math.abs(cysIdx - serThr2Idx) < 3
								|| Math.abs(cysIdx - serThr2Idx) > 6) 
							continue permutationLoop;
						
						// make sure cysteines are in order 
						if (!(permutation.indexOf(cys) > permutation.indexOf(serThr1)
								&& permutation.indexOf(cys) > permutation.indexOf(serThr2))) 
							continue cysPermutationLoop;
						
						// make sure all previous cysteines are in order, too
						for (int j = 0; j < i; j += 2) {
							Module previousCys = cysteines.get(cysPermutation[j/2] - 1);
							if (permutation.indexOf(previousCys) > permutation.indexOf(serThr1)
									|| permutation.indexOf(previousCys) > permutation.indexOf(serThr2)) 
								continue cysPermutationLoop;
						}
						
						substrate.add(serThr2);
						substrate.add(cys);
						substrate.add(serThr1);
					}
					substrates.add(substrate);
				}
			}
		}
		
		// add an empty substrate, to account for no labionin formation 
		if (substrates.size() == 0)
			substrates.add(new SubstrateSet());
		
		return substrates;
	}
	
	public static boolean ordered(int[] intArray) {
		boolean flag = true;
		for (int i = 0; i < intArray.length; i++) {
			int currentIdx = intArray[i];
			int j = i;
			while (j > 0) {
				j--;
				int prevIdx = intArray[j];
				if (prevIdx > currentIdx) {
					flag = false;
				}
			}
		}
		return flag;
	}

}
