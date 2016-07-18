package ca.mcmaster.magarveylab.prism.cluster.analysis;

import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.enums.Codons;
import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.domains.GenericDomains;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;
import ca.mcmaster.magarveylab.prism.combinatorialization.Combinations;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.Orf;
import ca.mcmaster.magarveylab.prism.data.Substrate;
import ca.mcmaster.magarveylab.prism.util.exception.RibosomalScaffoldGenerationException;

/**
 * Performs analysis of the biosynthetic gene clusters for ribosomally
 * synthesized and post-translationally modified peptides (RiPPs).
 * 
 * @author skinnider
 *
 */
public class RibosomalClusterAnalyzer {
		
	/**
	 * Get all ribosomal modules of a certain amino acid.
	 * @param aminoAcid		amino acid to query
	 * @param permutation	a module permutation to analyze 
	 * @return				all ribosomal modules of the query amino acid type 
	 */
	public static List<Module> getModules(ProteinogenicAminoAcids aminoAcid, List<Module> permutation) {
		List<Module> modules = new ArrayList<Module>();
		for (Module module : permutation) {
			if (module.type() == ModuleTypes.RIBOSOMAL) {
				Domain scaffold = module.scaffold();
				if (scaffold == null) {
					System.out.println("Error: could not get scaffold domain for ribosomal module");
					continue;					
				} else if (scaffold.topSubstrate() == null) {
					System.out.println("Error: could not get substrate for ribosomal module");
					continue;										
				}
				if (scaffold.topSubstrate().type() == aminoAcid) 
					modules.add(module);
			}
		}
		return modules;
	}

	/**
	 * Generate a list of ribosomal module permutations from the cleaved
	 * propeptide. Multiple permutations may be generated when one or more
	 * nucleotides within the propeptide is unknown (i.e. 'N' in the nucleotide
	 * sequence).
	 * 
	 * @param propeptide
	 *            the sequence of the cleaved propeptide
	 * @param orf
	 *            the parent open reading frame
	 * @return a list of ribosomal module permutations for scaffold generation
	 * @throws RibosomalScaffoldGenerationException
	 */
	public static List<List<Module>> generateRibosomalModulePermutations(String propeptide, Orf orf) 
			throws RibosomalScaffoldGenerationException {
		// get domains from propeptide
		List<Domain> domains = getDomainsFromPropeptide(propeptide, orf);
		// combinatorialize substrates
		List<List<Substrate>> substrates = getSubstrateCombinations(domains);
		// convert to modules
		List<List<Module>> permutations = convertToModulePermutations(substrates, domains);
		return permutations;
	}
	
	/**
	 * Convert a list of substrate combinations into a list of ribosomal module
	 * permutations.
	 * 
	 * @param substrates
	 *            list of substrate combinations
	 * @param domains
	 *            list of domains from which the substrate combinations were
	 *            derived
	 * @return list of module permutations
	 */
	public static List<List<Module>> convertToModulePermutations(List<List<Substrate>> substrates, 
			List<Domain> domains) {
		List<List<Module>> modulePermutations = new ArrayList<List<Module>>();
		for (List<Substrate> substrateList : substrates) {
			List<Module> modulePermutation = new ArrayList<Module>();
			for (Substrate substrate : substrateList)
				for (Domain domain : domains)
					if (domain.substrates().indexOf(substrate) != -1) {
						Domain d = new Domain(domain);
						d.substrates().clear();
						d.addSubstrate(substrate);
						Module m = new Module(ModuleTypes.RIBOSOMAL);
						m.add(d);
						modulePermutation.add(m);
					}
			modulePermutations.add(modulePermutation);
			System.out.println("[RibosomalClusterAnalyzer] Generated " 
					+ modulePermutation.size() + "-module permutation");
		}
		System.out.println("[RibosomalClusterAnalyzer] Generated " 
				+ modulePermutations.size() + " module permutations");
		return modulePermutations;
	}
	
	/**
	 * Get all combinations of the ribosomal substrates (i.e., amino acids)
	 * associated with a list of domains.
	 * 
	 * @param domains
	 *            list of ribosomal domains
	 * @return all substrate combinations
	 */
	public static List<List<Substrate>> getSubstrateCombinations(List<Domain> domains) {
		List<List<Substrate>> input = new ArrayList<List<Substrate>>();
		for (Domain domain : domains)
			input.add(domain.substrates());
		List<List<Substrate>> output =  Combinations.getCombinations(input);
		System.out.println("[RibosomalClusterAnalyzer] Generated " + output.size() 
				+ " ribosomal substrate combinations");
		return output;
	}
	
	/**
	 * Get a list of ribosomal domains from a propeptide.
	 * 
	 * @param propeptide
	 *            sequence of the cleaved propeptide
	 * @param orf
	 *            parent open reading frame
	 * @return a list of ribosomal domains
	 * @throws RibosomalScaffoldGenerationException
	 */
	public static List<Domain> getDomainsFromPropeptide(String propeptide, Orf orf) 
			throws RibosomalScaffoldGenerationException {
		List<Domain> domains = new ArrayList<Domain>();
		for (int i = 0; i < propeptide.length(); i++) {
			String name = propeptide.charAt(i) + "_" + (i+1);
			Domain domain = new Domain(i, i, 0.0d, name);
			domain.setType(GenericDomains.AMINO_ACID);
			if (propeptide.charAt(i) == '*') {
				String codon = getCodon(i, propeptide, orf);
				List<Codons> codons = getPossibleCodons(codon);
				for (Codons c : codons) {
					Substrate substrate = new Substrate(c.getAminoAcid());
					domain.addSubstrate(substrate);
				}
			} else {
				for (ProteinogenicAminoAcids aa : ProteinogenicAminoAcids.values()) {
					if (aa.abbreviation().equals(propeptide.charAt(i) + "")) {
						Substrate substrate = new Substrate(aa);
						domain.addSubstrate(substrate);
					} 
				}
				if (domain.substrates().size() == 0)
					throw new RibosomalScaffoldGenerationException("Could not get substrate for " + domain.name());
			}
			domains.add(domain);
		}
		System.out.println("[RibosomalClusterAnalyzer] Generated " + domains.size() 
				+ " domains from propeptide " + propeptide);
		return domains;
	}
	
	/**
	 * Get the DNA codon from an amino acid of the propeptide, used when
	 * the DNA codon contains an 'N' (unknown nucleotide).  
	 * @param i				index of 
	 * @param propeptide	the cleaved propeptide sequence 
	 * @param orf			the parent open reading frame
	 * @return				the DNA codon which corresponds to this character
	 * @throws RibosomalScaffoldGenerationException
	 */
	public static String getCodon(int i, String propeptide, Orf orf) 
			throws RibosomalScaffoldGenerationException {
		String dnaSequence = orf.dnaSequence();
		String aaSequence = orf.sequence();
		int aaStart = aaSequence.indexOf(propeptide);
		if (aaStart == -1) 
			throw new RibosomalScaffoldGenerationException("Error: could not get propeptide "
					+ "sequence from parent open reading frame " + orf.name());
		int dnaStart = aaStart * 3;
		String codon = dnaSequence.substring(dnaStart, dnaStart + 3);
		System.out.println("Got codon " + codon + " for open reading frame " + orf.name());
		return codon;
	}

	/**
	 * Get a list of all possible codons that could correspond to a codon in a
	 * user-input sequence where one or more nucleotides is 'N' (unknown).
	 * 
	 * @param codon	N-containing codon detected in a user-input sequence 
	 * @return		all possible codons that could correspond to the detected codon 
	 * @throws RibosomalScaffoldGenerationException
	 */
	public static List<Codons> getPossibleCodons(String codon) 
			throws RibosomalScaffoldGenerationException {
		if (codon.length() != 3)
			throw new RibosomalScaffoldGenerationException("Error: could not get possible "
					+ "codons for string " + codon + " with length != 3");
		List<Codons> codons = new ArrayList<Codons>();
		char c1 = codon.toLowerCase().charAt(0);
		char c2 = codon.toLowerCase().charAt(1);
		char c3 = codon.toLowerCase().charAt(2);
		for (Codons c : Codons.values()) {
			String comparison = c.toString().toLowerCase();
			if (comparison.charAt(0) == c1 || c1 == 'n'
					&& comparison.charAt(1) == c2 || c2 == 'n'
					&& comparison.charAt(2) == c3 || c3 == 'n')
				codons.add(c);
		}
		System.out.println("Found " + codons.size() + " possible codons "
				+ "for input codon " + codon);
		return codons;
	}

}
