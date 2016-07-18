package ca.mcmaster.magarveylab.prism.cluster.scaffold;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.enums.CStarterSubstrates;
import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.enums.interfaces.SubstrateType;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;
import ca.mcmaster.magarveylab.enums.substrates.TypeIIPolyketideStarters;
import ca.mcmaster.magarveylab.prism.cluster.analysis.DomainAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.analysis.TypeIIPolyketideAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Atoms;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Substrates;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.enums.hmms.AcyltransferaseHmms;
import ca.mcmaster.magarveylab.prism.util.SmilesIO;
import ca.mcmaster.magarveylab.prism.util.exception.BondFormationException;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;


/**
 * Generate scaffold residues.
 * @author skinnider
 *
 */
public class ResidueGenerator {
	
	/**
	 * Generate a residue for use in scaffold generation from a biosynthetic module. 
	 * @param module	module to use for residue creation
	 * @return			the generated residue
	 * @throws IOException 
	 * @throws NoResidueException 
	 * @throws ScaffoldGenerationException 
	 * @throws CDKException 
	 */
	public static Residue residue(Module module) throws IOException, NoResidueException, 
			ScaffoldGenerationException, CDKException {
		Residue residue = null;
		if (module.type() == ModuleTypes.ADENYLATION || module.type() == ModuleTypes.TRANS_ADENYLATION
				|| module.type() == ModuleTypes.ACYL_ADENYLATE) {
			residue = adenylationResidue(module);
		} else if (module.type() == ModuleTypes.ACYLTRANSFERASE) {
			residue = acyltransferaseResidue(module);
		} else if (module.type() == ModuleTypes.TRANS_ADENYLATION_INSERTION) {
			residue = transAdenylationResidue(module);
		} else if (module.type() == ModuleTypes.C_STARTER) {
			residue = cStarterResidue(module);
		} else if (module.type() == ModuleTypes.TYPE_II_PKS) {
			residue = typeIIPolyketideResidue(module);
		} else if (module.type() == ModuleTypes.TYPE_II_PKS_STARTER) {
			residue = typeIIPolyketideStarterResidue(module);
		} else if (module.type() == ModuleTypes.RIBOSOMAL) { 
			residue = ribosomalResidue(module);
		}
		return residue;
	}
	
	private static Residue adenylationResidue(Module module) throws IOException, NoResidueException, CDKException {
		// get substrate molecule
		Domain domain = module.scaffold();
		if (domain == null || domain.topSubstrate() == null)
			throw new NoResidueException("Error: could not generate residue for " 
					+ module.type() + " module with no scaffold!");
		SubstrateType type = domain.topSubstrate().type();
		String smiles = type.smiles();
		IAtomContainer molecule = SmilesIO.molecule(smiles);
		
		// get ketone and alpha carbon
		IAtom nitrogen = Substrates.getExtenderAtom(molecule);
		IAtom ketone = Substrates.getStarterAtom(molecule);
		IAtom alphaCarbon = Substrates.getAlphaCarbonFromKetone(ketone, molecule);
		
		// remove halogens
		UtilityReactions.removeIodine(molecule);
		UtilityReactions.removeFluorine(molecule);

		// create residue
		Residue residue = new Residue(module);
		residue.setNitrogen(nitrogen);
		residue.setKetone(ketone);
		residue.setAlphaCarbon(alphaCarbon);
		residue.setStructure(molecule);
		
		return residue;
	}
	
	private static Residue ribosomalResidue(Module module) throws IOException,
			NoResidueException, CDKException {
		String aa = module.first().name().substring(0, 1).toLowerCase();
		String smiles = null;
		for (ProteinogenicAminoAcids a : ProteinogenicAminoAcids.values()) {
			if (a.abbreviation().toLowerCase().equals(aa)) {
				smiles = a.smiles();
			}
		}
				
		IAtomContainer molecule = SmilesIO.molecule(smiles);
		
		// get ketone and alpha carbon
		IAtom nitrogen = Substrates.getExtenderAtom(molecule);
		IAtom ketone = Substrates.getStarterAtom(molecule);
		IAtom alphaCarbon = Substrates.getAlphaCarbonFromKetone(ketone, molecule);
		
		// remove halogens
		UtilityReactions.removeIodine(molecule);
		UtilityReactions.removeFluorine(molecule);

		// create residue
		Residue residue = new Residue(module);
		residue.setNitrogen(nitrogen);
		residue.setKetone(ketone);
		residue.setAlphaCarbon(alphaCarbon);
		residue.setStructure(molecule);
		
		return residue;
	}
	
	private static Residue acyltransferaseResidue(Module module) throws IOException, CDKException {
		// get substrate molecule
		Domain domain = module.scaffold();
		SubstrateType type = domain.topSubstrate().type();
		String smiles = type.smiles();
		IAtomContainer molecule = SmilesIO.molecule(smiles);
				
		// get ketone and alpha carbon
		IAtom alphaCarbon = Substrates.getExtenderAtom(molecule);
		IAtom ketone = Substrates.getStarterAtom(molecule);
		
		// remove halogens
		UtilityReactions.removeIodine(molecule);
		UtilityReactions.removeFluorine(molecule);
		
		// create residue
		Residue residue = new Residue(module);
		residue.setKetone(ketone);
		residue.setAlphaCarbon(alphaCarbon);
		residue.setStructure(molecule);
		
		return residue;
	}
	
	private static Residue transAdenylationResidue(Module module) throws IOException,
	 		ScaffoldGenerationException, NoResidueException, BondFormationException, CDKException {
		Residue residue = null;
		List<Domain> adenylation = new ArrayList<Domain>();
		for (Domain domain : module.domains()) 
			if (domain.type() == ThiotemplatedDomains.ADENYLATION) 
				adenylation.add(domain);
		
		if (adenylation.size() == 0) {
			System.out.println("[ResidueGenerator] Error: no adenylation domain in trans-adenylation insertion module!");
		} else if (adenylation.size() == 1) {
			residue = adenylationResidue(module);
		} else if (adenylation.size() >= 2) {
			// only consider 2 adenylation domains; start by getting molecules
			Domain starterDomain = adenylation.get(0);
			Domain extenderDomain = adenylation.get(1);
			SubstrateType starterType = starterDomain.topSubstrate().type();
			SubstrateType extenderType = extenderDomain.topSubstrate().type();
			String starterSmiles = starterType.smiles();
			String extenderSmiles = extenderType.smiles();
			IAtomContainer starter = SmilesIO.molecule(starterSmiles);
			IAtomContainer extender = SmilesIO.molecule(extenderSmiles);
			
			// get nitrogen, ketone, and alpha carbon
			IAtom starterKetone = Substrates.getStarterAtom(starter);
			IAtom starterNitrogen = Substrates.getExtenderAtom(starter);
			IAtom extenderKetone = Substrates.getStarterAtom(extender);
			IAtom extenderNitrogen = Substrates.getExtenderAtom(extender);
			IAtom extenderAlphaCarbon = Substrates.getAlphaCarbonFromKetone(extenderKetone, extender);

			// remove halogens
			UtilityReactions.removeIodine(starter);
			UtilityReactions.removeFluorine(starter);
			UtilityReactions.removeIodine(extender);
			UtilityReactions.removeFluorine(extender);

			// create didomain residue
			IAtomContainer molecule = new AtomContainer();
			molecule.add(starter);
			molecule.add(extender);
			UtilityReactions.addBond(starterKetone, extenderNitrogen, molecule);
			
			// create residue
			residue = new Residue(module);
			residue.setNitrogen(starterNitrogen);
			residue.setKetone(extenderKetone);
			residue.setAlphaCarbon(extenderAlphaCarbon);
			residue.setStructure(molecule);
		}
		
		return residue;
	}
	
	private static Residue cStarterResidue(Module module) throws IOException, CDKException {
		// get substrate molecule
		Domain domain = module.domains().get(0);
		CStarterSubstrates type = DomainAnalyzer.cStarterType(domain);
		String smiles = type.smiles();
		IAtomContainer molecule = SmilesIO.molecule(smiles);
		
		// get ketone and remove iodine
		IAtom ketone = Substrates.getStarterAtom(molecule);
		UtilityReactions.removeIodine(molecule);
		
		// get alpha carbon
		IAtom alphaCarbon = Atoms.getConnectedCarbon(ketone, molecule);
		
		// create residue
		Residue residue = new Residue(module);
		residue.setKetone(ketone);
		residue.setAlphaCarbon(alphaCarbon);
		residue.setStructure(molecule);

		return residue;
	}
	
	private static Residue typeIIPolyketideResidue(Module module) throws IOException, CDKException {
		// get default acetate molecule		
		SubstrateType type = AcyltransferaseHmms.MALONYL_COA_1;
		String smiles = type.smiles();
		IAtomContainer molecule = SmilesIO.molecule(smiles);
				
		// get ketone and alpha carbon
		IAtom alphaCarbon = Substrates.getExtenderAtom(molecule);
		IAtom ketone = Substrates.getStarterAtom(molecule);
		
		// remove halogens
		UtilityReactions.removeIodine(molecule);
		UtilityReactions.removeFluorine(molecule);
		
		// create residue
		Residue residue = new Residue(module);
		residue.setKetone(ketone);
		residue.setAlphaCarbon(alphaCarbon);
		residue.setStructure(molecule);
		
		return residue;
	}
	
	private static Residue typeIIPolyketideStarterResidue(Module module) 
			throws CDKException, IOException, NoResidueException {
		TypeIIPolyketideStarters starter = TypeIIPolyketideAnalyzer.getStarter(module);
		if (starter == null)
			throw new NoResidueException("Could not get starter unit for type II polyketide starter residue!");
		
		String smiles = starter.smiles();
		IAtomContainer molecule = SmilesIO.molecule(smiles);
				
		// get ketone
		IAtom ketone = Substrates.getStarterAtom(molecule);
		UtilityReactions.removeIodine(molecule);
		
		// get alpha carbon
		IAtom alphaCarbon = Atoms.getConnectedCarbon(ketone, molecule);
		
		// create residue
		Residue residue = new Residue(module);
		residue.setKetone(ketone);
		residue.setAlphaCarbon(alphaCarbon);
		residue.setStructure(molecule);
		
		return residue;
	}
	
}
