package ca.mcmaster.magarveylab.prism.cluster.reactions;

import java.io.IOException;

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.ResidueGenerator;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.SmilesIO;

/**
 * Executes common virtual reactions in natural product scaffold biosynthesis.
 * @author skinnider
 *
 */
public class ScaffoldReactions {
	
	/**
	 * Create a new natural product scaffold. 
	 * @param starter	the first scaffold module
	 * @return			the new scaffold
	 * @throws IOException 
	 * @throws NoResidueException 
	 * @throws ScaffoldGenerationException 
	 * @throws CDKException 
	 */
	public static Scaffold createScaffold(Module module)
			throws IOException, NoResidueException, ScaffoldGenerationException, CDKException {
		if (!module.isActive())
			return null;
		
		Scaffold scaffold = new Scaffold();
		IAtomContainer molecule = new AtomContainer();
		Residue residue = ResidueGenerator.residue(module);
		
		if (residue == null) {
			String s = module.scaffold() ==  null ? "with no scaffold" : module.scaffold().name(); 
			System.out.println("Could not generate residue for " + module.type() + " module " + s);
			return null;
		}
		
		molecule.add(residue.structure());
		
		scaffold.setMolecule(molecule);
		scaffold.addResidue(module, residue);
		System.out.println("[ScaffoldReactions] Started scaffold with SMILES " + 
				SmilesIO.smiles(scaffold.molecule()));
		return scaffold;
	}
	
	/**
	 * Extend a growing natural product scaffold.
	 * @param module	the module to extend with
	 * @param scaffold	the current scaffold
	 * @throws IOException 
	 * @throws NoResidueException 
	 * @throws ScaffoldGenerationException 
	 * @throws CDKException 
	 */
	public static void extendScaffold(Module module, Scaffold scaffold) 
			throws ScaffoldGenerationException, IOException, NoResidueException, CDKException {
		if (!module.isActive() || !module.canExtend())
			return;

		IAtomContainer molecule = scaffold.molecule();
		Residue starter = scaffold.getLastResidue();
		Residue extender = ResidueGenerator.residue(module);
		if (extender == null) {
			String s = module.scaffold() ==  null ? "with no scaffold" : module.scaffold().name(); 
			System.out.println("Could not generate residue for " + module.type() + " module " + s);
			return;
		}

		molecule.add(extender.structure());

		Module extenderModule = extender.module();
		IAtom extenderAtom = (extenderModule.isAdenylationModule() || extenderModule
				.type() == ModuleTypes.RIBOSOMAL) ? extender.nitrogen() : extender.alphaCarbon();
		UtilityReactions.addBond(starter.ketone(), extenderAtom, molecule);
		
		scaffold.addResidue(module, extender);
		System.out.println("[ScaffoldReactions] Extended scaffold with SMILES " + 
				SmilesIO.smiles(scaffold.molecule()));				
	}
		
	/**
	 * Finish a scaffold by adding a hydroxyl to its terminal ketone, creating a carboxylic acid.
	 * @param scaffold	the scaffold to finish
	 * @throws ScaffoldGenerationException 
	 * @throws CDKException 
	 */
	public static void finishScaffold(Scaffold scaffold) throws ScaffoldGenerationException, CDKException {
		Residue last = scaffold.getLastResidue();
		IAtomContainer molecule = scaffold.molecule();
		IAtom ketone = last.ketone();
		
		// create hydroxyl group
		IAtomContainer hydroxyl = new AtomContainer();
		IAtom oxygen = new Atom("O");
		oxygen.setImplicitHydrogenCount(1);
		hydroxyl.addAtom(oxygen);
	
		// add to scaffold & add bond
		molecule.add(hydroxyl);
		last.structure().add(hydroxyl);
		UtilityReactions.addBond(oxygen, ketone, molecule);
		scaffold.setMolecule(molecule);
		
		System.out.println("[ScaffoldReactions] Finished scaffold with SMILES " + 
				SmilesIO.smiles(scaffold.molecule()));				
	}
	
}
