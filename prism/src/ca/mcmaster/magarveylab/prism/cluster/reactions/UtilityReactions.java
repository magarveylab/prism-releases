package ca.mcmaster.magarveylab.prism.cluster.reactions;

import java.util.List;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Atoms;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Bonds;
import ca.mcmaster.magarveylab.prism.util.exception.BondFormationException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.SmilesIO;

/**
 * Perform common virtual reactions on chemical structures.
 * @author skinnider
 *
 */
public class UtilityReactions {

	/**
	 * Create a single bond between two atoms. 
	 * @param atom1		the first atom
	 * @param atom2		the second atom
	 * @param container	the container containing both atoms
	 * @throws ScaffoldGenerationException 
	 */
	public static void addBond(IAtom atom1, IAtom atom2, IAtomContainer container) 
			throws BondFormationException {
		int atom1Idx = container.getAtomNumber(atom1);
		int atom2Idx = container.getAtomNumber(atom2);
		if (atom1Idx == -1)
			throw new BondFormationException("Could not form single bond: could not locate atom 1!");
		if (atom2Idx == -1)
			throw new BondFormationException("Could not form single bond: could not locate atom 2!");
		if (container.getBond(atom1, atom2) != null)
			throw new BondFormationException("Could not add bond: bond already exists!");
		container.addBond(atom1Idx, atom2Idx, IBond.Order.SINGLE);
	}
	
	/**
	 * Create a bond with a specified order between two atoms.
	 * @param atom1		the first atom
	 * @param atom2		the second atom
	 * @param molecule	the container containing both atoms
	 * @param order		the order
	 * @throws ScaffoldGenerationException 
	 */
	public static void addBond(IAtom atom1, IAtom atom2, IAtomContainer molecule, IBond.Order order) 
			throws BondFormationException {
		int atom1Idx = molecule.getAtomNumber(atom1);
		int atom2Idx = molecule.getAtomNumber(atom2);
		if (atom1Idx == -1)
			throw new BondFormationException("Could not form " + order + " bond: could not locate atom 1!");
		if (atom2Idx == -1)
			throw new BondFormationException("Could not form " + order + " bond: could not locate atom 2!");
		if (molecule.getBond(atom1, atom2) != null)
			throw new BondFormationException("Could not add bond: bond already exists!");
		molecule.addBond(atom1Idx, atom2Idx, order);
	}

	/**
	 * Remove a bond between two atoms.
	 * 
	 * @param atom1
	 *            the first atom
	 * @param atom2
	 *            the second atom
	 * @param molecule
	 *            the molecule containing both atoms
	 * @throws ScaffoldGenerationException
	 */
	public static void removeBond(IAtom atom1, IAtom atom2,
			IAtomContainer molecule) throws ScaffoldGenerationException {
		IBond bond = molecule.getBond(atom1, atom2);
		if (bond == null)
			throw new ScaffoldGenerationException("Could not remove bond: bond is null");
		molecule.removeBond(bond);
	}
	
	/**
	 * Set the order of the bond between two atoms.
	 * 
	 * @param atom1
	 *            the first atom
	 * @param atom2
	 *            the second atom
	 * @param molecule
	 *            the molecule containing both atoms
	 * @param order
	 *            the order
	 * @throws ScaffoldGenerationException
	 */
	public static void setBondOrder(IAtom atom1, IAtom atom2,
			IAtomContainer molecule, IBond.Order order)
			throws ArrayIndexOutOfBoundsException, ScaffoldGenerationException {
		IBond bond = molecule.getBond(atom1, atom2);
		if (bond == null)
			throw new ScaffoldGenerationException("Could not set bond order: bond is null");
		bond.setOrder(order);
	}

	/**
	 * Remove a fluorine atom from this container.
	 * @param container	the container to remove fluorine from
	 */
	public static void removeFluorine(IAtomContainer container) {
		IAtom fluorine = Atoms.getFluorine(container);
		if (fluorine != null) {
			IBond starterFluorineBond = Bonds.getFluorineBond(container);
			container.removeBond(starterFluorineBond);
			container.removeAtom(fluorine);
		}
	}
	
	/**
	 * Remove an iodine atom from this container.
	 * @param container	the container to remove iodine from
	 */
	public static void removeIodine(IAtomContainer container) {
		IAtom iodine = Atoms.getIodine(container);
		if (iodine != null) {
			IBond extenderIodineBond = Bonds.getIodineBond(container);
			container.removeBond(extenderIodineBond);
			container.removeAtom(iodine);
		}
	}
	
	/**
	 * Remove the hydroxyl group from a carboxylic acid to form an aldehyde (and
	 * automatically partition).
	 * 
	 * @param cTerminus
	 *            the carboxylic acid carbon
	 * @param container
	 *            the container to remove carboxylic -OH from
	 */
	public static void removeCarboxylAlcohol(IAtom cTerminus, IAtomContainer container) {
		IBond carboxylAlcoholBond = Bonds.getCarboxylAlcoholBond(container, cTerminus);
		IAtom carboxylAlcohol = Atoms.getCarboxylAlcohol(cTerminus, container);
		if (carboxylAlcoholBond == null)
			System.out.println("Error: could not detect carboxyl alcohol bond!");
		if (carboxylAlcohol == null)
			System.out.println("Error: could not detect carboxyl alcohol!");
		container.removeBond(carboxylAlcoholBond);
		container.removeAtom(carboxylAlcohol);
	}
	
	/**
	 * Reduce a double-bonded oxygen to a single-bonded oxygen.
	 * @param carbon		carbon atom
	 * @param molecule		parent molecule 
	 */
	public static void reduceKetone(IAtom carbon, IAtomContainer molecule) {
		List<IBond> bonds = molecule.getConnectedBondsList(carbon);
		for (IBond bond : bonds) {
			IAtom atom = bond.getConnectedAtom(carbon);
			if (atom.getSymbol().equals("O") && bond.getOrder() == IBond.Order.DOUBLE) {
				bond.setOrder(IBond.Order.SINGLE);
				break;
			}
		}
	}
	
	/**
	 * Oxidize a single-bonded oxygen to a double-bonded oxygen.
	 * @param carbon		carbon atom
	 * @param molecule		parent molecule 
	 */
	public static void oxidizeAlcohol(IAtom carbon, IAtomContainer molecule) {
		List<IBond> bonds = molecule.getConnectedBondsList(carbon);
		for (IBond bond : bonds) {
			IAtom atom = bond.getConnectedAtom(carbon);
			if (atom.getSymbol().equals("O") && bond.getOrder() == IBond.Order.SINGLE) {
				bond.setOrder(IBond.Order.DOUBLE);
				break;
			}
		}
	}

	/**
	 * Remove a single double-bonded oxygen from a carbon atom.
	 * @param carbon		carbon atom
	 * @param molecule		parent molecule
	 */
	public static void removeKetone(IAtom carbon, IAtomContainer molecule) {
		List<IBond> bonds = molecule.getConnectedBondsList(carbon);
		for (IBond bond : bonds) {
			IAtom atom = bond.getConnectedAtom(carbon);
			if (atom.getSymbol().equals("O") && bond.getOrder() == IBond.Order.DOUBLE) {
				molecule.removeBond(bond);
				molecule.removeAtom(atom);
				break;
			}
		}	
	}
	
	/**
	 * Remove a single-bonded oxygen from a carbon atom.
	 * @param carbon		carbon atom
	 * @param molecule		parent molecule 
	 */
	public static void removeAlcohol(IAtom carbon, IAtomContainer molecule) {
		List<IBond> bonds = molecule.getConnectedBondsList(carbon);
		for (IBond bond : bonds) {
			IAtom atom = bond.getConnectedAtom(carbon);
			if (atom.getSymbol().equals("O") && bond.getOrder() == IBond.Order.SINGLE) {
				molecule.removeBond(bond);
				molecule.removeAtom(atom);
				break;
			}
		}
	}
	
	/**
	 * Remove an oxygen from a carbon atom.
	 * @param carbon	carbon atom
	 * @param molecule	parent molecule
	 */
	public static void removeOxygen(IAtom carbon, IAtomContainer molecule) {
		List<IBond> bonds = molecule.getConnectedBondsList(carbon);
		for (IBond bond : bonds) {
			IAtom atom = bond.getConnectedAtom(carbon);
			if (atom.getSymbol().equals("O")) {
				molecule.removeBond(bond);
				molecule.removeAtom(atom);
				break;
			}
		}
	}
	
	/**
	 * Remove an atom, and all bonds connected to it, from a molecule.
	 * @param atom		atom to remove
	 * @param molecule	parent molecule 
	 */
	public static void removeAtom(IAtom atom, IAtomContainer molecule) {
		List<IBond> bonds = molecule.getConnectedBondsList(atom);
		for (IBond bond : bonds) 
			molecule.removeBond(bond);
		molecule.removeAtom(atom);
	}
	
	/**
	 * Perform the virtual decarboxylation of a chemical scaffold.
	 * @param atom1		the carboxyl group carbon
	 * @param atom2		the atom connected to the carboxyl group 
	 * @param molecule	parent molecule
	 * @throws ScaffoldGenerationException
	 */
	public static void decarboxylate(IAtom atom1, IAtom atom2,
			IAtomContainer molecule) throws ScaffoldGenerationException {
		IBond c1c2Bond = molecule.getBond(atom1, atom2);
		if (c1c2Bond == null)
			throw new ScaffoldGenerationException("Error: could not locate bond between atoms 1 and 2 in decarboxylation!");
		molecule.removeBond(c1c2Bond);
		for (IBond bond : molecule.getConnectedBondsList(atom1)) {
			IAtom atom = bond.getConnectedAtom(atom1);
			molecule.removeBond(bond);
			molecule.removeAtom(atom);
		}
		molecule.removeAtom(atom1);
	}
	
	/**
	 * Add a functional group to a chemical scaffold.
	 * @param smiles	SMILES of the functional group, with the attachment site labelled with an iodine 
	 * @param atom		atom to functionalize
	 * @param molecule	parent molecule 
	 * @throws ScaffoldGenerationException 
	 * @throws BondFormationException 
	 */
	public static void functionalize(String smiles, IAtom atom, IAtomContainer molecule) 
			throws ScaffoldGenerationException, BondFormationException {
		// create functional group
		IAtomContainer functionalGroup = null;
		try {
			functionalGroup = SmilesIO.molecule(smiles);
		} catch (Exception e) {
			throw new ScaffoldGenerationException("Error encountered building functional group SMILES");
		}
		
		// get reverse-prenylation site
		IAtom iodine = Atoms.getIodine(functionalGroup);
		IAtom carbon = Atoms.getConnectedCarbon(iodine, functionalGroup);
		UtilityReactions.removeAtom(iodine, functionalGroup);
		
		// get bond order sum
		int bonds = -1;
		if (atom.getSymbol().equals("N")) {
			bonds = 3;
		} else if (atom.getSymbol().equals("C")) {
			bonds = 4;
		} else if (atom.getSymbol().equals("O")
				|| atom.getSymbol().equals("S")) {
			bonds = 2;
		}
		if (bonds == -1)
			throw new ScaffoldGenerationException("Error: tried to functionalize"
					+ " atom that is not C, N, S, or O!");

		if (molecule.getConnectedBondsCount(atom) < bonds
				&& molecule.getBondOrderSum(atom) < bonds) {
			// add bond
			molecule.add(functionalGroup);
			UtilityReactions.addBond(carbon, atom, molecule);
		} else {
			throw new ScaffoldGenerationException("Error: tried to functionalize"
					+ " atom with too many bonds!");
		}
	}
	
}
