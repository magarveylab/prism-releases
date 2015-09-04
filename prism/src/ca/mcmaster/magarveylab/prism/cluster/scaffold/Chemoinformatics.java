package ca.mcmaster.magarveylab.prism.cluster.scaffold;

import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;

/**
 * Performs chemoinformatic analysis of virtual molecules.
 * @author skinnider
 */
public class Chemoinformatics {

	/**
	 * Performs chemoinformatic analysis of virtual adenylation or acyltransferase substrates.
	 * @author skinnider
	 *
	 */
	public static class Substrates {

		/**
		 * Get the alpha carbon of a peptidic substrate, given the nitrogen or carbon terminus. 
		 * @param terminus	the nitrogen or carbon terminus
		 * @param mol		molecule (substrate) to search
		 * @return			the alpha carbon
		 */
		public static IAtom getPeptideAlphaCarbon(IAtom terminus, IAtomContainer mol) {
			IAtom alphaCarbon = null;
			List<IBond> nitrogenBonds = mol.getConnectedBondsList(terminus);
			for (IBond bond : nitrogenBonds) {
				IAtom connectedAtom = bond.getConnectedAtom(terminus);
				if (connectedAtom.getSymbol().equals("C") && connectedAtom.getHybridization() == IAtomType.Hybridization.SP3)
					alphaCarbon = connectedAtom;
			}
			return alphaCarbon;
		}

		/**
		 * Get a SP2-hybridized gamma carbon of a peptide substrate, given either the N or C terminus of the molecule.
		 * @param terminus	the nitrogen or carbon terminus
		 * @param mol		molecule to search
		 * @return			the SP2 gamma carbon, if present
		 */
		public static IAtom getPeptideSp2GammaCarbon(IAtom terminus, IAtomContainer mol) {
			IAtom betaCarbon = null;
			IAtom alphaCarbon = getPeptideAlphaCarbon(terminus, mol);
			List<IBond> bonds = mol.getConnectedBondsList(alphaCarbon);
			for (IBond bond : bonds) {
				IAtom connectedAtom = bond.getConnectedAtom(alphaCarbon);
				if (connectedAtom.getSymbol().equals("C")) {
					List<IBond> connectedAtomBonds = mol.getConnectedBondsList(connectedAtom);
					for (IBond connectedAtomBond : connectedAtomBonds) {
						IAtom connectedAtom2 = connectedAtomBond.getConnectedAtom(connectedAtom);
						if (connectedAtom2.getSymbol().equals("C")
								&& connectedAtom2.getHybridization() == IAtomType.Hybridization.SP2)
							betaCarbon = connectedAtom2;
					}
				}
			}
			return betaCarbon;
		}

		/**
		 * Get the atom at which a substrate reacts to extend a growing chain. Uses iodine to locate extender atom.
		 * @param mol	the molecule to search
		 * @return		the starter atom
		 */
		public static IAtom getStarterAtom(IAtomContainer mol) {
			IAtom iodine = Atoms.getIodine(mol);
			IBond iodineBond = Bonds.getIodineBond(mol);
			IAtom start = iodineBond.getConnectedAtom(iodine);
			return start;
		}
		
		/**
		 * Get the atom at which a substrate is extended. Uses fluorine to locate the extension atom.
		 * @param mol	the molecule to search
		 * @return		the extending atom
		 */
		public static IAtom getExtenderAtom(IAtomContainer mol) {
			IAtom fluorine = Atoms.getFluorine(mol);
			IAtom end = null;
			if (fluorine != null) {
				IBond fluorineBond = Bonds.getFluorineBond(mol);
				end = fluorineBond.getConnectedAtom(fluorine);
			}
			return end;
		}

		/**
		 * Get the ketone carbon of a substrate, given the location of the alpha carbon.
		 * @param end	the alpha carbon
		 * @param mol	the substrate
		 * @return
		 */
		public static IAtom getKetoneFromAlphaCarbon(IAtom alphaCarbon, IAtomContainer mol) {
			IAtom ketoneCarbon = null;
			List<IBond> bonds = mol.getConnectedBondsList(alphaCarbon);
			for (IBond bond : bonds) {
				IAtom atom = bond.getConnectedAtom(alphaCarbon);
				if (atom.getAtomTypeName().equals("C.sp2"))
					ketoneCarbon = atom; 
			}
			return ketoneCarbon;
		}

		/**
		 * Get the alpha carbon of a substrate, given the location of the ketone carbon.
		 * @param end	the ketone carbon
		 * @param mol	the substrate
		 * @return
		 */
		public static IAtom getAlphaCarbonFromKetone(IAtom ketone, IAtomContainer mol) {
			IAtom alphaCarbon = null;
			List<IBond> bonds = mol.getConnectedBondsList(ketone);
			for (IBond bond : bonds) {
				IAtom atom = bond.getConnectedAtom(ketone);
				if (atom.getSymbol().equals("C"))
					alphaCarbon = atom; 
			}
			return alphaCarbon;
		}

	}
	
	/**
	 * Performs chemoinformatic analysis of virtual molecule atoms.
	 * @author skinnider
	 *
	 */
	public static class Atoms {

		/**
		 * Get the sulfur atom within a molecule.
		 * @param molecule	molecule to search
		 * @return			the (last) sulfur atom
		 */	
		public static IAtom getSulfur(IAtomContainer molecule) {
			IAtom sulfur = null;
			for (IAtom atom : molecule.atoms()) {
				if (atom.getSymbol().equals("S"))
					sulfur = atom;
			}
			return sulfur;
		}

		/**
		 * Get the oxygen atom within a molecule.
		 * @param molecule	molecule to search
		 * @return			the (last) oxygen atom
		 */	
		public static IAtom getOxygen(IAtomContainer molecule) {
			IAtom oxygen = null;
			for (IAtom atom : molecule.atoms()) {
				if (atom.getSymbol().equals("O"))
					oxygen = atom;
			}
			return oxygen;
		}

		/**
		 * Get the fluorine atom in a molecule.
		 * @param molecule	the molecule to search
		 * @return			the (last) fluorine atom
		 */
		public static IAtom getFluorine(IAtomContainer molecule) {
			IAtom fluorine = null;
			for (IAtom atom : molecule.atoms()) {
				if (atom.getSymbol().equals("F"))
					fluorine = atom;
			}
			return fluorine;
		}

		/**
		 * Get the iodine atom in a molecule.
		 * @param mol	the molecule to search
		 * @return		the (last) iodine atom 
		 */
		public static IAtom getIodine(IAtomContainer mol) {
			IAtom iodine = null;
			for (IAtom atom : mol.atoms()) {
				if (atom.getSymbol().equals("I"))
					iodine = atom;
			}
			return iodine;
		}
		
		/**
		 * Get the nitrogen atom in a molecule.
		 * @param mol	the molecule to search
		 * @return		the (last) nitrogen atom 
		 */
		public static IAtom getNitrogen(IAtomContainer mol) {
			IAtom nitrogen = null;
			for (IAtom atom : mol.atoms()) {
				if (atom.getSymbol().equals("N"))
					nitrogen = atom;
			}
			return nitrogen;
		}

		/**
		 * Get the methylene carbon in a molecule, defined as a carbon with two
		 * hydrogens.
		 * 
		 * @param molecule
		 *            the molecule to search
		 * @return the (last) methylene carbon
		 */
		public static IAtom getMethyleneCarbon(IAtomContainer molecule) {
			IAtom carbon = null;
			for (IAtom atom : molecule.atoms())
				if (atom.getSymbol().equals("C")
						&& atom.getImplicitHydrogenCount() == 2)
					carbon = atom;
			return carbon;
		}

		/**
		 * Get the methyl carbon in a molecule, defined as a carbon with three
		 * hydrogens.
		 * 
		 * @param molecule
		 *            the molecule to search
		 * @return the (last) methyl carbon
		 */
		public static IAtom getMethylCarbon(IAtomContainer molecule) {
			IAtom carbon = null;
			for (IAtom atom : molecule.atoms())
				if (atom.getSymbol().equals("C")
						&& atom.getImplicitHydrogenCount() == 3)
					carbon = atom;
			return carbon;
		}

		/**
		 * Get the hydroxyl group within the substructure of a molecule.
		 * 
		 * @param substructure
		 *            substructure to restrict search to
		 * @param molecule
		 *            the molecule to search
		 * @return the (last) hydroxyl oxygen
		 */
		public static IAtom getHydroxyl(IAtomContainer substructure,
				IAtomContainer molecule) {
			IAtom oxygen = null;
			for (IAtom atom : substructure.atoms()) {
				if (atom.getSymbol().equals("O")
						&& molecule.getConnectedBondsCount(atom) == 1
						&& molecule.getBondOrderSum(atom) == 1)
					oxygen = atom;
			}
			return oxygen;
		}

		/**
		 * Get the hydroxyl group in a molecular structure.
		 * 
		 * @param molecule
		 *            the molecule to search
		 * @return the (last) hydroxyl oxygen
		 */
		public static IAtom getHydroxyl(IAtomContainer molecule) {
			IAtom hydroxyl = null;
			for (IAtom atom : molecule.atoms())
				if (atom.getSymbol().equals("O")
						&& molecule.getBondOrderSum(atom) == 1)
					hydroxyl = atom;
			return hydroxyl;
		}
		
		/**
		 * Get all hydroxyl oxygens in a molecule.
		 * 
		 * @param molecule
		 *            the molecule to search
		 * @return all hydroxyl oxygens
		 */
		public static List<IAtom> getAllHydroxyls(IAtomContainer molecule) {
			List<IAtom> oxygens = new ArrayList<IAtom>();
			for (IAtom atom : molecule.atoms())
				if (atom.getSymbol().equals("O")
						&& molecule.getBondOrderSum(atom) == 1)
					oxygens.add(atom);
			return oxygens;
		}

		/**
		 * Get the (last) primary, non-amide amino group within a substructure of a molecule.
		 * @param structure		substructure to restrict search to
		 * @param molecule		the molecule to search
		 * @return				the (last) primary amine 
		 */
		public static IAtom getPrimaryAmine(IAtomContainer structure, IAtomContainer molecule) {
			IAtom nitrogen = null;
			for (IAtom atom : structure.atoms()) {
				if (!atom.getSymbol().equals("N")) //nitrogen
					continue;
				if (molecule.getBondOrderSum(atom) > 1) //one bond
					continue;
				List<IAtom> atoms = molecule.getConnectedAtomsList(atom);
				if (atoms.size() > 0) {
					IAtom atom2 = atoms.get(0); // not amide check
					for (IAtom atom3 : structure.getConnectedAtomsList(atom2))
						if (atom3.getSymbol().equals("O")
								&& structure.getBond(atom2, atom3).getOrder() == IBond.Order.DOUBLE)
							continue;
					nitrogen = atom;
				}
			}
			return nitrogen;
		}
		
		/**
		 * Get the hydroxyl group connected to a (reduced) ketone.
		 * @param ketone	ketone carbon
		 * @param molecule	molecule to search
		 * @return			connected hydroxyl
		 */
		public static IAtom getConnectedHydroxyl(IAtom ketone, IAtomContainer molecule) {
			IAtom oxygen = null;
			for (IAtom atom : molecule.getConnectedAtomsList(ketone)) 
				if (atom.getSymbol().equals("O") && molecule.getBond(ketone, atom).getOrder() == IBond.Order.SINGLE
						&& molecule.getConnectedBondsCount(atom) < 2) {
					System.out.println("Atom: " + atom.getSymbol() + ", bonds: " + molecule.getConnectedBondsCount(atom));
					oxygen = atom;
				}
			return oxygen;
		}
		
		/**
		 * Get the (last) carbon atom connected to an atom in question.
		 * @param atom		atom in question
		 * @param molecule	molecule to search
		 * @return			connected carbon
		 */
		public static IAtom getConnectedCarbon(IAtom atom, IAtomContainer molecule) {
			IBond bond = Bonds.getConnectedCarbonBond(molecule, atom);
			if (bond != null) {
				return bond.getConnectedAtom(atom);
			} else {
				System.out.println("Error: no connected carbon");
				return null;
			}
		}
	
		/**
		 * Get the (last) nitrogen atom connected to an atom in question.
		 * @param atom		atom in question
		 * @param molecule	molecule to search
		 * @return			connected nitrogen
		 */
		public static IAtom getConnectedNitrogen(IAtom atom, IAtomContainer molecule) {
			IBond bond = Bonds.getConnectedNitrogenBond(molecule, atom);
			return bond.getConnectedAtom(atom);
		}
		
		/**
		 * Get the (last) sulfur atom connected to an atom in question.
		 * @param atom		atom in question
		 * @param molecule	molecule to search
		 * @return			connected sulfur
		 */
		public static IAtom getConnectedSulfur(IAtom atom, IAtomContainer molecule) {
			IBond bond = Bonds.getConnectedBond(molecule, atom, "S");
			return bond.getConnectedAtom(atom);
		}

		/**
		 * Get the (last) oxygen atom connected to an atom in question.
		 * @param atom		atom in question
		 * @param molecule	molecule to search
		 * @return			connected oxygen
		 */
		public static IAtom getConnectedOxygen(IAtom atom, IAtomContainer molecule) {
			IBond bond = Bonds.getConnectedBond(molecule, atom, "O");
			return bond.getConnectedAtom(atom);
		}

		/**
		 * Find all carboxyl carbons in a molecule.
		 * @param molecule	the molecule to search
		 * @return			a list of all carboxyl carbons
		 */
		public static List<IAtom> getAllCarboxyls(IAtomContainer molecule) {
			List<IAtom> carboxyls = new ArrayList<IAtom>();
			for (int i = 0; i < molecule.getAtomCount(); i++) {
				boolean hasCarbonyl = false, hasAlcohol = false;
				IAtom atom = molecule.getAtom(i);
				if (!atom.getSymbol().equals("C")) continue;
				List<IBond> bonds = molecule.getConnectedBondsList(atom);
				for (IBond bond : bonds) {
					IAtom connectedAtom = bond.getConnectedAtom(atom);
					if (connectedAtom.getSymbol().equals("O") && bond.getOrder() == IBond.Order.DOUBLE) 
						hasCarbonyl = true;
					if (connectedAtom.getSymbol().equals("O") && bond.getOrder() == IBond.Order.SINGLE 
							&& connectedAtom.getImplicitHydrogenCount() == 1)
						hasAlcohol = true;
				}
				if (hasCarbonyl && hasAlcohol) {
					IAtom carboxyl = molecule.getAtom(i);
					carboxyls.add(carboxyl);
				}
			}
			return carboxyls;
		}

		/**
		 * Find all ketones in a molecule.
		 * @param molecule	the molecule to search
		 * @return		a list of all ketone carbons
		 */
		public static List<IAtom> getAllKetones(IAtomContainer molecule) {
			List<IAtom> ketones = new ArrayList<IAtom>();
			boolean hasCarbonyl = false;
			for (int i = 0; i < molecule.getAtomCount(); i++) {
				IAtom atom = molecule.getAtom(i);
				if (!atom.getSymbol().equals("C")) continue;
				List<IBond> bonds = molecule.getConnectedBondsList(atom);
				for (IBond bond : bonds) {
					IAtom connectedAtom = bond.getConnectedAtom(atom);
					if (connectedAtom.getSymbol().equals("O") && bond.getOrder() == IBond.Order.DOUBLE) 
						hasCarbonyl = true;
				}
				if (hasCarbonyl) ketones.add(atom);
			}
			return ketones;
		}
		
		/**
		 * Get the carbon atom at which histidine is chlorinated on the peptide antibiotic GE81112.
		 * @param molecule	molecule to search
		 * @return			the chlorinated histidine carbon atom
		 */
		public static IAtom getHistidineChlorinationAtom(IAtomContainer molecule) {
			IAtom chlorination = null;
			for (IAtom atom : molecule.atoms()) {
				if (atom.getSymbol().equals("N") && atom.getImplicitHydrogenCount() == 1) {
					List<IAtom> atoms = molecule.getConnectedAtomsList(atom);
					for (IAtom atom2 : atoms) 
						if (atom2.getSymbol().equals("C") && molecule.getConnectedBondsCount(atom2) == 2)
							chlorination = atom2;
				}
			}
			return chlorination;
		}
		
		/**
		 * Get the carbon atom at which tryptophan is chlorinated on natural products such as rebeccamycin and cladoniamide.
		 * @param molecule	molecule to search
		 * @return			the chlorinated tryptophan carbon atom
		 */
		public static IAtom getTryptophanChlorinationAtom(IAtomContainer molecule) {
			IAtom chlorination = null;
			for (IAtom atom : molecule.atoms()) {
				if (atom.getSymbol().equals("N")) {
					List<IAtom> atoms2 = molecule.getConnectedAtomsList(atom);
					for (IAtom atom2 : atoms2) {
						if (atom2.getSymbol().equals("C") 
								&& molecule.getConnectedBondsCount(atom2) == 3) {
							List<IAtom> atoms3 = molecule.getConnectedAtomsList(atom);
							for (IAtom atom3 : atoms3)
								if (atom3.getSymbol().equals("C") 
										&& molecule.getConnectedBondsCount(atom3) == 2)
									chlorination = atom3;
						}
					}
				}
			}
			return chlorination;
		}
		
		/**
		 * Get the carbon atom at which proline is chlorinated on natural products such as hormaomycin,
		 * pyoluteorin, pyralomycin and pyrrolomicin. 
		 * @param molecule	molecule to search		
		 * @return			the chlorinated proline carbon atom
		 */
		public static IAtom getProlineChlorinationAtom(IAtomContainer molecule) {
			IAtom chlorination = null;
			for (IAtom atom : molecule.atoms())
				if (atom.getSymbol().equals("C") && molecule.getConnectedBondsCount(atom) == 2)
					chlorination = atom;
			return chlorination;
		}
		
		/**
		 * Get the meta carbon of an aromatic amino acid.
		 * @param molecule	molecule (AA) to search
		 * @return			the carbon at the meta position of this amino acid
		 */
		public static IAtom getMetaCarbon(IAtomContainer molecule) {
			IAtom ortho = null, meta = null;
			atomLoop:
			for (IAtom atom : molecule.atoms()) {
				if (atom.getSymbol().equals("C") && molecule.getConnectedBondsCount(atom) == 3
						&& atom.getHybridization() == IAtomType.Hybridization.SP2) {
					for (IAtom atom2 : molecule.getConnectedAtomsList(atom)) {
						// if bond to -OH, this is not the R carbon
						if (!atom2.getSymbol().equals("C"))
							continue atomLoop;
						// get adjacent carbon with only 2 bonds
						if (atom2.getSymbol().equals("C") && molecule.getConnectedAtomsCount(atom2) == 2
								&& atom2.getHybridization() == IAtomType.Hybridization.SP2)
							ortho = atom2;
					}
					if (ortho != null) {
						for (IAtom atom3 : molecule.getConnectedAtomsList(ortho)) {
							if (!atom3.getSymbol().equals("C") || atom3 == atom)
								continue;
							if (molecule.getConnectedAtomsCount(atom3) == 2 
									&& atom3.getHybridization() == IAtomType.Hybridization.SP2)
								meta = atom3;
						}
						
					}
				}
			}
			return meta;
		}

		/**
		 * Get an aromatic carbon with only two bonds. 
		 * @param molecule	molecule to search
		 * @return			an aromatic carbon with two bonds
		 */
		public static IAtom getAromaticCarbon(IAtomContainer molecule) {
			IAtom para = null;
			for (IAtom atom : molecule.atoms()) {
				if (atom.getSymbol().equals("C") && molecule.getConnectedBondsCount(atom) == 2
						&& atom.getHybridization() == IAtomType.Hybridization.SP2) 
					para = atom;
			}
			return para;
		}
		
		/**
		 * Get the aromatic carbon single-bonded to a hydroxyl.
		 * @param molecule	molecule to search
		 * @return			aromatic hydroxyl carbon
		 */
		public static IAtom getAromaticHydroxylCarbon(IAtomContainer molecule) {
			IAtom carbon = null;
			for (IAtom atom : molecule.atoms()) {
				if (atom.getSymbol().equals("O")) {
					List<IBond> bonds = molecule.getConnectedBondsList(atom);
					for (IBond bond : bonds) {
						if (bond.getOrder() == IBond.Order.SINGLE) {
							IAtom atom2 = bond.getConnectedAtom(atom);
							if (atom2.getSymbol().equals("C") && atom2.getHybridization() == IAtomType.Hybridization.SP2)
								carbon = atom2;
						}
					}
				}
			}
			return carbon;
		}

		/**
		 * Get the hydroxyl group from a carboxylic acid, given the carboxyl
		 * carbon.
		 * 
		 * @param cTerminus
		 *            the carboxyl carbon
		 * @param molecule
		 *            the molecule to search
		 * @return
		 */
		public static IAtom getCarboxylAlcohol(IAtom cTerminus,
				IAtomContainer molecule) {
			IAtom alcohol = null;
			List<IBond> cTerminusBonds = molecule
					.getConnectedBondsList(cTerminus);
			for (IBond bond : cTerminusBonds) {
				IAtom connectedAtom = bond.getConnectedAtom(cTerminus);
				if (connectedAtom.getSymbol().equals("O")
						&& bond.getOrder() == IBond.Order.SINGLE
						&& connectedAtom.getImplicitHydrogenCount() == 1)
					alcohol = connectedAtom;
			}
			return alcohol;
		}
		
	}
	
	/**
	 * Performs chemoinformatic analysis of virtual molecule bonds.
	 * @author skinnider
	 *
	 */
	public static class Bonds {

		/**
		 * Get the bond between and atom and a connected atom of a given type.
		 * @param molecule	molecule to search	
		 * @param atom		atom in question
		 * @param symbol	type of atom to search for
		 * @return			the bond between the given atom and an atom of the given type
		 */
		public static IBond getConnectedBond(IAtomContainer molecule, IAtom atom, String symbol) {
			IBond connectedBond = null;
			List<IBond> bonds = molecule.getConnectedBondsList(atom);
			for (IBond bond : bonds) {
				IAtom connectedAtom = bond.getConnectedAtom(atom);
				if (connectedAtom.getSymbol().equals(symbol))
					connectedBond = bond;
			}
			if (connectedBond == null)
				System.out.println("[SubstrateUtil] Error: no connected " + symbol + " bond!");
			return connectedBond;
		}

		/**
		 * Get the last bond between an atom and any other atom it is connected to. 
		 * @param molecule	molecule to search	
		 * @param atom		atom in question
		 * @return			the bond between the given atom and any other atom
		 */
		public static IBond getConnectedBond(IAtomContainer molecule, IAtom atom) {
			IBond connectedBond = null;
			List<IBond> bonds = molecule.getConnectedBondsList(atom);
			for (IBond bond : bonds) {
				IAtom connectedAtom = bond.getConnectedAtom(atom);
				if (connectedAtom != null)
					connectedBond = bond;
			}
			if (connectedBond == null)
				System.out.println("[SubstrateUtil] Error: no connected bond!");
			return connectedBond;
		}

		/**
		 * Get the bond between an atom and a connected oxygen atom
		 * @param molecule	molecule to search
		 * @param atom		atom in question
		 * @return			the bond between the given atom and an oxygen atom
		 */
		public static IBond getConnectedOxygenBond(IAtomContainer molecule, IAtom atom) {
			String oxygen = "O";
			return getConnectedBond(molecule, atom, oxygen);
		}

		/**
		 * Get the bond between an atom and a connected carbon atom.
		 * @param molecule	molecule to search
		 * @param atom		atom in question
		 * @return			the bond between the given atom and a carbon atom
		 */
		public static IBond getConnectedCarbonBond(IAtomContainer molecule, IAtom atom) {
			String carbon = "C";
			return getConnectedBond(molecule, atom, carbon);
		}

		/**
		 * Get the oxygen atom connected to a given atom.
		 * @param molecule	molecule to search
		 * @param atom		atom in question
		 * @return			oxygen atom connected to the given atom
		 */
		public static IAtom getConnectedOxygen(IAtomContainer molecule, IAtom atom) {
			List<IBond> bonds = molecule.getConnectedBondsList(atom);
			IAtom oxygen = null;
			for (IBond bond : bonds) {
				IAtom connectedAtom = bond.getConnectedAtom(atom);
				if (connectedAtom.getSymbol().equals("O"))
					oxygen = connectedAtom;
			}
			if (oxygen == null)
				System.out.println("[SubstrateUtil] Error: no connected oxygen!");
			return oxygen;
		
		}

		/**
		 * Get the bond between an atom and a connected nitrogen atom.	
		 * @param molecule	molecule to search
		 * @param atom		atom in question
		 * @return			the bond between the given atom and a nitrogen atom
		 */
		public static IBond getConnectedNitrogenBond(IAtomContainer molecule, IAtom atom) {
			String nitrogen = "N";
			return getConnectedBond(molecule, atom, nitrogen);
		}

		/**
		 * Get the fluorine-carbon bond used to indicate the point at which the polyketide unit is extended.
		 * @param molecule	the molecule to search
		 * @return			the fluorine-carbon bond
		 */
		public static IBond getFluorineBond(IAtomContainer molecule) {
			IAtom fluorine = Atoms.getFluorine(molecule);
			return getConnectedBond(molecule, fluorine);
		}

		/**
		 * Get the iodine-carbon bond used to indicate the point at which the polyketide unit extends a growing chain.
		 * @param molecule	the molecule to search
		 * @return			the iodine-carbon bond
		 */
		public static IBond getIodineBond(IAtomContainer molecule) {
			IAtom iodine = Atoms.getIodine(molecule);
			return getConnectedCarbonBond(molecule, iodine);
		}

		/**
		 * Get the bond between a carboxyl carbon and its hydroxyl group, given the carboxyl carbon. 
		 * @param molecule		the molecule to search
		 * @param cTerminus		the carboxyl carbon
		 * @return
		 */
		public static IBond getCarboxylAlcoholBond(IAtomContainer molecule, IAtom cTerminus) {
			IBond alcoholBond = null;
			List<IBond> cTerminusBonds = molecule.getConnectedBondsList(cTerminus);
			for (IBond bond : cTerminusBonds) {
				IAtom connectedAtom = bond.getConnectedAtom(cTerminus);
				if (connectedAtom.getSymbol().equals("O") && bond.getOrder() == IBond.Order.SINGLE 
						&& connectedAtom.getImplicitHydrogenCount() == 1)
					alcoholBond = bond;
			}
			return alcoholBond;
		}
		
	}
	
}
