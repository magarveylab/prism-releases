package ca.mcmaster.magarveylab.prism.util.fragment;


import java.util.Arrays;
import java.util.List;

import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;

/**
 * 
 * @author Lian
 * @version 0.1, 12/26/2009
 */
public class NRPAnalyser {
	protected IMolecule molecule;
	protected int atomCount;
	protected int bondCount;
	protected boolean atom_isN[], atom_isNonAlphaC[], atom_isAlphaC[];
	protected boolean bond_is_branchCCbond[];

	protected IBond atom_next_backbone_bond[];
	protected boolean atom_is_backbone[];

	public NRPAnalyser(IMolecule iMol) {
		this.molecule = iMol;
		atomCount = iMol.getAtomCount();
		bondCount = iMol.getBondCount();
		atom_isN = new boolean[atomCount];
		atom_isNonAlphaC = new boolean[atomCount];
		atom_isAlphaC = new boolean[atomCount];
		atom_is_backbone = new boolean[atomCount];
		atom_next_backbone_bond = new IBond[atomCount];
		bond_is_branchCCbond = new boolean[bondCount];

		Arrays.fill(atom_isN, false);
		Arrays.fill(atom_isNonAlphaC, false);
		Arrays.fill(atom_isAlphaC, false);
		Arrays.fill(atom_is_backbone, false);
		Arrays.fill(atom_next_backbone_bond, null);
		Arrays.fill(bond_is_branchCCbond, false);

		findBackboneC();
		findBackboneAlphaC();
		findBackboneN();
		findBranchCCBond();
	}

	private void findBackboneN() {
		for (int i = 0; i < atomCount; i++) {
			IAtom atom = molecule.getAtom(i);

			/* check if the type of the atom is N.sp3 */
			if (!atom.getAtomTypeName().equals("N.sp3") && !atom.getAtomTypeName().equals("N.amide")) {
				continue;
			}
			/* the number of hydrogen >= 1 */
			/*
			 * if (atom.getHydrogenCount() < 1) { continue; }
			 */
			List<IBond> bond_list = molecule.getConnectedBondsList(atom);
			boolean is_used[] = new boolean[bond_list.size()];
			Arrays.fill(is_used, false);
			/* check if there is a single bond to a C atom */
			int pnt = -1;
			/*
			 * if (atom.getHydrogenCount() == 1) { for (int j = 0; j <
			 * bond_list.size(); j++) { IBond bond = bond_list.get(j); if
			 * (bond.getOrder() == IBond.Order.SINGLE &&
			 * ((NRPAtom)bond.getConnectedAtom(atom)).getIsC()) { pnt = j;
			 * break; }
			 * 
			 * } if (pnt == -1) { continue; } is_used[pnt] = true; }
			 */

			/* check if there is a single bond to C alpha */
			pnt = -1;
			for (int j = 0; j < bond_list.size(); j++) {
				if (is_used[j]) {
					continue;
				}
				IBond bond = bond_list.get(j);
				if (bond.getOrder() == IBond.Order.SINGLE
						&& atom_isAlphaC[molecule.getAtomNumber(bond.getConnectedAtom(atom))] == true) {
					pnt = j;
					break;
				}
			}
			if (pnt == -1) {
				continue;
			}
			is_used[pnt] = true;
			atom_isN[i] = true;
			atom_next_backbone_bond[i] = bond_list.get(pnt);
		}
	}

	private void findBackboneAlphaC() {
		for (int i = 0; i < atomCount; i++) {
			IAtom atom = molecule.getAtom(i);

			// System.out.print("\n" + atom.getAtomTypeName());
			/* check if the type of the atom is C.sp3 */
			if (!atom.getAtomTypeName().equals("C.sp3")) {
				continue;
			}
			/*
			 * the number of hydrogen ?? some non-genomic amino acids do not a
			 * hydrogen bond
			 */
			/*
			 * if (atom.getHydrogenCount() < 1) { continue; }
			 */
			List<IBond> bond_list = molecule.getConnectedBondsList(atom);
			boolean is_used[] = new boolean[bond_list.size()];
			Arrays.fill(is_used, false);
			/* check if there is a single bond to an nitrogen */
			int pnt = -1;
			for (int j = 0; j < bond_list.size(); j++) {
				IBond bond = bond_list.get(j);
				if (bond.getOrder() == IBond.Order.SINGLE
						&& (bond.getConnectedAtom(atom).getAtomTypeName().equals("N.sp3") || bond
								.getConnectedAtom(atom).getAtomTypeName().equals("N.amide"))) {
					pnt = j;
					break;
				}

			}
			if (pnt == -1) {
				continue;
			}
			is_used[pnt] = true;
			/* check if there is a single bond to a C atom */
			pnt = -1;
			for (int j = 0; j < bond_list.size(); j++) {
				if (is_used[j]) {
					continue;
				}

				IBond bond = bond_list.get(j);
				if (bond.getOrder() == IBond.Order.SINGLE
						&& atom_isNonAlphaC[molecule.getAtomNumber(bond.getConnectedAtom(atom))] == true) {
					pnt = j;
					break;
				}

			}
			if (pnt == -1) {
				continue;
			}
			is_used[pnt] = true;
			atom_isAlphaC[i] = true;
			atom_next_backbone_bond[i] = bond_list.get(pnt);
		}
	}

	private void findBackboneC() {
		for (int i = 0; i < atomCount; i++) {
			IAtom atom = molecule.getAtom(i);

			/* check if the type of the atom is C.sp2 */
			if (!atom.getAtomTypeName().equals("C.sp2")) {
				continue;
			}
			/* the number of hydrogen is 0 */
			if (atom.getImplicitHydrogenCount() != 0) {
				continue;
			}
			/* check if there is three bonds */
			List<IBond> bond_list = molecule.getConnectedBondsList(atom);
			if (bond_list.size() != 3) {
				continue;
			}
			boolean is_used[] = new boolean[3];
			Arrays.fill(is_used, false);
			/* check if there is a double bond with oxygen */
			int pnt = -1;
			for (int j = 0; j < bond_list.size(); j++) {
				IBond bond = bond_list.get(j);
				if (bond.getOrder() == IBond.Order.DOUBLE
						&& bond.getConnectedAtom(atom).getAtomTypeName().equals("O.sp2")) {
					pnt = j;
					break;
				}

			}
			if (pnt == -1) {
				continue;
			}
			is_used[pnt] = true;
			/* check if there is a single bond with C-alpha */
			pnt = -1;
			for (int j = 0; j < bond_list.size(); j++) {
				if (is_used[j]) {
					continue;
				}
				IBond bond = bond_list.get(j);
				if (bond.getOrder() == IBond.Order.SINGLE
						&& bond.getConnectedAtom(atom).getAtomTypeName().equals("C.sp3")) {
					pnt = j;
					break;
				}

			}
			if (pnt == -1) {
				continue;
			}
			is_used[pnt] = true;

			/* check if there is a single bond with N */
			pnt = -1;
			for (int j = 0; j < bond_list.size(); j++) {
				if (is_used[j]) {
					continue;
				}
				IBond bond = bond_list.get(j);
				IAtom nghb = bond.getConnectedAtom(atom);
				if (bond.getOrder() == IBond.Order.SINGLE && nghb.getAtomTypeName().equals("N.amide")) {
					pnt = j;
					break;
				}
				/*
				 * if (bond.getOrder() == IBond.Order.SINGLE &&
				 * nghb.getAtomTypeName().equals("O.sp3") &&
				 * nghb.getHydrogenCount() == 1) { pnt = j; break; }
				 */
			}
			/*
			 * if (pnt == -1) { continue; } is_used[pnt] = true;
			 */
			atom_isNonAlphaC[i] = true;
			if (pnt >= 0) {
				atom_next_backbone_bond[i] = bond_list.get(pnt);
			}
		}
	}

	private void findBranchCCBond() {
		for (int i = 0; i < bondCount; i++) {
			IBond bond = molecule.getBond(i);

			// only SINGLE bonds are considered
			// DOUBLE bonds are filtered out as they are hard to break
			if (bond.getOrder() != IBond.Order.SINGLE)
				continue;
			// bond should connect two atoms
			if (bond.getAtomCount() != 2)
				continue;

			IAtom atom1 = bond.getAtom(0);
			IAtom atom2 = bond.getAtom(1);

			// check if both atoms in the bond are C
			if (!atom1.getAtomTypeName().equals("C.sp2") && !atom1.getAtomTypeName().equals("C.sp3"))
				continue;
			if (!atom2.getAtomTypeName().equals("C.sp2") && !atom2.getAtomTypeName().equals("C.sp3"))
				continue;

			// check this C-C bond is not in peptide backbone
			if (atom_isNonAlphaC[molecule.getAtomNumber(atom1)] && atom_isAlphaC[molecule.getAtomNumber(atom2)])
				continue;
			if (atom_isNonAlphaC[molecule.getAtomNumber(atom2)] && atom_isAlphaC[molecule.getAtomNumber(atom1)])
				continue;

			bond_is_branchCCbond[i] = true;
		}
	}

}
