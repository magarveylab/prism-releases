package ca.mcmaster.magarveylab.prism.util.fragment;

import java.util.ArrayList;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IMoleculeSet;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;

/**
 * 
 * @author Lian
 */
public class NRPFragmenter extends NRPAnalyser {

	private ArrayList<IMolecule> break1BFragment, break1YFragment, break1CircleFragment;
	private ArrayList<IMolecule> break1BranchCFragment, break2Fragment;

	public NRPFragmenter(IMolecule iMol) {
		super(iMol);
		generatePeaks();
	}

	public ArrayList<IMolecule> getAllFragments() {
		ArrayList<IMolecule> fragmentList = new ArrayList<IMolecule>();
		fragmentList.addAll(break1BFragment);
		fragmentList.addAll(break1YFragment);
		fragmentList.addAll(break1BranchCFragment);
		fragmentList.addAll(break1CircleFragment);
		fragmentList.addAll(break2Fragment);
		return fragmentList;
		// return duplicateFilter(fragmentList);
	}

	public ArrayList<IMolecule> getAll1BreakFragments() {
		ArrayList<IMolecule> fragmentList = new ArrayList<IMolecule>();
		fragmentList.addAll(break1BFragment);
		fragmentList.addAll(break1YFragment);
		fragmentList.addAll(break1BranchCFragment);
		fragmentList.addAll(break1CircleFragment);
		return fragmentList;
		// return duplicateFilter(fragmentList);
	}

	public ArrayList<IMolecule> get1BreakBFragments() {
		return break1BFragment;
		// return duplicateFilter(break1BFragment);
	}

	public ArrayList<IMolecule> get1BreakYFragments() {
		return break1YFragment;
		// return duplicateFilter(break1YFragment);
	}

	public ArrayList<IMolecule> get1BreakCircleFragments() {
		return break1CircleFragment;
		// return duplicateFilter(break1CircleFragment);
	}

	public ArrayList<IMolecule> get1BreakBranchCFragments() {
		return break1BranchCFragment;
		// return duplicateFilter(break1BranchCFragment);
	}

	public ArrayList<IMolecule> get2BreakFragments() {
		return break2Fragment;
		// return duplicateFilter(break2Fragment);
	}

	private void generatePeaks() {
		generateBreak1Peaks();
		generateBreak2Peaks();
	}

	private void generateBreak1Peaks() {
		break1BFragment = new ArrayList<IMolecule>();
		break1YFragment = new ArrayList<IMolecule>();
		break1CircleFragment = new ArrayList<IMolecule>();
		break1BranchCFragment = new ArrayList<IMolecule>();

		// one amide bond in backbone breaks
		for (int i = 0; i < atomCount; i++) {
			if (atom_isNonAlphaC[i]) { // break at C-N
				IAtom atom = molecule.getAtom(i);
				IBond bond = atom_next_backbone_bond[i];
				if (bond != null) {
					molecule.removeBond(bond);

					IMoleculeSet mol_set = ConnectivityChecker.partitionIntoMolecules(molecule);

					// one bond in a ring breaks, molecule breaks into one piece
					if (mol_set.getAtomContainerCount() == 1) {
						IMolecule mol = mol_set.getMolecule(0);
						break1CircleFragment.add(mol);
					}

					// one bond in non-cyclic backbone breaks, molecule breaks
					// into two pieces
					if (mol_set.getAtomContainerCount() == 2) {
						for (int j = 0; j < mol_set.getAtomContainerCount(); j++) {
							IMolecule mol = mol_set.getMolecule(j);
							if (mol.contains(atom)) {
								// b-ion
								break1BFragment.add(mol);
							} else {
								// y-ion
								break1YFragment.add(mol);
							}

						}
					}
					molecule.addBond(bond);
				}
			}
		}
		// one C-C bond in branches breaks
		for (int i = 0; i < bondCount; i++) {
			if (bond_is_branchCCbond[i]) {
				IBond bond = molecule.getBond(i);
				molecule.removeBond(bond);

				IMoleculeSet mol_set = ConnectivityChecker.partitionIntoMolecules(molecule);
				// molecule breaks into one piece
				if (mol_set.getAtomContainerCount() == 1) {
					IMolecule mol = mol_set.getMolecule(0);
					break1BranchCFragment.add(mol);
				}

				// molecule breaks into two pieces
				if (mol_set.getAtomContainerCount() == 2) {
					for (int j = 0; j < mol_set.getAtomContainerCount(); j++) {
						IMolecule mol = mol_set.getMolecule(j);
						break1BranchCFragment.add(mol);
					}
				}
				molecule.addBond(bond);
			}
		}
	}

	private void generateBreak2Peaks() {
		break2Fragment = new ArrayList<IMolecule>();

		// two bonds break
		for (int i = 0; i < atomCount - 1; i++)
			for (int j = i + 1; j < atomCount; j++) {
				if (atom_isNonAlphaC[i] && atom_isNonAlphaC[j]) { // break at
																	// C-N
					IBond bond1 = atom_next_backbone_bond[i];
					IBond bond2 = atom_next_backbone_bond[j];
					if (bond1 != null && bond2 != null) {
						IAtom atomC1 = molecule.getAtom(i);
						IAtom atomC2 = molecule.getAtom(j);

						IAtom atomN1 = bond1.getConnectedAtom(atomC1);
						IAtom atomN2 = bond2.getConnectedAtom(atomC2);

						molecule.removeBond(bond1);
						molecule.removeBond(bond2);

						IMoleculeSet mol_set = ConnectivityChecker.partitionIntoMolecules(molecule);
						for (int k = 0; k < mol_set.getAtomContainerCount(); k++) {
							IMolecule mol = mol_set.getMolecule(k);
							int countC = 0, countN = 0;
							if (mol.contains(atomC1))
								countC++;
							if (mol.contains(atomC2))
								countC++;
							if (mol.contains(atomN1))
								countN++;
							if (mol.contains(atomN2))
								countN++;
							// exclude y-ion or b-ion in Break1Peaks
							if (countN + countC == 1)
								continue;
							break2Fragment.add(mol);
						}

						molecule.addBond(bond1);
						molecule.addBond(bond2);
					}
				}
			}
	}
}
