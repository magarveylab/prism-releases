package ca.mcmaster.magarveylab.prism.cluster.reactions.impl;

import java.util.List;

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.AtomType;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Substrates;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

public class TryptophanDioxygenaseReaction extends GenericReaction implements Reaction {

	public TryptophanDioxygenaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { TailoringDomains.TRYPTOPHAN_DIOXYGENASE };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
	
		Module module = plan.get(0);
		Residue residue = scaffold.residue(module);
		IAtom ketone = residue.ketone();

		if (ketone != null) {
			IAtom gammaCarbon = Substrates.getPeptideSp2GammaCarbon(ketone, molecule);
			List<IBond> bonds = molecule.getConnectedBondsList(gammaCarbon);
			for (IBond bond : bonds) {
				IAtom connectedAtom = bond.getConnectedAtom(gammaCarbon);
				if (connectedAtom.getSymbol().equals("C")
						&& connectedAtom.getHybridization() == IAtomType.Hybridization.SP2
						&& connectedAtom.getImplicitHydrogenCount() == 1) {
					List<IBond> connectedAtomBonds = molecule.getConnectedBondsList(connectedAtom);

					// remove nitrogen bond
					for (IBond connectedAtomBond : connectedAtomBonds) {
						IAtom nitrogen = connectedAtomBond.getConnectedAtom(connectedAtom);
						if (nitrogen.getSymbol().equals("N")) {
							molecule.removeBond(connectedAtomBond);
							nitrogen.setImplicitHydrogenCount(2);
							AtomTypeManipulator.configure(nitrogen, new AtomType("N"));
							IBond otherNitrogenBond = molecule.getConnectedBondsList(nitrogen).get(0);
							otherNitrogenBond.setOrder(IBond.Order.SINGLE);
						}
					}
					
					// remove aromatic carbon
					molecule.removeBond(bond);
					molecule.removeAtom(connectedAtom);
				} else if (connectedAtom.getSymbol().equals("C") 
						&& connectedAtom.getHybridization() == IAtomType.Hybridization.SP2
						&& connectedAtom.getImplicitHydrogenCount() == 0) {
					// create single bond
					bond.setOrder(IBond.Order.SINGLE);

					// create ketone and add
					IAtomContainer newKetone = new AtomContainer();
					IAtom oxygen = new Atom("O");
					newKetone.addAtom(oxygen);
					molecule.add(newKetone);
					UtilityReactions.addBond(gammaCarbon, oxygen, molecule, IBond.Order.DOUBLE);
				}
			}
		} else {
			throw new ScaffoldGenerationException("Error: tryptophan dioxygenase site is null");
		}
	}
	
}
