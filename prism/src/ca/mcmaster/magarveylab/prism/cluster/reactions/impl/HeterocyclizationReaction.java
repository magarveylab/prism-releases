package ca.mcmaster.magarveylab.prism.cluster.reactions.impl;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Atoms;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Bonds;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;
import ca.mcmaster.magarveylab.wasp.util.SmilesIO;

/**
 * Execute the heterocyclization of a serine, threonine, or cysteine residue to
 * an oxazoline or thiazoline, respectively.
 * 
 * @author skinnider
 *
 */
public class HeterocyclizationReaction extends GenericReaction implements
		Reaction {

	public HeterocyclizationReaction(ReactionPlan plan, Scaffold scaffold,
			Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { ThiotemplatedDomains.CONDENSATION };
	}

	public void execute() throws NoResidueException,
			TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		
		Module last = plan.get(0);
		Module module = plan.get(1);
		Residue residue = scaffold.residue(module);
		Residue lastResidue = scaffold.residue(last);
		if (residue == null || lastResidue == null)
			throw new ScaffoldGenerationException(
					"Error: could not heterocyclize on null residue!");

		IAtom ketone = residue.ketone();
		if (ketone == null)
			throw new ScaffoldGenerationException("Could not find ketone for ser/thr/cys residue");
		IAtom lastKetone = lastResidue.ketone();
		if (lastKetone == null)
			throw new ScaffoldGenerationException("Could not find ketone for ser/thr/cys-adjacent residue");
		IAtom alphaCarbon = residue.alphaCarbon();
		if (alphaCarbon == null)
			throw new ScaffoldGenerationException("Could not find alpha carbon for ser/thr/cys residue");
		IAtom betaCarbon = null;
		IAtom terminalAtom = null;
		
		for (IBond bond : molecule.getConnectedBondsList(alphaCarbon)) {
			IAtom connectedAtom = bond.getConnectedAtom(alphaCarbon);
			if (connectedAtom.getSymbol().equals("C") && connectedAtom != ketone)
				betaCarbon = connectedAtom;
		}
		
		if (betaCarbon == null)
			throw new ScaffoldGenerationException("Could not generate thiazoline/oxazoline:"
					+ " could not find beta carbon!");
		
		for (IBond bond : molecule.getConnectedBondsList(betaCarbon)) {
			IAtom connectedAtom = bond.getConnectedAtom(betaCarbon);
			if (connectedAtom.getSymbol().equals("O") 
					|| connectedAtom.getSymbol().equals("S") 
					&& connectedAtom != alphaCarbon)
				terminalAtom = connectedAtom;
		}
		
		if (terminalAtom == null)
			throw new ScaffoldGenerationException("Could not generate thiazoline/oxazoline: "
					+ "could not find terminal oxygen or sulfur atom!");

		// create bond between terminal O/S and ketone carbon
		UtilityReactions.addBond(terminalAtom, lastKetone, molecule);

		// remove ketone oxygen
		IAtom ketoneOxygen = Atoms.getConnectedOxygen(lastKetone, molecule);
		IBond oxygenBond = Bonds.getConnectedOxygenBond(molecule, lastKetone);
		molecule.removeBond(oxygenBond);
		molecule.removeAtom(ketoneOxygen);
		
		// create double bond between ketone carbon and nitrogen
		IBond nitrogenBond = Bonds.getConnectedNitrogenBond(molecule, lastKetone);
		nitrogenBond.setOrder(IBond.Order.DOUBLE);
		
		System.out.println("After azoline reaction: " + SmilesIO.smiles(molecule));
	}

}
