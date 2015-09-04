package ca.mcmaster.magarveylab.prism.cluster.reactions.impl;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Atoms;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

public class ProlineDehydrogenaseReaction extends GenericReaction implements Reaction {
	
	public ProlineDehydrogenaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { TailoringDomains.PROLINE_DEHYDROGENASE };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
	
		Module module = plan.get(0);
		Residue residue = scaffold.residue(module);
		IAtom ketone = residue.ketone();
		IAtom nitrogen = residue.nitrogen();
		if (nitrogen == null) 
			nitrogen = Atoms.getNitrogen(residue.structure());

		if (ketone == null)
			throw new ScaffoldGenerationException("Error: could not find ketone for proline dehydrogenase reaction");
		if (nitrogen == null)
			throw new ScaffoldGenerationException("Error: could not find nitrogen for proline dehydrogenase reaction");
		
		// get alpha carbon
		IAtom alphaCarbon = null;
		for (IAtom atom : molecule.getConnectedAtomsList(ketone))
			if (atom.getSymbol().equals("C") && molecule.getConnectedAtomsCount(atom) == 3)
				for (IAtom atom2 : molecule.getConnectedAtomsList(atom))
					if (atom2.getSymbol().equals("N"))
						alphaCarbon = atom;

		// get beta carbon
		IAtom betaCarbon = null;
		for (IAtom atom : molecule.getConnectedAtomsList(alphaCarbon)) 
			if (atom.getSymbol().equals("C") && atom != ketone)
				betaCarbon = atom;
		
		// create double bond b/w beta and alpha carbon
		if (betaCarbon == null) {
			System.out.println("[TailoringReactions] Unable to locate proline beta carbon");
			return;
		} else {
			IBond bond = molecule.getBond(alphaCarbon, betaCarbon);
			bond.setOrder(IBond.Order.DOUBLE);
		}
		
		// get gamma carbon
		IAtom gammaCarbon = null;
		for (IAtom atom : molecule.getConnectedAtomsList(betaCarbon))
			if (atom.getSymbol().equals("C") && atom != alphaCarbon)
				gammaCarbon = atom;

		// get delta carbon
		IAtom deltaCarbon = null;
		for (IAtom atom : molecule.getConnectedAtomsList(gammaCarbon))
			if (atom.getSymbol().equals("C") && atom != betaCarbon)
				deltaCarbon = atom;
		
		// add double bond
		if (deltaCarbon == null || gammaCarbon == null) {
			System.out.println("[TailoringReactions] Unable to locate proline carbon 3 and/or 4");
			return;
		} else {
			IBond bond = molecule.getBond(gammaCarbon, deltaCarbon);
			bond.setOrder(IBond.Order.DOUBLE);
		}

	}
	
}
