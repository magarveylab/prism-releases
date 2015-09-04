package ca.mcmaster.magarveylab.prism.cluster.reactions.impl;

import java.util.List;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Atoms;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

public class P450CReaction extends GenericReaction implements Reaction {

	public P450CReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { TailoringDomains.P450C };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();

		Module firstModule = plan.get(0);
		Module secondModule = plan.get(1);
		Residue firstResidue = scaffold.residue(firstModule);
		Residue secondResidue = scaffold.residue(secondModule);
		IAtomContainer first = firstResidue.structure();
		IAtomContainer second = secondResidue.structure();

		// get hpg meta carbon
		IAtom paraCarbon = Atoms.getAromaticHydroxylCarbon(first);
		IAtom c1 = Atoms.getConnectedCarbon(paraCarbon, first);

		// get dhpg ortho carbon
		List<IAtom> ketones = Atoms.getAllKetones(second);
		IAtom ketone = ketones.get(0);
		IAtom alphaCarbon = Atoms.getConnectedCarbon(ketone, second);
		IAtom rCarbon = null;
		IAtom orthoCarbon = null;
		for (IAtom atom : second.getConnectedAtomsList(alphaCarbon))
			if (atom.getSymbol().equals("C") && atom != ketone)
				rCarbon = atom;
		for (IAtom atom : second.getConnectedAtomsList(rCarbon))
			if (atom.getSymbol().equals("C") && atom != alphaCarbon)
				orthoCarbon = atom;		

		// add bond
		if (molecule.getBond(c1, orthoCarbon) == null)
			UtilityReactions.addBond(c1, orthoCarbon, molecule);
	}

}
