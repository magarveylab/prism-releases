package ca.mcmaster.magarveylab.prism.cluster.reactions.impl;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

/**
 * Execute the oxidation of thiazolines and oxazolines to thiazoles and
 * oxazoles, respectively.
 * 
 * @author skinnider
 *
 */
public class NitroreductaseReaction extends GenericReaction implements Reaction {

	public NitroreductaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { ThiotemplatedDomains.NITROREDUCTASE };
	}

	public void execute() throws NoResidueException,
			TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();

		Module module = plan.get(0);
		Residue residue = scaffold.residue(module);
		IAtom ketone = residue.ketone();
		IAtom alphaCarbon = residue.alphaCarbon();

		IAtom betaCarbon = null;
		for (IAtom atom : molecule.getConnectedAtomsList(alphaCarbon))
			if (atom.getSymbol().equals("C") && atom != ketone)
				betaCarbon = atom;

		if (betaCarbon == null)
			throw new TailoringSubstrateException("Error: could not execute thiazoline/oxazoline"
					+ " oxidation: could not find alpha/beta carbon bond");
		
		UtilityReactions.setBondOrder(alphaCarbon, betaCarbon, molecule, IBond.Order.DOUBLE);
	}

}
