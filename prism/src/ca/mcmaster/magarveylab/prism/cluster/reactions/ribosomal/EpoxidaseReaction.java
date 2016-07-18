package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.RibosomalUtil;
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

/**
 * Execute isoleucine epoxidation, as catalyzed by the TpdJ1 or 2 enzyme.
 * 
 * @author skinnider
 *
 */
public class EpoxidaseReaction extends GenericReaction implements Reaction {
	
	public EpoxidaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module module = plan.get(0);
		Residue residue = scaffold.residue(module);
		if (residue == null) 
			throw new TailoringSubstrateException("Error: could not get residue for epoxidation!");
		IAtomContainer structure = residue.structure();
		
		IAtom alphaCarbon = residue.alphaCarbon();
		IAtom betaCarbon = RibosomalUtil.getBetaCarbon(residue, structure);
		IAtom gammaCarbon = null;
		for (IAtom atom : structure.getConnectedAtomsList(betaCarbon))
			if (atom.getSymbol().equals("C")
					&& molecule.getConnectedBondsCount(atom) == 2
					&& atom != alphaCarbon)
				gammaCarbon = atom;
		if (gammaCarbon == null)
			throw new TailoringSubstrateException("Error: could not get gamma carbon for epoxidation!");
		IAtom deltaCarbon = null;
		for (IAtom atom : structure.getConnectedAtomsList(gammaCarbon))
			if (atom.getSymbol().equals("C")
					&& atom != betaCarbon)
				deltaCarbon = atom;
		if (deltaCarbon == null)
			throw new TailoringSubstrateException("Error: could not get delta carbon for epoxidation!");
		
		// make sure there isn't already a gamma-hydroxylation
		boolean hasGammaHydroxylation = false;
		for (IAtom atom : molecule.getConnectedAtomsList(gammaCarbon))
			if (atom.getSymbol().equals("O")
					&& molecule.getBond(atom, gammaCarbon).getOrder() == IBond.Order.DOUBLE)
				hasGammaHydroxylation = true;
		if (hasGammaHydroxylation)
			throw new TailoringSubstrateException(
					"Error: could not epoxidize a residue "
							+ "which is already gamma-hydroxylated!");
		
		// hydroxylate beta carbon
		String smiles = "OI";
		UtilityReactions.functionalize(smiles, gammaCarbon, molecule);
		
		// add epoxide bond
		IAtom oxygen = Atoms.getConnectedOxygen(gammaCarbon, molecule);
		UtilityReactions.addBond(oxygen, deltaCarbon, molecule);
		oxygen.setImplicitHydrogenCount(0);
	}

}
