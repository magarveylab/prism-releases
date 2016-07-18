package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.Atom;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.RibosomalUtil;
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
 * Execute gamma-hydroxylation of an amino acid substrate.
 * 
 * @author skinnider
 *
 */
public class GammaHydroxylaseReaction extends GenericReaction implements Reaction {
	
	public GammaHydroxylaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.NosB };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module module = plan.get(0);
		Residue residue = scaffold.residue(module);
		if (residue == null) 
			throw new TailoringSubstrateException("Error: could not get residue for beta-hydroxylation!");
		IAtomContainer structure = residue.structure();
		
		IAtom alphaCarbon = residue.alphaCarbon();
		IAtom betaCarbon = RibosomalUtil.getBetaCarbon(residue, structure);
		IAtom gammaCarbon = null;
		for (IAtom atom : structure.getConnectedAtomsList(betaCarbon))
			if (atom.getSymbol().equals("C") && atom != alphaCarbon)
				gammaCarbon = atom;
		
		IAtom oxygen = new Atom("O");
		structure.addAtom(oxygen);
		molecule.addAtom(oxygen);
		UtilityReactions.addBond(gammaCarbon, oxygen, molecule);
	}

}
