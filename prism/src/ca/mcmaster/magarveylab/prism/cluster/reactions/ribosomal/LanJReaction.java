package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

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
 * Execute dehydroalanine oxidation to D-Ala (no chirality included in this
 * reaction as PRISM structures are achiral).
 * 
 * @author skinnider
 *
 */
public class LanJReaction extends GenericReaction implements Reaction {

	public LanJReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.LanJ };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		for (int i = 0; i < plan.modules().size(); i++) {
			Module module = plan.get(i);
			Residue residue = scaffold.residue(module);
			IAtomContainer structure = residue.structure();
			IAtom alphaCarbon = residue.alphaCarbon();
			IAtom betaCarbon = RibosomalUtil.getBetaCarbon(residue, structure);
			
			// must be dehydroalanine
			if (!RibosomalUtil.isDehydrated(residue, structure, molecule))
				throw new TailoringSubstrateException("Error: could not oxidize serine "
						+ "residue which has not been dehydrated!");
			
			// there cannot be anything attached to the beta carbon! 
			for (IAtom atom : molecule.getConnectedAtomsList(betaCarbon))
				if (atom != alphaCarbon)
					throw new TailoringSubstrateException("Error: could not oxidize "
							+ "Dha residue with more than 1 bond to beta carbon!");
			
			// oxidize double bond
			UtilityReactions.setBondOrder(alphaCarbon, betaCarbon, molecule, IBond.Order.SINGLE);
		}
	}
	
}
