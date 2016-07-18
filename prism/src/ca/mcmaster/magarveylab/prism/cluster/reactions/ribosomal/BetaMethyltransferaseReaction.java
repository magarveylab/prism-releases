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
 * Execute beta-C-methylation of an amino acid substrate, as catalyzed by 
 * bottromycin C-methyltransferases and the proteusin PoyB/PoyC methyltransferases.
 * 
 * @author skinnider
 *
 */
public class BetaMethyltransferaseReaction extends GenericReaction implements Reaction {
	
	public BetaMethyltransferaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.BotRMT1,
				RibosomalDomains.BotRMT2, RibosomalDomains.BotRMT3,
				RibosomalDomains.PoyBC };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
				
		for (Module module : plan.modules().getAllModules()) {
			Residue residue = scaffold.residue(module);
			if (residue == null) 
				throw new TailoringSubstrateException(
						"Error: could not get residue for beta-methylation!");
			IAtomContainer structure = residue.structure();			
			IAtom betaCarbon = RibosomalUtil.getBetaCarbon(residue, structure);
			
			// make sure a methyl group can be added 
			if (molecule.getConnectedBondsCount(betaCarbon) >= 4)
				throw new TailoringSubstrateException(
						"Error: could not methylate beta carbon with 4 or more bonds!");
				
			IAtom c = new Atom("C");
			structure.addAtom(c);
			molecule.addAtom(c);
			UtilityReactions.addBond(betaCarbon, c, molecule);
		}
	}

}
