package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
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
 * Execute the reaction catalyzed by the proteusin radical SAM
 * N-methyltransferase PoyE, which catalyzes iterative asparagine side chain
 * N-methylation.
 * 
 * @author skinnider
 *
 */
public class PoyEReaction extends GenericReaction implements Reaction {
	
	public PoyEReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.PoyE };
	}
	
	public void execute() throws NoResidueException,
			TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		
		for (int i = 0; i < plan.size(); i++) {
			Module module = plan.get(i);
			Residue residue = scaffold.residue(module);
			if (residue == null)
				throw new TailoringSubstrateException(
						"Error: could not get residue for PoyE!");
			IAtomContainer structure = residue.structure();
			
			IAtom sideChainNitrogen = null;
			IAtom nitrogen = residue.nitrogen();
			for (IAtom atom : structure.atoms())
				if (atom.getSymbol().equals("N") && atom != nitrogen)
					sideChainNitrogen = atom;
			if (sideChainNitrogen == null)
				throw new TailoringSubstrateException(
						"Error: could not get residue for PoyE!");
			
			UtilityReactions.functionalize("CI", sideChainNitrogen, molecule);
		}
	}

}
