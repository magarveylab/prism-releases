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
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

/**
 * Execute oxazole/thiazole formation.
 * 
 * @author skinnider
 *
 */
public class AzoleReaction extends GenericReaction implements Reaction {

	public AzoleReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { 
				RibosomalDomains.McbB,
				RibosomalDomains.LazF,
				RibosomalDomains.YmC1,
				RibosomalDomains.YmBC_b,
				RibosomalDomains.GodE // arbitarily chosen over GodD
			};
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		SubstrateSet substrates = plan.modules();
		
		for (int i = 0; i < substrates.size() - 1; i += 2) {
			IAtomContainer molecule = scaffold.molecule();

			Module module = plan.get(i);
			Residue residue = scaffold.residue(module);
			IAtomContainer structure = residue.structure();
			
			IAtom alphaCarbon = residue.alphaCarbon();
			IAtom betaCarbon = RibosomalUtil.getBetaCarbon(residue, structure);
			if (betaCarbon == null)
				throw new TailoringSubstrateException("Error: could not execute thiazoline/oxazoline"
						+ " oxidation: could not find alpha/beta carbon bond");
			
			IAtom terminalAtom = null;
			for (IBond bond : molecule.getConnectedBondsList(betaCarbon)) {
				IAtom connectedAtom = bond.getConnectedAtom(betaCarbon);
				if (connectedAtom.getSymbol().equals("O") 
						|| connectedAtom.getSymbol().equals("S"))
					terminalAtom = connectedAtom;
			}
			
			if (terminalAtom == null)
				throw new ScaffoldGenerationException("Could not execute thiazoline/oxazoline: "
						+ "oxidation: could not find terminal oxygen or sulfur atom!");

			if (molecule.getConnectedBondsCount(terminalAtom) == 1)
				continue;
			
			UtilityReactions.setBondOrder(alphaCarbon, betaCarbon, molecule,
					IBond.Order.DOUBLE);
		}
	}

}
