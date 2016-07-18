package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

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
 * Execute dehydrothreonine O-methylation, as catalyzed by the nocathiacin NocQ
 * enzyme.
 * 
 * @author skinnider
 *
 */
public class NocQReaction extends GenericReaction implements Reaction {

	public NocQReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.NocQ };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module module = plan.get(0);
		Residue residue = scaffold.residue(module);
		IAtomContainer structure = residue.structure();
		IAtom ketone = residue.ketone();
		
		// must be dehydrated
		IAtom alphaCarbon = residue.alphaCarbon();
		IAtom betaCarbon = null;
		for (IAtom atom : structure.getConnectedAtomsList(alphaCarbon))
			if (atom.getSymbol().equals("C") && atom != ketone)
				betaCarbon = atom;
		IBond bond = structure.getBond(alphaCarbon, betaCarbon);
		if (bond == null)
			throw new TailoringSubstrateException("Could not O-methylate dehydrothreonine: "
					+ "could not get alpha/beta carbon bond!");
		if (bond.getOrder() != IBond.Order.DOUBLE)
			throw new TailoringSubstrateException("Could not O-methylate dehydrothreonine: "
					+ "alpha/beta carbon bond order is not double!");
		
		UtilityReactions.functionalize("COI", betaCarbon, molecule);
	}
	
}
