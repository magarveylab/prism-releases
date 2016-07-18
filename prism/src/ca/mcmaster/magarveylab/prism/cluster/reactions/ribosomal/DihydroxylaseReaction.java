package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

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
 * Execute isoleucine dihydroxylation, as catalyzed by the TsrR/SioR enzymes.
 * 
 * @author skinnider
 *
 */
public class DihydroxylaseReaction extends GenericReaction implements Reaction {
	
	public DihydroxylaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.TsrR };
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
		
		String smiles = "OI";
		UtilityReactions.functionalize(smiles, betaCarbon, molecule);
		UtilityReactions.functionalize(smiles, gammaCarbon, molecule);
	}

}
