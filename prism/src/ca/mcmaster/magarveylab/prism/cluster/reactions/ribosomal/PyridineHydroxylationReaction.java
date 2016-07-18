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
 * Execute pyridine ring hydroxylation, as catalyzed by the NosB/NosC enzymes.
 * 
 * @author skinnider
 *
 */
public class PyridineHydroxylationReaction extends GenericReaction implements Reaction {
	
	public PyridineHydroxylationReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.NosC };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module m1 = plan.get(0);
		Residue r1 = scaffold.residue(m1);
		if (r1 == null) 
			throw new TailoringSubstrateException("Error: could not get first Dha residue "
					+ "for pyridine hydroxylation!");
		IAtomContainer s1 = r1.structure();
		
		Module m2 = plan.get(1);
		Residue r2 = scaffold.residue(m2);
		if (r2 == null) 
			throw new TailoringSubstrateException("Error: could not get second Dha residue "
					+ "for pyridine hydroxylation!");
		IAtomContainer s2 = r2.structure();
		
		// both residues must be dehydrated
		if (!RibosomalUtil.isDehydrated(r1, s1, molecule)) 
			throw new TailoringSubstrateException("Error: could not execute pyridine "
					+ "hydroxylation with first non-dehydrated serine!");
		if (!RibosomalUtil.isDehydrated(r2, s2, molecule)) 
			throw new TailoringSubstrateException("Error: could not execute pyridine "
					+ "hydroxylation with second non-dehydrated serine!");
		
		// must be in a ring
		IAtom betaCarbon1 = RibosomalUtil.getBetaCarbon(r1, s1);
		IAtom betaCarbon2 = RibosomalUtil.getBetaCarbon(r2, s2);
		if (molecule.getBond(betaCarbon1, betaCarbon2) == null)
			throw new TailoringSubstrateException("Error: could not execute pyridine "
					+ "hydroxylation--residues are not in a ring!");
		
		// hydroxylate
		String smiles = "OI";
		UtilityReactions.functionalize(smiles, betaCarbon2, molecule);
	}
	
}
