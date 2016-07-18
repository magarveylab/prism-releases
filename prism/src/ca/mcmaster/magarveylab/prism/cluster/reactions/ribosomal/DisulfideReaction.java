package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
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
 * Execute disulfide bond formation. 
 * 
 * @author skinnider
 *
 */
public class DisulfideReaction extends GenericReaction implements Reaction {

	public DisulfideReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.Lasso_precursors,
				RibosomalDomains.BdbB, RibosomalDomains.SkfH,
				RibosomalDomains.SunA };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		for (int i = 0; i < plan.modules().size() - 1; i += 2) {
			Module m1 = plan.get(i);
			Module m2 = plan.get(i+1);
			Residue r1 = scaffold.residue(m1);
			Residue r2 = scaffold.residue(m2);
			IAtomContainer s1 = r1.structure();
			IAtomContainer s2 = r2.structure();
			
			// cysteines cannot be cyclized
			IAtom sulfur1 = Atoms.getSulfur(s1);
			if (molecule.getConnectedBondsCount(sulfur1) != 1)
				throw new ScaffoldGenerationException("Could not execute disulfide bond formation: "
						+ "cysteine #1 has already been cyclized!");
			IAtom sulfur2 = Atoms.getSulfur(s2);
			if (molecule.getConnectedBondsCount(sulfur2) != 1)
				throw new ScaffoldGenerationException("Could not execute disulfide bond formation: "
						+ "cysteine #2 has already been cyclized!");

			// add bond between alpha carbon and sulfur
			UtilityReactions.addBond(sulfur1, sulfur2, molecule);
		}
	}
	
}
