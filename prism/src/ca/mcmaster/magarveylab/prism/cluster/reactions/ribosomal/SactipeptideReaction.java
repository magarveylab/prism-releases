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
 * Execute sactipeptide thioether bond formation between a cysteine sulfhydryl
 * group and the alpha carbon of another residue.
 * 
 * @author skinnider
 *
 */
public class SactipeptideReaction extends GenericReaction implements Reaction {

	public SactipeptideReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.AlbA };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		for (int i = 0; i < plan.modules().size() - 1; i += 2) {
			Module m1 = plan.get(i);
			Module m2 = plan.get(i+1);
			Residue r1 = scaffold.residue(m1);
			Residue r2 = scaffold.residue(m2);
			IAtomContainer s1 = r1.structure();
			
			// cysteine cannot be cyclized
			IAtom sulfur = Atoms.getSulfur(s1);
			if (sulfur == null)
				throw new ScaffoldGenerationException("Could not execute sactipeptide bond formation: "
						+ "could not locate sulfur atom!");
			if (molecule.getConnectedBondsCount(sulfur) != 1)
				throw new ScaffoldGenerationException("Could not execute sactipeptide bond formation: "
						+ "cysteine has already been cyclized!");

			// add bond between alpha carbon and sulfur
			IAtom alphaCarbon = r2.alphaCarbon();
			UtilityReactions.addBond(sulfur, alphaCarbon, molecule);
		}
	}
	
}
