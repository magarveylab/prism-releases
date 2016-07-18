package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.RibosomalUtil;
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
 * Execute cysteine S-glucosylation as catalyzed by the glycocin
 * glycosyltransferase SunS, or O-glycosylation when the hairpin cysteine is
 * glycosylated.
 * 
 * @author skinnider
 *
 */
public class SunSReaction extends GenericReaction implements Reaction {

	public SunSReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.SunS };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		for (int i = 0; i < plan.size(); i++) {
			Module module = plan.get(i);
			Residue residue = scaffold.residue(module);
			IAtomContainer structure = residue.structure();

			// get sulfur
			IAtom site = Atoms.getSulfur(structure);
			
			// if sulfur is null, it might be a serine 
			if (site == null)
				site = RibosomalUtil.getSerineOrThreonineHydroxyl(residue, structure, molecule);
			
			// if it's still null, throw an error 
			if (site == null)
				throw new TailoringSubstrateException("Could not execute S-glycosylation: "
						+ "could not locate cysteine sulfur!");
			if (molecule.getConnectedBondsCount(site) != 1)
				throw new TailoringSubstrateException("Could not execute S-glycosylation: "
						+ "cysteine sulfur has >1 bond!");

			// glycosylate
			String glucose = "OCC1OC(C(O)C(O)C1O)I";
			UtilityReactions.functionalize(glucose, site, molecule);
		}
	}
	
}
