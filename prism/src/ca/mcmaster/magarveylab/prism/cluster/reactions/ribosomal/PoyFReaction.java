package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.Atom;
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
 * Execute the reaction catalyzed by the proteusin dehydratase PoyF, which
 * catalyzes the dehydration of N-terminal threonine (or, maybe, serine)
 * residues, with subsequent leader cleavage resulting in beta-ketone formation.
 * 
 * @author skinnider
 *
 */
public class PoyFReaction extends GenericReaction implements Reaction {
	
	public PoyFReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.PoyF };
	}
	
	public void execute() throws NoResidueException,
			TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		
		Module module = plan.get(0);
		Residue residue = scaffold.residue(module);
		if (residue == null)
			throw new TailoringSubstrateException(
					"Error: could not get residue for PoyF!");
		IAtomContainer structure = residue.structure();
		
		// remove nitrogen
		IAtom nitrogen = residue.nitrogen();
		UtilityReactions.removeAtom(nitrogen, molecule);
		
		// add =O to alpha carbon 
		IAtom alphaCarbon = residue.alphaCarbon();
		IAtom oxygen = new Atom("O");
		molecule.addAtom(oxygen);
		UtilityReactions.addBond(alphaCarbon, oxygen, molecule, IBond.Order.DOUBLE);

		// remove Ser/Thr -OH 
		IAtom hydroxyl = RibosomalUtil.getSerineOrThreonineHydroxyl(residue, structure, molecule);
		UtilityReactions.removeAtom(hydroxyl, molecule);
	}

}
