package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.Atom;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

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
 * Execute thioamide formation as catalyzed by the thioviridamide enzyme TvaH.
 * 
 * @author skinnider
 *
 */
public class ThioamideReaction extends GenericReaction implements Reaction {

	public ThioamideReaction(ReactionPlan plan, Scaffold scaffold,
			Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.TvaH };
	}

	public void execute() throws NoResidueException,
			TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();

		for (int i = 0; i < plan.modules().size(); i++) {
			Module module = plan.get(i);
			Residue residue = scaffold.residue(module);
			IAtomContainer structure = residue.structure();

			// remove =O
			IAtom ketone = residue.ketone();
			IAtom ketoneOxygen = Atoms.getConnectedOxygen(ketone, structure);
			UtilityReactions.removeAtom(ketoneOxygen, molecule);

			// add =S
			IAtom sulfur = new Atom("S");
			molecule.addAtom(sulfur);
			UtilityReactions.addBond(ketone, sulfur, molecule,
					IBond.Order.DOUBLE);
		}
	}

}
