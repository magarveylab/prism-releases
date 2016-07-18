package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import java.util.List;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;
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
 * Execute the reaction catalyzed by homologs of the TP-1161 TpaJ enzmye, which
 * is speculated to catalyze conversion of a C-terminal threonine residue to an
 * amino acetone group by oxidation.
 * 
 * @author skinnider
 *
 */
public class TpaJReaction extends GenericReaction implements Reaction {

	public TpaJReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.TpaJ };
	}

	public void execute() throws NoResidueException,
			TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();

		Module module = plan.get(0);
		Residue residue = scaffold.residue(module);
		if (residue == null)
			throw new TailoringSubstrateException(
					"Error: could not get residue for TpaJ!");
		IAtomContainer structure = residue.structure();

		if (module.type() == ModuleTypes.RIBOSOMAL
				&& module.scaffold() != null
				&& module.scaffold().topSubstrate() != null
				&& module.scaffold().topSubstrate().type() == ProteinogenicAminoAcids.THREONINE) {
			// get carboxyl
			IAtom carboxyl = null;
			List<IAtom> carboxyls = Atoms.getAllCarboxyls(molecule);
			for (IAtom atom : structure.atoms())
				if (carboxyls.contains(atom))
					carboxyl = atom;
			if (carboxyl == null)
				throw new TailoringSubstrateException(
						"Error: could not get carboxyl group for TpaJ!");

			// decarboxylate
			IAtom alphaCarbon = residue.alphaCarbon();
			UtilityReactions.decarboxylate(carboxyl, alphaCarbon, molecule);

			// oxidize
			IAtom betaCarbon = RibosomalUtil.getBetaCarbon(residue, structure);
			IAtom oxygen = RibosomalUtil.getSerineOrThreonineHydroxyl(residue,
					structure, molecule);
			if (oxygen == null)
				throw new TailoringSubstrateException(
						"Error: could not get hydroxyl group for TpaJ!");
			IBond bond = molecule.getBond(betaCarbon, oxygen);
			if (bond == null)
				throw new TailoringSubstrateException(
						"Error: could not get beta carbon-hydroxyl bond for TpaJ!");
			UtilityReactions.setBondOrder(betaCarbon, oxygen, molecule,
					IBond.Order.DOUBLE);
		} else {
			throw new TailoringSubstrateException(
					"Error: could not get threonine residue for TpaJ!");
		}
	}

}
