package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType.Hybridization;
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
import ca.mcmaster.magarveylab.prism.util.SmilesIO;

/**
 * Execute tryptophan isoprenylation as catalyzed by the ComX peptide
 * isoprenyltransferase ComQ.
 * 
 * @author skinnider
 *
 */
public class ComQReaction extends GenericReaction implements Reaction {

	public ComQReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.ComQ };
	}

	public void execute() throws NoResidueException,
			TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();

		Module module = plan.get(0);
		Residue residue = scaffold.residue(module);
		IAtomContainer structure = residue.structure();

		IAtom nitrogen = residue.nitrogen();

		// find isoprenylation site, carbon #2
		IAtom isoprenylationSite = null;
		IAtom carbon2 = null;
		for (IAtom atom : structure.atoms())
			if (atom.getSymbol().equals("N") && atom != nitrogen)
				for (IAtom atom2 : structure.getConnectedAtomsList(atom))
					if (structure.getConnectedBondsCount(atom2) == 2)
						for (IAtom atom3 : structure
								.getConnectedAtomsList(atom2))
							if (atom3 != atom) {
								isoprenylationSite = atom3;
								carbon2 = atom2;
							}

		if (carbon2 == null)
			throw new ScaffoldGenerationException("Error:"
					+ " couldn't get tryptophan carbon 2 for isoprenylation!");

		// reduce all bonds
		for (IBond bond : molecule.getConnectedBondsList(carbon2)) {
			IAtom connectedAtom = bond.getConnectedAtom(carbon2);
			UtilityReactions.setBondOrder(connectedAtom, carbon2, molecule,
					IBond.Order.SINGLE);
			connectedAtom.setHybridization(Hybridization.SP3);
		}
		carbon2.setHybridization(Hybridization.SP3);

		// create isoprenyl group
		String smiles = "IC/C=C(C)/CC/C=C(C)/C";
		IAtomContainer isoprenyl = null;
		try {
			isoprenyl = SmilesIO.molecule(smiles);
		} catch (Exception e) {
			throw new ScaffoldGenerationException(
					"Error encountered building isoprenyl group SMILES");
		}

		// get isoprenylation site
		IAtom iodine = Atoms.getIodine(isoprenyl);
		IAtom carbon = Atoms.getConnectedCarbon(iodine, isoprenyl);
		UtilityReactions.removeAtom(iodine, isoprenyl);

		// isoprenylate
		molecule.add(isoprenyl);
		structure.add(isoprenyl);
		UtilityReactions.addBond(carbon, isoprenylationSite, molecule);

		// cyclize indole again
		UtilityReactions.addBond(carbon2, nitrogen, molecule);

		// oxidize double bond
		UtilityReactions.setBondOrder(carbon2, isoprenylationSite, molecule,
				IBond.Order.SINGLE);
	}

}
