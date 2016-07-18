package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Atoms;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Bonds;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

/**
 * Execute oxazoline/thiazoline formation in ribosomal peptides.
 * 
 * @author skinnider
 *
 */
public class AzolineReaction extends GenericReaction implements Reaction {

	public AzolineReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { 
				RibosomalDomains.McbC, 
				RibosomalDomains.McbD,
				RibosomalDomains.PatD,
				RibosomalDomains.LazE,
				RibosomalDomains.YmB1,
				RibosomalDomains.YmBC_a,
				RibosomalDomains.BotCD,
				RibosomalDomains.GodD  //arbitrarily chosen over GodE
		};
	}

	public void execute() throws NoResidueException,
			TailoringSubstrateException, ScaffoldGenerationException,
			CDKException {
		SubstrateSet substrates = plan.modules();
				
		for (int i = 0; i < substrates.size() - 1; i += 2) {
			Module m1 = plan.get(i);
			Module m2 = plan.get(i + 1);
			
			IAtomContainer molecule = scaffold.molecule();
			
			Module last = m2;
			Module module = m1;
			Residue residue = scaffold.residue(module);
			Residue lastResidue = scaffold.residue(last);
			if (residue == null || lastResidue == null)
				throw new ScaffoldGenerationException(
						"Error: could not heterocyclize on null residue!");

			IAtom ketone = residue.ketone();
			if (ketone == null)
				throw new ScaffoldGenerationException("Could not find ketone for ser/thr/cys residue");
			IAtom lastKetone = lastResidue.ketone();
			if (lastKetone == null)
				throw new ScaffoldGenerationException("Could not find ketone for ser/thr/cys-adjacent residue");
			IAtom alphaCarbon = residue.alphaCarbon();
			if (alphaCarbon == null)
				throw new ScaffoldGenerationException("Could not find alpha carbon for ser/thr/cys residue");
			IAtom betaCarbon = null;
			IAtom terminalAtom = null;
						
			for (IBond bond : molecule.getConnectedBondsList(alphaCarbon)) {
				IAtom connectedAtom = bond.getConnectedAtom(alphaCarbon);
				if (connectedAtom.getSymbol().equals("C") && connectedAtom != ketone)
					betaCarbon = connectedAtom;
			}
			
			if (betaCarbon == null)
				throw new ScaffoldGenerationException("Could not generate thiazoline/oxazoline:"
						+ " could not find beta carbon!");
			
			for (IBond bond : molecule.getConnectedBondsList(betaCarbon)) {
				IAtom connectedAtom = bond.getConnectedAtom(betaCarbon);
				if (connectedAtom.getSymbol().equals("O") 
						|| connectedAtom.getSymbol().equals("S"))
					terminalAtom = connectedAtom;
			}
			
			if (terminalAtom == null) 
				throw new ScaffoldGenerationException("Could not generate thiazoline/oxazoline: "
						+ "could not find terminal oxygen or sulfur atom!");

			// terminalAtom can't have 2 bonds (e.g. methylation, esterification)
			if (molecule.getConnectedBondsCount(terminalAtom) == 2) 
				throw new TailoringSubstrateException("Error: could not generate "
						+ "thiazoline/oxazoline: O/S has 2 bonds!");

			// remove ketone oxygen
			try {
				IAtom ketoneOxygen = Atoms.getConnectedOxygen(lastKetone, molecule); 
				if (ketoneOxygen == null)
					throw new ScaffoldGenerationException("Couldn't identify ketone oxygen!");
				
				IBond oxygenBond = Bonds.getConnectedOxygenBond(molecule, lastKetone);
				if (oxygenBond == null)
					throw new ScaffoldGenerationException("Couldn't get ketone oxygen bond!");
				if (oxygenBond.getOrder() != IBond.Order.DOUBLE) 
					throw new ScaffoldGenerationException("Couldn't remove single-bonded ketone oxygen!");
				
				molecule.removeBond(oxygenBond);
				molecule.removeAtom(ketoneOxygen);
			} catch (NullPointerException e) {
				throw new ScaffoldGenerationException("Could not remove ketone oxygen!");
			}

			// create double bond between ketone carbon and nitrogen
			try {
				IBond nitrogenBond = Bonds.getConnectedNitrogenBond(molecule, lastKetone);
				nitrogenBond.setOrder(IBond.Order.DOUBLE);
			} catch (NullPointerException e) {
				throw new ScaffoldGenerationException("Could not add double bond between ketone carbon and nitrogen!");
			}
			
			// create bond between terminal O/S and ketone carbon
			UtilityReactions.addBond(lastKetone, terminalAtom, molecule);
		}
	}

}
