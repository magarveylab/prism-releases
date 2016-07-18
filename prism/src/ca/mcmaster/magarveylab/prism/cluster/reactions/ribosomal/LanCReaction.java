package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

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
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

/**
 * Execute the reaction catalyzed by the lantibiotic cyclase LanC, which
 * catalyzes lantithione bond formation between a cysteine residue and a
 * dehydroalanine/dehydroaminobutyric acid residue.
 * 
 * @author skinnider
 *
 */
public class LanCReaction extends GenericReaction implements Reaction {
	
	public LanCReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.LanC };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		for (int i = 0; i < plan.modules().size() - 1; i += 2) {
			Module m1 = plan.get(i);
			Residue dhaDhb = scaffold.residue(m1);
			if (dhaDhb == null) 
				throw new TailoringSubstrateException("Error: could not get Dha/Dhb residue for LanC!");
			IAtomContainer dhaDhbStructure = dhaDhb.structure();
			
			Module m2 = plan.get(i+1);
			Residue cysteine = scaffold.residue(m2);
			if (cysteine == null) 
				throw new TailoringSubstrateException("Error: could not get Cys residue for LanC!");
			IAtomContainer cysteineStructure = cysteine.structure();
			
			// get cysteine sulfur
			IAtom sulfur = Atoms.getSulfur(cysteineStructure);
			if (sulfur == null)
				throw new TailoringSubstrateException("Error: could not get cysteine sulfur!");
			
			// cysteine cannot already be participating in another reaction 
			if (molecule.getConnectedBondsCount(sulfur) > 1) 
				continue;

			// get Dha/Dhb beta carbon
			IAtom alphaCarbon = dhaDhb.alphaCarbon();
			IAtom betaCarbon = RibosomalUtil.getBetaCarbon(dhaDhb, dhaDhbStructure);
			
			// undo tautomerization
			for (IAtom atom : molecule.getConnectedAtomsList(alphaCarbon))
				if (atom.getSymbol().equals("O")) {
					IBond bond = molecule.getBond(atom, alphaCarbon);
					if (bond.getOrder() == IBond.Order.DOUBLE) {
						UtilityReactions.removeAtom(atom, molecule);
						UtilityReactions.setBondOrder(betaCarbon, alphaCarbon, molecule, IBond.Order.DOUBLE);
						UtilityReactions.functionalize("NI", alphaCarbon, molecule);
					}
				}
			
			// can't already be lanthionine
			for (IAtom atom : molecule.getConnectedAtomsList(betaCarbon))
				if (atom.getSymbol().equals("S"))
					throw new TailoringSubstrateException(
							"Error: could not form lanthionine: residue is already the site of lanthionine!");

			if (!RibosomalUtil.isDehydrated(dhaDhb, dhaDhbStructure, molecule))
				if (plan.domain().type() == RibosomalDomains.LanKC) {
					// dehydration has not yet been performed for a lanthionine (NOT labionin)
					SubstrateSet substrate = new SubstrateSet(m1);
					ReactionPlan lanbPlan = new ReactionPlan(plan.domain(), substrate,
							plan.reaction());
					LanBReaction lanb = new LanBReaction(lanbPlan, scaffold, cluster);
					lanb.execute();
				} else {
					// Dhb/Dha cannot still have Ser/Thr oxygen 
					throw new TailoringSubstrateException("Error: serine/threonine residue has not been dehydrated!");
				}
			
			// add bond
			UtilityReactions.addBond(betaCarbon, sulfur, molecule);
			
			// set alpha/beta bond order
			UtilityReactions.setBondOrder(alphaCarbon, betaCarbon, molecule, IBond.Order.SINGLE);
		}
	}

}
