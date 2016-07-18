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
import ca.mcmaster.magarveylab.prism.cluster.reactions.impl.HeterocyclizationReaction;
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
 * Execute trifolitoxin chromophore formation.
 * 
 * @author skinnider
 *
 */
public class TfxBReaction extends GenericReaction implements Reaction {

	public TfxBReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.TfxB };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module m1 = plan.get(0);
		Module m2 = plan.get(1);
		Module m3 = plan.get(2);
		Module m4 = plan.get(3);
		Residue r1 = scaffold.residue(m1);
		Residue r2 = scaffold.residue(m2);
		IAtomContainer s2 = r2.structure();
		
		// oxidize 1st residue alpha carbon-nitrogen bond
		IAtom alphaCarbon1 = r1.alphaCarbon();
		IAtom nitrogen1 = r1.nitrogen();
		UtilityReactions.setBondOrder(nitrogen1, alphaCarbon1, molecule, IBond.Order.DOUBLE);
		
		// oxidize 2nd residue (Q) alpha carbon-nitrogen bond
		IAtom alphaCarbon2 = r2.alphaCarbon();
		IAtom nitrogen2 = r2.nitrogen();
		UtilityReactions.setBondOrder(nitrogen2, alphaCarbon2, molecule, IBond.Order.DOUBLE);
		
		// remove 2nd residue (Q) ketone oxygen
		IAtom ketone2 = r2.ketone();
		IAtom ketoneOxygen2 = Atoms.getConnectedOxygen(ketone2, s2);
		UtilityReactions.removeAtom(ketoneOxygen2, molecule);
		
		// get 2nd residue (Q) terminal nitrogen
		IAtom terminalNitrogen2 = null;
		for (IAtom atom : s2.atoms())
			if (atom.getSymbol().equals("N") && atom != nitrogen2)
				terminalNitrogen2 = atom;
		
		// add double bond between Q terminal nitrogen and ketone
		UtilityReactions.addBond(terminalNitrogen2, ketone2, molecule, IBond.Order.DOUBLE);
		
		// get Q gamma carbon
		IAtom gammaCarbon = null;
		IAtom betaCarbon = RibosomalUtil.getBetaCarbon(r2, s2);
		for (IAtom atom : s2.getConnectedAtomsList(betaCarbon))
			if (atom.getSymbol().equals("C") && atom != betaCarbon)
				gammaCarbon = atom;
		
		// add ketone to Q gamma carbon
		IAtom newKetone = new Atom("O");
		molecule.addAtom(newKetone);
		s2.addAtom(newKetone);
		UtilityReactions.addBond(gammaCarbon, newKetone, molecule, IBond.Order.DOUBLE);
		
		// heterocyclization of cysteine to thiazoline
		SubstrateSet newSubstrate = new SubstrateSet(m3, m4);
		ReactionPlan newPlan = new ReactionPlan(plan.domain(), newSubstrate, plan.reaction());
		HeterocyclizationReaction newReaction = new HeterocyclizationReaction(newPlan, scaffold, cluster);
		newReaction.execute();
	}
	
}
