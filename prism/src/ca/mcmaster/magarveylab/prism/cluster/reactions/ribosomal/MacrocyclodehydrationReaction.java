package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

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
 * Execute the macrocyclodehydration of a bottromycin-family ribosomal peptide
 * as putatively catalyzed by the BotC enzyme.
 * 
 * @author skinnider
 *
 */
public class MacrocyclodehydrationReaction extends GenericReaction implements Reaction {

	public MacrocyclodehydrationReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.BotC };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module m1 = plan.get(0);
		Residue r1 = scaffold.residue(m1);
		Module m2 = plan.get(1);
		Residue r2 = scaffold.residue(m2);
		Module m3 = plan.get(2);
		Residue r3 = scaffold.residue(m3);
		if (r1 == null)
			throw new TailoringSubstrateException("Error: could not get first "
					+ "residue for macrocyclodehydration!");
		if (r2 == null)
			throw new TailoringSubstrateException("Error: could not get fourth "
					+ "residue for macrocyclodehydration!");
		if (r3 == null)
			throw new TailoringSubstrateException("Error: could not get fifth "
					+ "residue for macrocyclodehydration!");
		IAtomContainer s2 = r2.structure();

		// check cyclization atoms
		IAtom nitrogen1 = r1.nitrogen();
		if (molecule.getConnectedBondsCount(nitrogen1) > 2)
			throw new ScaffoldGenerationException("Could not macrocyclize: "
					+ "nitrogen atom has >2 bonds");

		// remove ketone
		IAtom ketone = r2.ketone();
		IAtom ketoneOxygen = Atoms.getConnectedOxygen(ketone, s2);
		UtilityReactions.removeAtom(ketoneOxygen, molecule);
		
		// add bond
		UtilityReactions.addBond(ketone, nitrogen1, molecule);

		// set double bond
		IAtom nitrogen2 = r3.nitrogen();
		UtilityReactions.setBondOrder(ketone, nitrogen2, molecule, IBond.Order.DOUBLE);
	}
	
}
