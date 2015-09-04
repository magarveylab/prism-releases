package ca.mcmaster.magarveylab.prism.cluster.reactions.typeII;

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.TypeIIPolyketideDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

public class AntibioticMonooxygenaseReaction extends GenericReaction implements Reaction {

	public AntibioticMonooxygenaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { TypeIIPolyketideDomains.ABM, 
				TypeIIPolyketideDomains.C6OX };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		Module first = plan.get(0);
		Module second = plan.get(1);
		
		Residue firstResidue = scaffold.residue(first);
		Residue secondResidue = scaffold.residue(second);
		IAtomContainer firstStructure = firstResidue.structure();

		// create ketone oxygen
		IAtomContainer oxygen = new AtomContainer();
		IAtom o = new Atom("O");
		oxygen.addAtom(o);
		firstStructure.add(oxygen);
		molecule.add(oxygen);
		
		// add to first residue
		IAtom alphaCarbon = firstResidue.alphaCarbon();
		UtilityReactions.addBond(alphaCarbon, o, molecule, IBond.Order.DOUBLE);
		
		// reduce alpha carbon/ketone bond 
		IAtom firstKetone = firstResidue.ketone();
		UtilityReactions.setBondOrder(firstKetone, alphaCarbon, molecule, IBond.Order.SINGLE);

		// oxidize ring closing bond
		IAtom secondAlphaCarbon = secondResidue.alphaCarbon();
		UtilityReactions.setBondOrder(firstKetone, secondAlphaCarbon, molecule, IBond.Order.DOUBLE);
		
		// oxidize second alcohol
		IAtom secondKetone = secondResidue.ketone();
		UtilityReactions.oxidizeAlcohol(secondKetone, molecule); 
		
		// reduce alpha carbon/ketone bond
		UtilityReactions.setBondOrder(secondKetone, secondAlphaCarbon, molecule, IBond.Order.SINGLE);
	}

}
