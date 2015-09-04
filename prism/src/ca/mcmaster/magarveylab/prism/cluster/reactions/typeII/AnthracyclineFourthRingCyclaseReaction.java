package ca.mcmaster.magarveylab.prism.cluster.reactions.typeII;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.TypeIIPolyketideDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Atoms;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

public class AnthracyclineFourthRingCyclaseReaction extends GenericReaction implements Reaction {
	
	public AnthracyclineFourthRingCyclaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { TypeIIPolyketideDomains.CYCLASE_CLADE_6a, 
				TypeIIPolyketideDomains.CYCLASE_CLADE_6b_SUBTYPE_1 };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();

		Residue r1 = scaffold.residue(plan.get(0));
		Residue r2 = scaffold.residue(plan.get(1));
		Residue r3 = scaffold.residue(plan.get(2));
		
		IAtom c1 = r1.ketone();
		IAtom c2 = r1.alphaCarbon();
		IAtom c17 = r2.ketone();
		IAtom c19 = r3.ketone();
		
		// add bond between C19-C2
		UtilityReactions.addBond(c2, c19, molecule);
		
		// reduce C19 ketone
		UtilityReactions.reduceKetone(c19, molecule);
		
		// reduce C17 ketone
		UtilityReactions.reduceKetone(c17, molecule);
		
		// if C1 is COOH, decarboxylate; if methylated, do not
		boolean methyl = false;
		IAtom alcohol = Atoms.getCarboxylAlcohol(c1, molecule);
		for (IAtom atom : molecule.getConnectedAtomsList(alcohol))
			if (atom.getSymbol().equals("C") && atom != c1) 
				methyl = true;
		if (!methyl)
			UtilityReactions.decarboxylate(c1, c2, molecule);
	}

}
