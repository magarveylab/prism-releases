package ca.mcmaster.magarveylab.prism.cluster.reactions.typeII;

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
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

public class C2MethyltransferaseReaction extends GenericReaction implements Reaction {

	public C2MethyltransferaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { TypeIIPolyketideDomains.C2MT };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		Residue residue = scaffold.residue(plan.get(0));
		
		IAtom oxygen = null;
		IAtom nitrogen = Atoms.getPrimaryAmine(residue.structure(), molecule);
		if (nitrogen == null)
			oxygen = Atoms.getHydroxyl(residue.structure(), molecule);
		if (oxygen == null && nitrogen == null)
			throw new ScaffoldGenerationException("Error: could not find atom for D-ring N-/O-methyltransferase!");
		
		// create methyl group
		IAtomContainer methane = new AtomContainer();
		IAtom carbon = new Atom("C");
		methane.addAtom(carbon);

		// add to scaffold & add bond
		if (nitrogen != null) {
			molecule.add(methane);
			residue.structure().add(methane);
			UtilityReactions.addBond(carbon, nitrogen, molecule);
			
			// create 2nd methyl group
			IAtomContainer methane2 = new AtomContainer();
			IAtom carbon2 = new Atom("C");
			methane2.addAtom(carbon2);
			molecule.add(methane2);
			residue.structure().add(methane2);
			UtilityReactions.addBond(carbon2, nitrogen, molecule);
		} else if (oxygen != null) {
			molecule.add(methane);
			UtilityReactions.addBond(carbon, oxygen, molecule);
		}
	}

}
