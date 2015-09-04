package ca.mcmaster.magarveylab.prism.cluster.reactions.typeII;

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.TypeIIPolyketideDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

public class AminotransferaseReaction extends GenericReaction implements Reaction {

	public AminotransferaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) { 
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { TypeIIPolyketideDomains.C2AMT };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		Residue residue = scaffold.residue(plan.get(0));
		IAtom alphaCarbon = residue.alphaCarbon();

		// create methyl group
		IAtomContainer nitrogen = new AtomContainer();
		IAtom n = new Atom("N");
		n.setImplicitHydrogenCount(2);
		nitrogen.addAtom(n);
	
		// add to scaffold & add bond
		if (molecule.getConnectedBondsCount(alphaCarbon) < 4) {
			molecule.add(nitrogen);
			residue.structure().add(nitrogen);
			UtilityReactions.addBond(n, alphaCarbon, molecule);
		} else {
			throw new ScaffoldGenerationException("Error: tried to aminate alpha carbon with 4 or more bonds!");
		}
	}
	
}
