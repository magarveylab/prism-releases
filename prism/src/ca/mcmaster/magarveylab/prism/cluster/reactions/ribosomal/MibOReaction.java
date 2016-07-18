package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.Atom;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.RibosomalUtil;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

/**
 * Execute proline 3,4-dihydroxylation, as catalyzed by the MibO enzyme and
 * homologs.
 * 
 * @author skinnider
 *
 */
public class MibOReaction extends GenericReaction implements Reaction {

	public MibOReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.MibO };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module module = plan.get(0);
		Residue residue = scaffold.residue(module);
		IAtomContainer structure = residue.structure();
		IAtom alphaCarbon = residue.alphaCarbon();
		IAtom betaCarbon = RibosomalUtil.getBetaCarbon(residue, structure);
				
		// get gamma carbon 
		IAtom gammaCarbon = null;
		for (IAtom atom : structure.getConnectedAtomsList(betaCarbon)) 
			if (atom != alphaCarbon)
				gammaCarbon = atom; 
		
		// create hydroxyl groups 
		IAtom oh1 = new Atom("O");
		IAtom oh2 = new Atom("O");
		structure.addAtom(oh1);
		structure.addAtom(oh2);
		molecule.addAtom(oh1);
		molecule.addAtom(oh2);
		
		// add to beta, gamma carbons 
		UtilityReactions.addBond(oh1, betaCarbon, molecule);
		UtilityReactions.addBond(oh2, gammaCarbon, molecule);
		UtilityReactions.addBond(oh1, betaCarbon, structure);
		UtilityReactions.addBond(oh2, gammaCarbon, structure);
	}
	
}
