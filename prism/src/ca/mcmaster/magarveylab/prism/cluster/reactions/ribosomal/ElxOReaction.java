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
 * Convert an N-terminal pyruvate or 2-oxobutyrate residue to lactate or
 * 2-hydroxybutyrate, as catalyzed by the epilancin short-chain dehydrogenase
 * ElxO.
 * 
 * @author skinnider
 *
 */
public class ElxOReaction extends GenericReaction implements Reaction {

	public ElxOReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.ElxO };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		Module module = plan.get(0);
		Residue residue = scaffold.residue(module);
		IAtomContainer molecule = scaffold.molecule();
		
		// must have already been transformed to pyruvate or 2-oxobutyrate
		IAtom nitrogen = residue.nitrogen();
		if (molecule.getAtomNumber(nitrogen) != -1)
			throw new TailoringSubstrateException("Error: could not execute ElxO reaction on "
					+ "N-terminal Dha/Dhb which has not been converted to pyruvate/2-oxobutyrate");
		
		// get 2-ketone
		IAtom alphaCarbon = residue.alphaCarbon();
		IAtom oxygen = Atoms.getConnectedOxygen(alphaCarbon, molecule);
	
		// reduce bond
		UtilityReactions.setBondOrder(alphaCarbon, oxygen, molecule, IBond.Order.SINGLE);
	}

}
