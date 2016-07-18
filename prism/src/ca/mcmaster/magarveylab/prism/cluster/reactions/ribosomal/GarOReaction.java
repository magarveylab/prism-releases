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
 * Execute lantithione S-oxidation as catalyzed by the actagardine GarO enzyme.
 * 
 * @author skinnider
 *
 */
public class GarOReaction extends GenericReaction implements Reaction {

	public GarOReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.GarO };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module m1 = plan.get(0);
		Residue dhaDhb = scaffold.residue(m1);
		if (dhaDhb == null) 
			throw new TailoringSubstrateException("Error: could not get Dha/Dhb residue for GarO!");
		IAtomContainer dhaDhbStructure = dhaDhb.structure();
		
		Module m2 = plan.get(1);
		Residue cysteine = scaffold.residue(m2);
		if (cysteine == null) 
			throw new TailoringSubstrateException("Error: could not get Cys residue for GarO!");
		IAtomContainer cysteineStructure = cysteine.structure();

		// Dhb/Dha cannot still have Ser/Thr oxygen 
		if (!RibosomalUtil.isDehydrated(dhaDhb, dhaDhbStructure, molecule))
			throw new TailoringSubstrateException("Error: serine/threonine residue "
					+ "which has not been dehydrated!");
		
		// get cysteine sulfur
		IAtom sulfur = Atoms.getSulfur(cysteineStructure);
		if (sulfur == null)
			throw new TailoringSubstrateException("Error: could not get cysteine sulfur!");
		
		// check that there is a bond between cysteine sulfur and dha/dhb beta carbon
		IAtom dhaBetaCarbon = RibosomalUtil.getBetaCarbon(dhaDhb, dhaDhbStructure);
		if (molecule.getBond(dhaBetaCarbon, sulfur) == null)
			throw new TailoringSubstrateException("Error: could not oxidize sulfur that is "
					+ "not in a lantithione bond!");
		
		// create =O
		IAtom oxygen = new Atom("O");
		molecule.addAtom(oxygen);
		cysteineStructure.addAtom(oxygen);
		
		// add double bond
		UtilityReactions.addBond(sulfur, oxygen, molecule, IBond.Order.DOUBLE);
	}
	
}
