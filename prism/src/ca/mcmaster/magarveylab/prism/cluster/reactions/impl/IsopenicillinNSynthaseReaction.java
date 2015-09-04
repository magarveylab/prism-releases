package ca.mcmaster.magarveylab.prism.cluster.reactions.impl;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;

import ca.mcmaster.magarveylab.enums.domains.BetaLactamDomains;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
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

public class IsopenicillinNSynthaseReaction extends GenericReaction implements Reaction {

	public IsopenicillinNSynthaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { BetaLactamDomains.IPNS };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		
		Module cysteineModule = plan.get(0);
		Module valineModule = plan.get(1);
		Residue cysteine = scaffold.residue(cysteineModule);
		Residue valine = scaffold.residue(valineModule);

		IAtom sulfur = Atoms.getSulfur(cysteine.structure());
		IAtom alphaCarbon = valine.alphaCarbon();

		// get cysteine beta carbon
		IAtom cysteineBetaCarbon = Atoms.getConnectedCarbon(sulfur, molecule);
		
		// get valine beta carbon
		IAtom valineBetaCarbon = null;
		for (IAtom atom : molecule.getConnectedAtomsList(alphaCarbon)) {
			if (atom.getSymbol().equals("C") && atom.getHybridization() == IAtomType.Hybridization.SP3)
				valineBetaCarbon = atom;
		}
		
		// get valine nitrogen
		IAtom valineNitrogen = Atoms.getConnectedNitrogen(alphaCarbon, molecule);
		
		// add cysteine beta carbon-nitrogen bond
		UtilityReactions.addBond(cysteineBetaCarbon, valineNitrogen, molecule);
		
		// add sulfur-valine beta carbon bond
		UtilityReactions.addBond(valineBetaCarbon, sulfur, molecule);
	}
	
}
