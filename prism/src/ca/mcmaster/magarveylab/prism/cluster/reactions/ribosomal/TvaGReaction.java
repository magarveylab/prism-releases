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
 * Execute N1,N3-dimethylation catalyzed by the thioviridamide TvaG enzyme. 
 * 
 * @author skinnider
 *
 */
public class TvaGReaction extends GenericReaction implements Reaction {

	public TvaGReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.TvaG };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		Module module = plan.get(0);
		Residue residue = scaffold.residue(module);
		IAtomContainer structure = residue.structure();
		IAtomContainer molecule = scaffold.molecule();
		IAtom nitrogen = residue.nitrogen();
		for (IAtom atom : structure.atoms())
			if (atom.getSymbol().equals("N") && atom != nitrogen) {
				IAtom carbon = new Atom("C");
				structure.addAtom(carbon);
				molecule.addAtom(carbon);
				UtilityReactions.addBond(carbon, atom, molecule);
			}

		// set formal charge
		IAtom alphaCarbon = residue.alphaCarbon();
		IAtom betaCarbon = RibosomalUtil.getBetaCarbon(residue, structure);
		IAtom gammaCarbon = null;
		for (IAtom atom : molecule.getConnectedAtomsList(betaCarbon))
			if (atom.getSymbol().equals("C") && atom != alphaCarbon)
				gammaCarbon = atom;
		if (gammaCarbon == null)
			throw new TailoringSubstrateException("Error: could not locate gamma carbon for TvaG reaction!");
		for (IAtom atom : structure.atoms())
			if (atom.getSymbol().equals("N") && atom != nitrogen
					&& molecule.getBond(atom, gammaCarbon	) == null)
				atom.setFormalCharge(1);
	}
	
}
