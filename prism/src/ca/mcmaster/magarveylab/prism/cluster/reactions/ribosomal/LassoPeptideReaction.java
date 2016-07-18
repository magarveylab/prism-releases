package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import java.util.List;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
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

/**
 * Execute macrolactam formation catalyzed by an asparagine synthase homolog in
 * the lasso peptide family of RiPPs.
 * 
 * @author skinnider
 *
 */
public class LassoPeptideReaction extends GenericReaction implements Reaction {
	
	public LassoPeptideReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.Asparagine_synthase };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module m1 = plan.get(0);
		Module m2 = plan.get(1);
		Residue r1 = scaffold.residue(m1);
		Residue r2 = scaffold.residue(m2);
		IAtomContainer s2 = r2.structure();
		
		// get nitrogen
		IAtom nitrogen = r1.nitrogen();
		
		// get Asp/Glu terminal -COOH
		IAtom carboxyl = null;		
		for (IAtom a : s2.atoms()) {
			boolean hydroxyl = false, ketone = false;
			if (a.getSymbol().equals("C") && a != r2.ketone()
					&& s2.getConnectedBondsCount(a) == 3) {
				List<IAtom> atoms = s2.getConnectedAtomsList(a);
				for (IAtom atom : atoms)
					if (atom.getSymbol().equals("O")) {
						if (s2.getBond(atom, a).getOrder() == IBond.Order.DOUBLE)
							ketone = true;
						if (s2.getBond(atom, a).getOrder() == IBond.Order.SINGLE)
							hydroxyl = true;
					}
			}
			if (ketone && hydroxyl)
				carboxyl = a;
		}
		
		if (carboxyl == null)
			throw new TailoringSubstrateException("Could not identify COOH for lasso peptide reaction!");
		
		// remove carboxyl alcohol
		UtilityReactions.removeAlcohol(carboxyl, s2);
		UtilityReactions.removeAlcohol(carboxyl, molecule);
		
		// connected carboxyl carbon and N-terminal nitrogen
		UtilityReactions.addBond(nitrogen, carboxyl, molecule);
	}
	
}
