package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

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
 * Execute the attachment of a modified quinaldic acid to a thiopeptide
 * scaffold, as in the thiostrepton and siomycin clusters.
 * 
 * @author skinnider
 *
 */
public class TsrIReaction extends GenericReaction implements Reaction {

	public TsrIReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.TsrI };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module m1 = plan.get(0); // first residue 
		Module m2 = plan.get(1); // ser/thr residue 
		Residue r1 = scaffold.residue(m1);
		Residue r2 = scaffold.residue(m2);
		if (r1 == null)
			throw new TailoringSubstrateException("Could not get first residue for TsrI reaction!");
		if (r2 == null)
			throw new TailoringSubstrateException("Could not get second residue for TsrI reaction!");
		IAtomContainer s2 = r2.structure();

		// check first residue
		if (RibosomalUtil.isDehydrated(r2, s2, molecule))
			throw new TailoringSubstrateException("Could not esterify a dehydrated serine/threonine residue!");
		
		IAtom hydroxyl = RibosomalUtil.getSerineOrThreonineHydroxyl(r2, s2, molecule);
		if (hydroxyl == null)
			throw new TailoringSubstrateException("Could not get esterification hydroxyl for TsrI reaction!");
		
		// add the indolic acid
		String smiles = "O=C(C1=NC2=C(C(C(C)O)=C1)C=CC(F)C2O)I";
		UtilityReactions.functionalize(smiles, hydroxyl, molecule);
		
		// get the nitrogen
		IAtom nitrogen = r1.nitrogen();

		// get amination site
		IAtom fluorine = Atoms.getFluorine(molecule);
		if (fluorine == null)
			throw new TailoringSubstrateException("Could not get fluorine for TsrI reaction!");
		IAtom carbon = Atoms.getConnectedCarbon(fluorine, molecule);
		if (carbon == null)
			throw new TailoringSubstrateException("Could not get quinaldic acid amination "
					+ "carbon for TsrI reaction!");

		// remove F
		UtilityReactions.removeAtom(fluorine, molecule);

		// add the last bond
		UtilityReactions.addBond(nitrogen, carbon, molecule);
	}
	
}
