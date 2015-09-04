package ca.mcmaster.magarveylab.prism.cluster.reactions.impl;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;

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
import ca.mcmaster.magarveylab.prism.util.SmilesIO;

public class IsopenicillinNAcyltransferaseReaction extends GenericReaction implements Reaction {

	public IsopenicillinNAcyltransferaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { BetaLactamDomains.IAT };
	}
	
	public void execute() throws NoResidueException,
			TailoringSubstrateException, ScaffoldGenerationException,
			CDKException {
		IAtomContainer molecule = scaffold.molecule();
		
		Module module = plan.get(0);
		Residue residue = scaffold.residue(module);

		IAtomContainer adipicAcid = residue.structure();
		IAtom nitrogen = Atoms.getConnectedNitrogen(residue.ketone(), scaffold.molecule());

		// generate phenylacetate
		String smiles = "O=C(I)Cc1ccccc1";
		IAtomContainer phenylacetate = null;
		try {
			phenylacetate = SmilesIO.molecule(smiles);
 		} catch (Exception e) {
 			throw new ScaffoldGenerationException("Error: could not generate "
 					+ "phenylacetate for isopenicillin N acyltransferase!");
 		}
		IAtom iodine = Atoms.getIodine(phenylacetate);
		IAtom ketone = Atoms.getConnectedCarbon(iodine, phenylacetate);	
		UtilityReactions.removeIodine(phenylacetate);
		
		// remove 2-amino-adipic acid
		for (IBond bond : molecule.getConnectedBondsList(nitrogen)) {
			IAtom atom = bond.getConnectedAtom(nitrogen);
			if (atom.getSymbol().equals("C") && atom.getHybridization() == IAtomType.Hybridization.SP2)
				molecule.removeBond(bond);
		}
		molecule.remove(adipicAcid);
		System.out.println("[TailoringReactions] Removed amino adipic acid, new SMILES: " 
				+ SmilesIO.smiles(molecule));

		// add new bond
		molecule.add(phenylacetate);
		UtilityReactions.addBond(nitrogen, ketone, molecule);
	}
	
}
