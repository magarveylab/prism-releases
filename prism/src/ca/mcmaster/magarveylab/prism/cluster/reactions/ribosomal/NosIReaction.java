package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import java.util.List;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

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
 * Execute the attachment of a modified indolic acid to a thiopeptide scaffold,
 * as in the nocathiacin and nosiheptide clusters.
 * 
 * @author skinnider
 *
 */
public class NosIReaction extends GenericReaction implements Reaction {

	public NosIReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.NosI };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module m1 = plan.get(0);
		Module m2 = plan.get(1);
		Residue r1 = scaffold.residue(m1);
		Residue r2 = scaffold.residue(m2);
		if (r1 == null)
			throw new TailoringSubstrateException("Could not get first residue for NosI reaction!");
		if (r2 == null)
			throw new TailoringSubstrateException("Could not get second residue for NosI reaction!");
		IAtomContainer s1 = r1.structure();
		IAtomContainer s2 = r2.structure();

		// check first residue
		IAtom site = null;
		IAtom ketone = r1.ketone();
		List<IAtom> carboxyls = Atoms.getAllCarboxyls(molecule);
		IAtom ketoneOxygen = Atoms.getConnectedOxygen(ketone, s1);
		for (IAtom atom : s1.atoms())
			if (atom.getSymbol().equals("S") 
					&& molecule.getConnectedBondsCount(atom) == 1) {
				site = atom;
			} else if (atom.getSymbol().equals("O") 
					&& !carboxyls.contains(atom)
					&& molecule.getConnectedBondsCount(atom) == 1
					&& atom != ketoneOxygen) {
				site = atom;
			}
		if (site == null)
			throw new TailoringSubstrateException("Could not ester/thioester site for NosI reaction!");
		
		// add the indolic acid
		String smiles = "CC1=C(C(I)=O)NC2=CC=CC(COF)=C21";
		UtilityReactions.functionalize(smiles, site, molecule);
		
		// get the D/E side-chain carboxyl
		IAtom site2 = null;
		IAtom ketone2 = r2.ketone();
		carboxyls = Atoms.getAllCarboxyls(s2);
		for (IAtom atom : carboxyls)
			if (atom != ketone2)
				site2 = atom;
		if (site2 == null)
			throw new TailoringSubstrateException("Could not get Asp/Glu carboxylic acid for NosI reaction!");
		
		// remove carboxyl alcohol
		UtilityReactions.removeCarboxylAlcohol(site2, molecule);
		
		// remove F
		IAtom fluorine = Atoms.getFluorine(molecule);
		if (fluorine == null)
			throw new TailoringSubstrateException("Could not get fluorine for NosI reaction!");
		IAtom oxygen = Atoms.getConnectedOxygen(fluorine, molecule);
		if (oxygen == null)
			throw new TailoringSubstrateException("Could not get indolic acid oxygen for NosI reaction!");
		UtilityReactions.removeAtom(fluorine, molecule);
		
		// add the last bond
		UtilityReactions.addBond(site2, oxygen, molecule);
	}
	
}
