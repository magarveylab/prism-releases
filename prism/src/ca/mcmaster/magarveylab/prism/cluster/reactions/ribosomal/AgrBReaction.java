package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;
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
 * Execute ester or thioester bond formation between the C-terminal residue and
 * a conserved serine or cysteine residue at position -5 in auto-inducing
 * peptide biosynthesis.
 * 
 * @author skinnider
 *
 */
public class AgrBReaction extends GenericReaction implements Reaction {

	public AgrBReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.AgrB };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module m1 = plan.get(0); // C-terminus
		Module m2 = plan.get(1); // cysteine/threonine
		Residue r1 = scaffold.residue(m1);
		Residue r2 = scaffold.residue(m2);
		IAtomContainer s2 = r2.structure();
		
		// get cysteine -SH or serine -OH
		IAtom cyclization = null;
		if (m2.scaffold().topSubstrate().type() == ProteinogenicAminoAcids.CYSTEINE) {
			cyclization = Atoms.getSulfur(s2);
		} else if (m2.scaffold().topSubstrate().type() == ProteinogenicAminoAcids.SERINE) {
			cyclization = RibosomalUtil.getSerineOrThreonineHydroxyl(r2, s2, molecule);
		}
		if (cyclization == null)
			throw new TailoringSubstrateException("Could not get serine hydroxyl or cysteine sulfhydryl atom!");
		if (molecule.getConnectedBondsCount(cyclization) != 1)
			throw new TailoringSubstrateException("Serine hydroxyl or cysteine sulfhydryl has >1 bond!");
		
		// get C-terminal ketone and remove -COOH alcohol
		IAtom ketone = r1.ketone();
		UtilityReactions.removeCarboxylAlcohol(ketone, molecule);
		
		// form macrocycle 
		UtilityReactions.addBond(cyclization, ketone, molecule);
	}
	
}
