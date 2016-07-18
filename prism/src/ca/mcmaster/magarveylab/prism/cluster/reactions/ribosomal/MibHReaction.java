package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import java.util.List;

import org.openscience.cdk.Atom;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

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
 * Execute tryptophan 5-halogenation, as catalyzed by the MibH enzyme and
 * homologs.
 * 
 * @author skinnider
 *
 */
public class MibHReaction extends GenericReaction implements Reaction {

	public MibHReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.MibH };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module module = plan.get(0);
		Residue residue = scaffold.residue(module);
		IAtomContainer structure = residue.structure();

		IAtom nitrogen = residue.nitrogen();

		// get carbon 5
		IAtom carbon = null;
		for (IAtom atom1 : structure.atoms()) {
			if (atom1.getSymbol().equals("N") && atom1 != nitrogen) {
				List<IAtom> connectedAtoms = structure.getConnectedAtomsList(atom1);
				for (IAtom atom2 : connectedAtoms) {
					if (structure.getConnectedBondsCount(atom2) == 2) {
						List<IAtom> connectedAtoms2 = structure.getConnectedAtomsList(atom2);
						for (IAtom atom3 : connectedAtoms2) {
							if (atom3 != atom1) {
								List<IAtom> connectedAtoms3 = structure.getConnectedAtomsList(atom3);
								for (IAtom atom4 : connectedAtoms3) {
									if (structure.getConnectedAtomsCount(atom4) == 3) {
										List<IAtom> connectedAtoms4 = structure.getConnectedAtomsList(atom4);
										for (IAtom atom5 : connectedAtoms4) {
											if (structure.getConnectedAtomsCount(atom5) == 2) {
												List<IAtom> connectedAtoms5 = structure.getConnectedAtomsList(atom5);
												for (IAtom atom6 : connectedAtoms5) {
													if (atom6 != atom4)
														carbon = atom6;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		if (carbon == null)
			throw new TailoringSubstrateException("Error: could not get carbon 5 for tryptophan halogenation!");

		IAtom chlorine = new Atom("Cl");
		structure.addAtom(chlorine);
		molecule.addAtom(chlorine);
		UtilityReactions.addBond(chlorine, carbon, structure);
		UtilityReactions.addBond(chlorine, carbon, molecule);
	}
	
}
