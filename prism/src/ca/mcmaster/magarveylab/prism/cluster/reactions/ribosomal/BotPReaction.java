package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IAtomContainerSet;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.impl.NitroreductaseReaction;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

/**
 * Execute oxidative C-terminal decarboxylation and thiazole formation as
 * catalyzed by the bottromycin-family enzyme BotP.
 * 
 * @author skinnider
 *
 */
public class BotPReaction extends GenericReaction implements Reaction {

	public BotPReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.BotP };
	}

	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {	
		IAtomContainer molecule = scaffold.molecule();
		
		Module module = plan.get(0);
		Residue residue = scaffold.residue(module);
		
		// create azole
		SubstrateSet substrate = new SubstrateSet(module);
		ReactionPlan newPlan = new ReactionPlan(plan.domain(), substrate, plan.reaction());
		NitroreductaseReaction reaction = new NitroreductaseReaction(newPlan, scaffold, cluster);
		reaction.execute();

		// oxidatively decarboxylate
		IAtom ketone = residue.ketone();
		IAtom alphaCarbon = residue.alphaCarbon();
		for (IBond bond : molecule.getConnectedBondsList(ketone)) {
			IAtom atom = bond.getConnectedAtom(ketone);
			if (atom != alphaCarbon) {
				molecule.removeBond(bond);
				molecule.removeAtom(atom);
			} 
		}
		IBond bond = molecule.getBond(ketone, alphaCarbon);
		molecule.removeBond(bond);
		molecule.removeAtom(ketone);
		
		// partition container 
		IAtomContainerSet molecules = ConnectivityChecker.partitionIntoMolecules(molecule);
		IAtomContainer newContainer = null;
		for (IAtomContainer m : molecules.atomContainers())
			if (newContainer == null || m.getAtomCount() > newContainer.getAtomCount())
				newContainer = m;
		scaffold.setMolecule(newContainer);
	}

}
