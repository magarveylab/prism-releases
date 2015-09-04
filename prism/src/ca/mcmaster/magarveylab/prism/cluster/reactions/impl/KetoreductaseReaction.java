package ca.mcmaster.magarveylab.prism.cluster.reactions.impl;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

public class KetoreductaseReaction extends GenericReaction implements Reaction {

	public KetoreductaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { ThiotemplatedDomains.KETOREDUCTASE };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		
		Residue residue = scaffold.residue(plan.get(0));
		if (residue == null)
			throw new TailoringSubstrateException("Error: could not get residue for ketoreductase!");
		IAtom ketone = residue.ketone();
		if (ketone == null) 
			throw new TailoringSubstrateException("Error: could not get ketone for ketoreductase!");
		
		boolean reduced = false;
		for (IBond bond : molecule.getConnectedBondsList(ketone)) {
			IAtom atom = bond.getConnectedAtom(ketone);
			if (atom.getSymbol().equals("O") && bond.getOrder() == IBond.Order.DOUBLE) {
				bond.setOrder(IBond.Order.SINGLE);
				reduced = true;
			}
		}
		
		if (!reduced)
			throw new ScaffoldGenerationException("Error: could not reduce ketone!");
	}

}
