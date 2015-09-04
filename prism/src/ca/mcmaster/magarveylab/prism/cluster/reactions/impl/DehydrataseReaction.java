package ca.mcmaster.magarveylab.prism.cluster.reactions.impl;

import java.util.List;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
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

public class DehydrataseReaction extends GenericReaction implements Reaction {
	
	public DehydrataseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { ThiotemplatedDomains.DEHYDRATASE };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException {
		IAtomContainer molecule = scaffold.molecule();
		
		Residue residue = scaffold.residue(plan.get(0));
		if (residue == null) 
			throw new TailoringSubstrateException("Error: could not get residue for dehydratase!");
		IAtom ketone = residue.ketone();
		IAtom alphaCarbon = residue.alphaCarbon();
		if (ketone == null) 
			throw new TailoringSubstrateException("Error: could not get ketone for dehydratase!");
		if (alphaCarbon == null) 
			throw new TailoringSubstrateException("Error: could not get alpha carbon for dehydratase!");

		List<IBond> bonds = molecule.getConnectedBondsList(ketone);
		for (IBond bond : bonds) {  
			IAtom connectedAtom = bond.getConnectedAtom(ketone);
			if (connectedAtom.getSymbol().equals("O")) {
				molecule.removeBond(bond);
				molecule.removeAtom(connectedAtom);
			}
			if (connectedAtom != alphaCarbon && !connectedAtom.getSymbol().equals("O") 
					&& molecule.getConnectedBondsCount(connectedAtom) < 4
					&& alphaCarbon != null && alphaCarbon.getHybridization() != IAtomType.Hybridization.SP2)
				bond.setOrder(IBond.Order.DOUBLE);
		}
	}

}
