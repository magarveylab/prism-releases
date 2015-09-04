package ca.mcmaster.magarveylab.prism.cluster.reactions.impl;

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.config.Elements;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.enums.substrates.AdenylationSubstrates;
import ca.mcmaster.magarveylab.enums.substrates.SubstrateType;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Atoms;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.Substrate;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.data.structure.Residue;
import ca.mcmaster.magarveylab.prism.data.structure.Scaffold;
import ca.mcmaster.magarveylab.prism.util.exception.NoResidueException;
import ca.mcmaster.magarveylab.prism.util.exception.ScaffoldGenerationException;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

public class ChlorinaseReaction extends GenericReaction implements Reaction {

	public ChlorinaseReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { TailoringDomains.CHLORINATION };
	}
	
	public void execute() throws NoResidueException, TailoringSubstrateException, ScaffoldGenerationException, CDKException {
		IAtomContainer molecule = scaffold.molecule();
		Module module = plan.get(0);
		Residue residue = scaffold.residue(module);
		Domain domain = plan.domain();
		
		IAtom atom = null;
		IAtomContainer structure = residue.structure();
		System.out.println("[ChlorinaseAnnotator] Finding chlorination atom in " + module.type() + " module");
	
		if (module.type() == ModuleTypes.ACYLTRANSFERASE) {
			atom = residue.alphaCarbon();
		} else if (module.isAdenylationModule()
				&& module.type() != ModuleTypes.ACYL_ADENYLATE) {
			Domain adenylation = module.scaffold();
			Substrate s = adenylation.topSubstrate();
			SubstrateType type = s.type();
			if (type == AdenylationSubstrates.THREONINE_1
					|| type == AdenylationSubstrates.THREONINE_2)
				atom = Atoms.getMethylCarbon(structure);
			if (type == AdenylationSubstrates.HISTIDINE_1
					|| type == AdenylationSubstrates.HISTIDINE_2)
				atom = Atoms.getHistidineChlorinationAtom(structure);
			if (type == AdenylationSubstrates.LEUCINE_1
					|| type == AdenylationSubstrates.LEUCINE_2
					|| type == AdenylationSubstrates.LEUCINE_3
					|| type == AdenylationSubstrates.VALINE_1
					|| type == AdenylationSubstrates.VALINE_3
					|| type == AdenylationSubstrates.VALINE_3)
				atom = Atoms.getMethylCarbon(structure);
			if (type == AdenylationSubstrates.TYROSINE_1 
					|| type == AdenylationSubstrates.TYROSINE_2 
					|| type == AdenylationSubstrates._4_HYDROXY_PHENYLGLYCINE 
					|| type == AdenylationSubstrates.BETA_HYDROXYTYROSINE)
				atom = Atoms.getMetaCarbon(structure);
			if (type == AdenylationSubstrates._3_5_DIHYDROXYPHENYLGLYCINE)
				atom = Atoms.getAromaticCarbon(structure);
			if (type == AdenylationSubstrates.TRYPTOPHAN)
				atom = Atoms.getTryptophanChlorinationAtom(structure);
			if (type == AdenylationSubstrates.PROLINE_1
					|| type == AdenylationSubstrates.PROLINE_2
					|| type == AdenylationSubstrates.PROLINE_3
					|| type == AdenylationSubstrates.METHYL_PROLINE)
				atom = Atoms.getProlineChlorinationAtom(structure);
		} else if (module.type() == ModuleTypes.C_STARTER
				|| module.type() == ModuleTypes.ACYL_ADENYLATE) {
			atom = Atoms.getMetaCarbon(structure);
			if (domain.blastResults().size() > 0 
					&& domain.blastResults().get(0).subject().contains("proline")) 
				atom = Atoms.getProlineChlorinationAtom(structure);
		} else if (module.type() == ModuleTypes.TYPE_II_PKS 
				|| module.type() == ModuleTypes.TYPE_II_PKS_STARTER) {
			atom = residue.alphaCarbon();
			if (molecule.getBondOrderSum(atom) != 3 || molecule.getConnectedBondsCount(atom) != 2)
				throw new ScaffoldGenerationException("Error: could not chlorinate type II polyketide module"
						+ " with incorrect bonding!");
		}
		
		if (atom == null) 
			throw new ScaffoldGenerationException("Error: chlorination site is null");
		if (scaffold.molecule().getBondOrderSum(atom) == 4)
			throw new ScaffoldGenerationException("Error: tried to chlorinate a site with four bonds");

		// create methyl group
		IAtomContainer chlorine = new AtomContainer();
		IAtom cl = new Atom(Elements.CHLORINE);
		chlorine.addAtom(cl);

		// add to scaffold & add bond
		molecule.add(chlorine);
		UtilityReactions.addBond(atom, cl, molecule);
	}
	
}
