package ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal;

import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.enums.interfaces.SubstrateType;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.Reaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.RibosomalUtil;
import ca.mcmaster.magarveylab.prism.cluster.reactions.UtilityReactions;
import ca.mcmaster.magarveylab.prism.cluster.scaffold.Chemoinformatics.Atoms;
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
 * Execute the reaction catalyzed by the lantibiotic LanKC fused cyclase/
 * dehydratase, which catalyzes labionin bond formation between three residues.
 * 
 * @author skinnider
 *
 */
public class LanKCReaction extends GenericReaction implements Reaction {
	
	public LanKCReaction(ReactionPlan plan, Scaffold scaffold, Cluster cluster) {
		super(plan, scaffold, cluster);
		this.domains = new DomainType[] { RibosomalDomains.LanKC };
	}

	public void execute() throws NoResidueException,
			TailoringSubstrateException, ScaffoldGenerationException,
			CDKException {
		// split up reactions 
		SubstrateSet set = plan.modules();
		List<Module> labionins = new ArrayList<Module>();
		List<Module> lanthionines = new ArrayList<Module>();
		List<Module> dehydrations = new ArrayList<Module>();
		int spacerPos = 0;
		for (Module module : set.getAllModules()) {
			// increment spacer
			if (module == null) {
				spacerPos++;
				continue;
			}

			if (spacerPos == 0) {
				labionins.add(module);
			} else if (spacerPos == 1) {
				lanthionines.add(module);
			} else if (spacerPos == 2) {
				dehydrations.add(module);
			}
		}

		// execute labionin formation 
		executeLabionins(labionins);

		// execute lanthionine formation
		SubstrateSet s1 = new SubstrateSet(lanthionines);
		ReactionPlan lancPlan = new ReactionPlan(plan.domain(), s1,
				plan.reaction());
		LanCReaction lanc = new LanCReaction(lancPlan, scaffold, cluster);
		lanc.execute();
		
		// execute dehydrations 
		SubstrateSet s2 = new SubstrateSet(dehydrations);
		ReactionPlan lanbPlan = new ReactionPlan(plan.domain(), s2,
				plan.reaction());
		LanBReaction lanb = new LanBReaction(lanbPlan, scaffold, cluster);
		lanb.execute();
	}
	
	private void executeLabionins(List<Module> modules)
			throws ArrayIndexOutOfBoundsException, ScaffoldGenerationException,
			TailoringSubstrateException {
		IAtomContainer molecule = scaffold.molecule();
		for (int i = 0; i < modules.size() - 2; i+= 3) {
			Module m1 = modules.get(i);
			Module m2 = modules.get(i+1);
			Module m3 = modules.get(i+2);

			SubstrateType type = m2.scaffold().topSubstrate().type(); 
			if (type == ProteinogenicAminoAcids.CYSTEINE) {
				// execute labionin formation 
				Residue r1 = scaffold.residue(m1);
				if (r1 == null) 
					throw new TailoringSubstrateException("Error: could not get residue 1 for LanKC!");
				IAtomContainer s1 = r1.structure();
				IAtom k1 = r1.ketone();
				IAtom kO1 = Atoms.getConnectedOxygen(k1, s1);
				if (k1 == null) 
					throw new TailoringSubstrateException("Error: could not get ketone carbon for LanKC!");
				if (kO1 == null)
					throw new TailoringSubstrateException("Error: could not get ketone oxygen for LanKC!");
				IAtom alphaCarbon = r1.alphaCarbon();
				if (alphaCarbon == null) 
					throw new TailoringSubstrateException("Error: could not get alpha carbon for LanKC!");

				IAtom o1 = null;
				for (IAtom atom : s1.atoms())
					if (atom.getSymbol().equals("O") && atom != kO1)
						o1 = atom;
				if (o1 == null)
					throw new TailoringSubstrateException("Error: could not get serine/threonine oxygen for LanKC!");
				
				// set double bond
				IAtom betaCarbon = Atoms.getConnectedCarbon(o1, s1);
				UtilityReactions.setBondOrder(alphaCarbon, betaCarbon, molecule, IBond.Order.DOUBLE);
				UtilityReactions.removeAtom(o1, molecule);
				UtilityReactions.removeAtom(o1, s1);
				
				SubstrateSet s2 = new SubstrateSet(m1, m2);
				ReactionPlan lancPlan = new ReactionPlan(plan.domain(), s2, plan.reaction());
				LanCReaction lanc = new LanCReaction(lancPlan, scaffold, cluster);
				lanc.execute();
				
				// get first Dha/Dhb
				Residue dhaDhb = scaffold.residue(m1);
				if (dhaDhb == null) 
					throw new TailoringSubstrateException("Error: could not get first Dha/Dhb residue for LanKC!");
				
				// get second Dha/Dhb
				Residue dhaDhb2 = scaffold.residue(m3);
				if (dhaDhb2 == null) 
					throw new TailoringSubstrateException("Error: could not get second Dha/Dhb residue for LanKC!");
				IAtomContainer dhaDhbStructure2 = dhaDhb2.structure();
				
				// remove Ser/Thr oxygen from second Dha/Dhb
				IAtom oxygen = null;
				IAtom ketone = dhaDhb2.ketone();
				IAtom ketoneOxygen = Atoms.getConnectedOxygen(ketone, dhaDhbStructure2);
				if (ketone == null) 
					throw new TailoringSubstrateException("Error: could not get ketone carbon for LanKC!");
				if (ketoneOxygen == null)
					throw new TailoringSubstrateException("Error: could not get ketone oxygen for LanKC!");
				for (IAtom atom : dhaDhbStructure2.atoms())
					if (atom.getSymbol().equals("O") && atom != ketoneOxygen)
						oxygen = atom;
				if (oxygen == null)
					throw new TailoringSubstrateException("Error: could not get serine/threonine oxygen for LanKC!");
				UtilityReactions.removeAtom(oxygen, molecule);
				UtilityReactions.removeAtom(oxygen, dhaDhbStructure2);
				
				// get first Dha/Dhb alpha carbon
				IAtom alphaCarbon1 = dhaDhb.alphaCarbon();
				
				// get second Dhb/Dha beta carbon
				IAtom betaCarbon2 = RibosomalUtil.getBetaCarbon(dhaDhb2, dhaDhbStructure2);
				
				// add bond between first Dhb/Dha alpha carbon and second Dhb/Dha beta carbon
				UtilityReactions.addBond(alphaCarbon1, betaCarbon2, molecule);
				
				// reduce second Dhb/Dha
				IAtom alphaCarbon2 = dhaDhb2.alphaCarbon();
				UtilityReactions.setBondOrder(alphaCarbon2, betaCarbon2, molecule, IBond.Order.SINGLE);
			} else if (type == ProteinogenicAminoAcids.SERINE
					|| type == ProteinogenicAminoAcids.THREONINE) {
				// dehydrate 
				SubstrateSet s1 = new SubstrateSet(m1);
				ReactionPlan lanbPlan1 = new ReactionPlan(plan.domain(), s1, plan.reaction());
				LanBReaction lanb1 = new LanBReaction(lanbPlan1, scaffold, cluster);
				lanb1.execute();
				SubstrateSet s2 = new SubstrateSet(m2);
				ReactionPlan lanbPlan2 = new ReactionPlan(plan.domain(), s2, plan.reaction());
				LanBReaction lanb2 = new LanBReaction(lanbPlan2, scaffold, cluster);
				lanb2.execute();
				SubstrateSet s3 = new SubstrateSet(m3);
				ReactionPlan lanbPlan3 = new ReactionPlan(plan.domain(), s3, plan.reaction());
				LanBReaction lanb3 = new LanBReaction(lanbPlan3, scaffold, cluster);
				lanb3.execute();
			}
		}
	}
	
}
