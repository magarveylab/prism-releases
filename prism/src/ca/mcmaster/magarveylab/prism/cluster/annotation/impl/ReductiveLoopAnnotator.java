package ca.mcmaster.magarveylab.prism.cluster.annotation.impl;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.enums.substrates.AcylAdenylatingSubstrates;
import ca.mcmaster.magarveylab.enums.substrates.AcyltransferaseSubstrates;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;

/**
 * Find sites at which enzyme domains of the polyketide reductive loop (KS, DH,
 * ER) could react.
 * 
 * @author skinnider
 *
 */
public class ReductiveLoopAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { ThiotemplatedDomains.KETOREDUCTASE,
				ThiotemplatedDomains.ENOYLREDUCTASE,
				ThiotemplatedDomains.DEHYDRATASE };
	}

	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster) 
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>(); 
		
		Module module = null;
		for (Module m : permutation)
			if (m.contains(domain))
				module = m;
		
		if (module != null) {
			// first, check some edge cases
			if (domain.type() == ThiotemplatedDomains.DEHYDRATASE 
					&& !module.contains(ThiotemplatedDomains.KETOREDUCTASE)) {
				// cannot have DH without KR
				return substrates;
			} else if (domain.type() == ThiotemplatedDomains.ENOYLREDUCTASE
					&& (!module.contains(ThiotemplatedDomains.KETOREDUCTASE) 
					|| !module.contains(ThiotemplatedDomains.DEHYDRATASE))) {
				// cannot have ER without KR, DH
				return substrates;
			}
			
			if (module.isAcyltransferaseModule()) {
				// the domain is in a module; target the previous domain
				int currentIdx = permutation.indexOf(module);
				int lastIdx = currentIdx - 1;
				if (lastIdx < 0)
					return substrates;
				Module last = permutation.get(lastIdx);
				
				// cannot be inactive
				if (!last.isActive())
					return substrates;

				// cannot have C-MT and no ER if DH on non-malonyl (so ignore)
				if (last.scaffold() != null 
						&& last.isAcyltransferaseModule()
						&& !(last.scaffold().topSubstrate().type() == AcyltransferaseSubstrates.MALONYL_COA_1
						|| last.scaffold().topSubstrate().type() == AcyltransferaseSubstrates.MALONYL_COA_2)
						&& module.contains(ThiotemplatedDomains.C_METHYLTRANSFERASE) 
						&& !module.contains(ThiotemplatedDomains.ENOYLREDUCTASE))
					return substrates;
				
				SubstrateSet substrate = new SubstrateSet(last);
				substrates.add(substrate);
			} else {
				if (module.type() == ModuleTypes.ACYL_ADENYLATE) {
					// if the module is actually an alpha-ketoacid, then the reaction has already been performed 
					Domain scaffold = module.scaffold();
					if (scaffold == null || scaffold.topSubstrate().type() == null)
						return substrates;
					if (scaffold.topSubstrate().type() == AcylAdenylatingSubstrates.ALPHA_KETOISOCAPROATE
							|| scaffold.topSubstrate().type() == AcylAdenylatingSubstrates.ALPHA_KETOISOVALERATE
							|| scaffold.topSubstrate().type() == AcylAdenylatingSubstrates.PYRUVATE
							|| scaffold.topSubstrate().type() == AcylAdenylatingSubstrates.PHENYLPYRUVATE
							|| scaffold.topSubstrate().type() == AcylAdenylatingSubstrates._3_METHYL_2_OXOPENTANOIC_ACID)
						return substrates;
				}
				SubstrateSet substrate = new SubstrateSet(module);
				substrates.add(substrate);
			}
		} else {
			System.out.println("[ModularAnnotator] Couldn't find module for domain " + domain.name());
			// the domain is not in a module
			//XXX Could we set up a rule to identify where it reacts, then? 
		}
		
		return substrates;
	}

}
