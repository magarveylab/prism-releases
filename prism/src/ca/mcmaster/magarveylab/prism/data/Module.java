package ca.mcmaster.magarveylab.prism.data;

import java.util.LinkedList;
import java.util.List;

import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.enums.substrates.AdenylationSubstrates;
import ca.mcmaster.magarveylab.prism.cluster.analysis.DomainAnalyzer;

/**
 * A generic biosynthetic module.
 * @author skinnider
 *
 */
public class Module {

	private boolean active = true;
	private boolean canExtend = true;
	private ModuleTypes type;
	private List<Domain> domains = new LinkedList<Domain>();

	/**
	 * Instantiate a new module of a given type.
	 * @param type	type of this module
	 */
	public Module(ModuleTypes type) {
		this.type = type;
	}

	/**
	 * Instantiate a new module by making a deep copy of an existing module. (Domains are not deep copied, but list
	 * structure is.)
	 * @param module
	 */
	public Module(Module module) {
		this.type = module.type();
		for (Domain domain : module.domains())
			domains.add(domain);
	}

	/**
	 * Add a domain to this module.
	 * @param domain	the domain to add
	 */
	public void add(Domain domain) {
		domains.add(domain);
	}

	/**
	 * Add domains to this module.
	 * @param domains	the domains to add
	 */
	public void add(List<Domain> domains) {
		this.domains.addAll(domains);
	}

	/**
	 * Get all domains within this module.
	 * @return	all domains in this module
	 */
	public List<Domain> domains() {
		return domains;
	}

	/**
	 * Get the type of this module
	 * @return	the module's type
	 */
	public ModuleTypes type() {
		return type;
	}

	/**
	 * Set the type of this module.
	 * @param type	the new type of this module
	 */
	public void setType(ModuleTypes type) {
		this.type = type;
	}

	/**
	 * Get the first domain of a given module.
	 * @return			module's first domain
	 */
	public Domain first() {
		return domains.get(0);
	}

	/**
	 * Get the last domain of a given module
	 * @return			module's last domain
	 */
	public Domain last() {
		return domains.get(domains.size() - 1);
	}

	/**
	 * Test whether this module contains a given domain type.
	 * @param type	the type to check for
	 * @return		true if the module contains a cluster of this type
	 */
	public boolean contains(DomainType type) {
		boolean flag = false;
		for (Domain domain : domains) {
			if (domain.type() == type)
				flag = true;
		}
		return flag;
	}

	/**
	 * Get the first index of a given domain type within this module.
	 * @param type	the type to check for
	 * @return		index of the first instance of this type of domain 
	 */
	public int firstIndexOf(DomainType type) {
		int index = -1;
		for (Domain t : domains)
			if (t.type() == type) {
				index = domains.indexOf(t);
				break;
			}
		return index;
	}

	/**
	 * Test whether this module contains a given domain.
	 * @param domain	the domain to check for
	 * @return		true if the module contains this domain
	 */
	public boolean contains(Domain d) {
		boolean contains = false;
		for (Domain domain : domains) {
			if (domain.equals(d)) {
				contains = true;
			}
		}
		return contains;
	}
	/**
	 * Get the scaffold domain in this module.
	 * @return	the A, AL, AT, or C-starter domain in this module
	 */
	public Domain scaffold() {
		Domain scaffold = null;
		if (type == ModuleTypes.ACYL_ADENYLATE) {
			for (Domain domain : domains)
				if (domain.type() == ThiotemplatedDomains.ACYL_ADENYLATING 
						|| (domain.type() == ThiotemplatedDomains.ADENYLATION 
						&& domain.topSubstrate().type() == AdenylationSubstrates.PROLINE_3)
						&& (scaffold == null || domain.score() > scaffold.score())) 
					scaffold = domain;
		} else if (isAdenylationModule()) {
			for (Domain domain : domains) 
				if (domain.type() == ThiotemplatedDomains.ADENYLATION
						&& (scaffold == null || domain.score() > scaffold.score())) 
					scaffold = domain;			
		} else if (type == ModuleTypes.ACYLTRANSFERASE) {
			for (Domain domain : domains)
				if (domain.type() == ThiotemplatedDomains.ACYLTRANSFERASE
						&& (scaffold == null || domain.score() > scaffold.score())) 
					scaffold = domain;				
		} else if (type == ModuleTypes.C_STARTER) {
			for (Domain domain : domains) 
				if (DomainAnalyzer.isCStarter(domain) 
						&& (scaffold == null || domain.score() > scaffold.score())) 
					scaffold = domain;
		} else if (type == ModuleTypes.RIBOSOMAL) {
			for (Domain domain : domains)
				if (domain.type() == RibosomalDomains.AMINO_ACID)
					scaffold = domain;
		}
		return scaffold;
	}

	/**
	 * Determine whether this module represents any type of adenylation module
	 * (adenylation, trans-adenylation, fatty acyl-AMP ligase, starter
	 * adenylation, trans-adenylation insertion, or integrated
	 * adenylation-ketoreductase).
	 * 
	 * @return true if this is an adenylation module
	 */
	public boolean isAdenylationModule() {
		return (type == ModuleTypes.ADENYLATION
				|| type == ModuleTypes.TRANS_ADENYLATION
				|| type == ModuleTypes.TRANS_ADENYLATION_INSERTION
				|| type == ModuleTypes.ACYL_ADENYLATE);
	}

	/**
	 * Determine whether this module represents any type of acyltransferase module (acyltransferase, trans-acyltransferase, 
	 * or trans-acyltransferase insertion).
	 * @return	true if this is an acyltransferase module
	 */
	public boolean isAcyltransferaseModule() {
		return (type == ModuleTypes.ACYLTRANSFERASE || type == ModuleTypes.TRANS_AT
				|| type == ModuleTypes.TRANS_AT_INSERTION);
	}

	/**
	 * Get the activity state of this module (i.e. is it implicated in biosynthesis, or biosynthetically inactive?)
	 * @return	activity of this module
	 */
	public boolean isActive() {
		return active;
	}

	/**
	 * Set the activity state of this module to false - i.e. assume it is not biosynthetically active. 
	 */
	public void inactivate() {
		active = false;
	}
	
	/**
	 * Get the extension state of this module: i.e., can it extend a growing natural product scaffold? 
	 * Fatty acids, for instance, can only start a natural product. Set, by default, to false.
	 * @return	extendability of this module 
	 */
	public boolean canExtend() {
		return canExtend;
	}
	
	/**
	 * Set the extension state of this module: i.e., can it extend a growing natural product scaffold?
	 * @param canExtend	extendability of this module 
	 */
	public void setCanExtend(boolean canExtend) {
		this.canExtend = canExtend;
	}

}
