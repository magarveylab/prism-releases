package ca.mcmaster.magarveylab.prism.data;

import java.util.LinkedList;
import java.util.List;

import ca.mcmaster.magarveylab.enums.DomainFamilies;
import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.OrfTypes;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.prism.cluster.analysis.DomainAnalyzer;

/**
 * A putative biosynthetic open reading frame.
 * 
 * @author skinnider
 *
 */
public class Orf {

	private int start;
	private int end;
	private String name;
	private String frame;
	private String dnaSequence;
	private String aaSequence;
	private List<Domain> domains = new LinkedList<Domain>();
	private List<Module> modules = new LinkedList<Module>();
	private OrfTypes type = null;
		
	/**
	 * Instantiate a new biosynthetic orf from a protein sequence.
	 * 
	 * @param name
	 *            name of the orf
	 * @param aaSequence
	 *            amino acid sequence of the orf
	 */
	public Orf(String name, String aaSequence) {
		this.name = name;
		this.aaSequence = aaSequence;
	}
	
	/**
	 * Instantiate a new orf by making a deep copy of an existing orf. Modules
	 * within the orf will be deep-copied, but not domains.
	 * 
	 * @param orf
	 *            orf to copy
	 */
	public Orf(Orf orf) {
		this.name = orf.name();
		this.start = orf.start();
		this.end = orf.end();
		this.frame = orf.frame();
		this.modules = new LinkedList<Module>();
		for (Module module : orf.modules()) {
			Module copy = new Module(module);
			modules.add(copy);
		}
		domains.addAll(orf.domains());
	}
	
	/**
	 * Get a string representation of this orf, in the format
	 * {@code start|end|frame}.
	 */
	@Override
	public String toString() {
		return start + "|" + end + "|" + frame;
	}

	/**
	 * Get the name of this orf.
	 * 
	 * @return the orf's name
	 */
	public String name() {
		return name;
	}
	
	/**
	 * Get the start point of this orf.
	 * 
	 * @return the index of the orf's start within a parent contig
	 */
	public int start() {
		return start;
	}
	
	/**
	 * Set the start point of this orf.
	 * 
	 * @param start
	 *            the index of the orf's start within a parent contig
	 */
	public void setStart(int start) {
		this.start = start;
	}
	
	/**
	 * Get the end point of this orf.
	 * 
	 * @return the index of the orf's end within a parent contig
	 */
	public int end() {
		return end;
	}

	/**
	 * Set the end point of this orf.
	 * 
	 * @param end
	 *            the index of the orf's end within a parent contig
	 */
	public void setEnd(int end) {
		this.end = end;
	}

	/**
	 * Get the length of this orf.
	 * 
	 * @return the orf's length, equivalent to (start - end) + 1
	 */
	public int length() {
		return Math.abs((start - end) + 1);
	}

	/**
	 * Get the DNA sequence of this orf.
	 * 
	 * @return the dna sequence of the orf
	 */
	public String dnaSequence() {
		return dnaSequence;
	}

	/**
	 * Get the amino acid sequence of this orf.
	 * 
	 * @return the amino acid sequence of the orf
	 */
	public String sequence() {
		return aaSequence;
	}

	/**
	 * Get all domains associated with this orf.
	 * 
	 * @return all orf domains
	 */
	public List<Domain> domains() {
		return domains;
	}

	/**
	 * Add a domain to this orf.
	 * 
	 * @param domain
	 *            domain to add
	 */
	public void add(Domain domain) {
		domains.add(domain);
	}

	/**
	 * Add domains to this orf.
	 * 
	 * @param domains
	 *            domains to add
	 */
	public void add(List<Domain> domains) {
		this.domains.addAll(domains);
	}

	/**
	 * Check whether this open reading frame has any hits to an enzymatic
	 * domain.
	 * 
	 * @return true if this orf contains a domain
	 */
	public boolean hasDomains() {
		return (domains.size() > 0);
	}

	/**
	 * Check whether this open reading frame has domains associated with it
	 * which should extend the putative cluster, i.e. which are associated with
	 * biosynthesis. Currently, resistance and regulatory genes do not extend a
	 * cluster.
	 * 
	 * @see ca.mcmaster.magarveylab.prism.cluster.ClusterFinder
	 * @return true if clusterable domains is larger than 0
	 */
	public boolean hasBiosyntheticDomains() {
		List<Domain> resistance = domains(DomainFamilies.RESISTANCE);
		List<Domain> regulator = domains(DomainFamilies.REGULATOR);
		return domains.size() > (resistance.size() + regulator.size());
	}

	/**
	 * Check whether this open reading frame contains a given domain.
	 * 
	 * @param domain
	 *            domain to test
	 * @return true if the open reading frame contains this particular domain
	 */
	public boolean contains(Domain domain) {
		return (domains.contains(domain));
	}

	/**
	 * Get all domains of a given family within this orf.
	 * 
	 * @param family
	 *            family to get
	 * @return all domains from that family in this orf
	 */
	public List<Domain> domains(DomainFamilies family) {
		List<Domain> familyDomains = new LinkedList<Domain>();
		for (Domain domain : domains)
			if (domain.family() == family)
				familyDomains.add(domain);
		return familyDomains;
	}

	/**
	 * Get all domains of a given type on this orf.
	 * 
	 * @param type
	 *            the type of domain to return
	 * @return all domains of this type
	 */
	public List<Domain> domains(DomainType type) {
		return DomainAnalyzer.domains(type, domains);
	}

	/**
	 * Get all biosynthetic modules associated with this orf.
	 * 
	 * @return all biosynthetic modules
	 */
	public List<Module> modules() {
		return modules;
	}

	/**
	 * Get all biosynthetic modules associated with this orf which have not been
	 * marked as biosynthetically inactive.
	 * 
	 * @return all active biosynthetic modules
	 */
	public List<Module> activeModules() {
		List<Module> active = new LinkedList<Module>();
		for (Module module : modules)
			if (module.isActive())
				active.add(module);
		return active;
	}

	/**
	 * Get the last module on this open reading frame.
	 * 
	 * @return the last module
	 */
	public Module getLastModule() {
		return modules.get(modules.size() - 1);
	}

	/**
	 * Set the biosynthetic modules present on this orf.
	 * 
	 * @param modules
	 *            list of biosynthetic modules associated with this orf
	 */
	public void setModules(List<Module> modules) {
		this.modules = modules;
	}

	/**
	 * Get the type of this orf.
	 * 
	 * @return the orf's type
	 */
	public OrfTypes type() {
		return type;
	}
	
	/**
	 * Set the type of this orf (PKS, NRPS, hybrid, modification, inactive, or null).
	 * @param type	the orf's type
	 */
	public void setType(OrfTypes type) {
		this.type = type;
	}
	
	/**
	 * Get the frame of this orf. Can either be '+' or '-'. 
	 * @return	the orf's frame
	 */
	public String frame() {
		return frame;
	}

	/**
	 * Set the frame of this orf.
	 * @param frame	the orf's frame
	 */
	public void setFrame(String frame) {
		this.frame = frame;
	}
	
	/**
	 * Get all modules of a given type within this orf.
	 * @param type		type of module to find
	 * @return			all modules of that type
	 */
	public List<Module> modules(ModuleTypes type) {
		List<Module> modules = new LinkedList<Module>();
		for (Module module : this.modules)
			if (module.type() == type)
				modules.add(module);
		return modules;
	}

	/**
	 * Determine whether this orf contains a given module.
	 * @param module	module to query
	 * @return			true if orf contains module
	 */
	public boolean contains(Module module) {
		boolean flag = false;
		for (Module m : modules) {
			if (m == module) 
				flag = true;
		}
		return flag;
	}

	/**
	 * Determine whether this orf contains a given module type.
	 * 
	 * @param module
	 *            module type to query
	 * @return true if orf contains module type
	 */
	public boolean contains(ModuleTypes module) {
		boolean flag = false;
		for (Module m : modules) {
			if (m.type() == module)
				flag = true;
		}
		return flag;
	}

	/**
	 * Replace one module with another in this orf.
	 * 
	 * @param module1
	 *            module to be replaced
	 * @param module2
	 *            module to replace with
	 */
	public void replace(Module module1, Module module2) {
		for (int i = 0; i < modules.size(); i++) {
			Module module = modules.get(i);
			if (module == module1) {
				modules.set(i, module2);
			}
		}
	}

	/**
	 * Determine whether this orf contains an active biosynthetic module.
	 * 
	 * @return true if this orf contains an active module
	 */
	public boolean hasActiveModule() {
		boolean flag = false;
		for (Module module : modules)
			if (module.isActive())
				flag = true;
		return flag;
	}

	/**
	 * Get all domains with substrates on this orf (i.e., adenylation,
	 * acyltransferase, and acyl-adenylating).
	 * 
	 * @return all substrate-containing domains
	 */
	public List<Domain> substrateDomains() {
		List<Domain> substrateDomains = new LinkedList<Domain>();
		substrateDomains.addAll(domains(ThiotemplatedDomains.ADENYLATION));
		substrateDomains.addAll(domains(ThiotemplatedDomains.ACYLTRANSFERASE));
		substrateDomains.addAll(domains(ThiotemplatedDomains.ACYL_ADENYLATING));
		return substrateDomains;
	}
	
}
