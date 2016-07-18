package ca.mcmaster.magarveylab.prism.data;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import ca.mcmaster.magarveylab.enums.DomainFamilies;
import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.OrfTypes;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.prism.cluster.analysis.DomainAnalyzer;
import ca.mcmaster.magarveylab.prism.orfs.GenePredictionModes;

/**
 * A putative biosynthetic open reading frame.
 * 
 * @author skinnider
 *
 */
public class Orf implements Serializable {

	private static final long serialVersionUID = -5713920469341661161L;
	
	private int start;
	private int end;
	private String name;
	private String frame;
	private String dnaSequence;
	private String aaSequence;
	private List<Domain> domains = new LinkedList<Domain>();
	private List<Module> modules = new LinkedList<Module>();
	private List<Propeptide> propeptides = new ArrayList<Propeptide>(); 
	private OrfTypes type = null;
	private GenePredictionModes mode = null;
	private String partial;
		
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
		return Math.abs(start - end) + 1;
	}

	/**
	 * Get the DNA sequence of this orf.
	 * 
	 * @return the dna sequence of the orf
	 */
	public String dnaSequence() {
		return dnaSequence;
	}
	
	public void setDnaSequence(String seq) {
		this.dnaSequence = seq;
	}
	
	public void setAASequence(String seq) {
		this.aaSequence = seq;
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
	
	/**
	 * Get all cleaved ribosomal propeptides associated with this orf. 
	 * @return	all cleaved ribosomal propeptides associated with this orf
	 */
	public List<Propeptide> propeptides() {
		return propeptides;
	}
	
	/**
	 * Check whether this orf has cleaved ribosomal propeptides associated with it.
	 * @return true if this orf has cleaved ribosomal propeptides associated with it
	 */
	public boolean hasPropeptides() {
		return propeptides.size() > 0;
	}
	
	/**
	 * Associate a new cleaved, ribosomal propeptide with this orf. 
	 * @param propeptide	a cleaved ribosomal propeptide from this orf 
	 */
	public void addPropeptide(Propeptide propeptide) {
		propeptides.add(propeptide);
	}
	
	/**
	 * Associate a list of cleaved, ribosomal propeptides with this orf. 
	 * @param propeptide	a list of cleaved ribosomal propeptides from this orf 
	 */
	public void addPropeptides(List<Propeptide> propeptides) {
		this.propeptides.addAll(propeptides);
	}

	/**
	 * Get the mode by which this orf was identified or predicted (for example,
	 * was this orf identified by scanning the input sequence for all possible
	 * coding sequences, or with an orf prediction software such as Prodigal or
	 * GeneMark?)
	 * 
	 * @return	the mode by which this orf was identified or predicted 
	 */
	public GenePredictionModes getMode() {
		return mode;
	}

	/**
	 * Set the mode by which this orf was identified or predicted (for example,
	 * by scanning the input sequence for all possible coding sequences, or with
	 * an orf prediction software such as Prodigal or GeneMark)
	 * 
	 * @param mode
	 *            the mode by which this orf was identified or predicted
	 */
	public void setMode(GenePredictionModes mode) {
		this.mode = mode;
	}
	
	/**
	 * Check whether another orf overlaps with this one.
	 * 
	 * @param o2
	 *            the second orf
	 * @return true if their start-end ranges overlap
	 */
	public boolean overlaps(final Orf orf2) {
		final int x1 = start;
		final int x2 = end;
		final int y1 = orf2.start();
		final int y2 = orf2.end();
		return (x1 <= y2 && y1 <= x2);
	}

	/**
	 * For open reading frames predicted by Prodigal, get the partial indicator.
	 * This is a string of length 2 where a value of 1 means the orf is partial
	 * at that boundary: e.g. the string "01" indicates the orf is partial at
	 * the right boundary.
	 * 
	 * @return the partial indicator set by Prodigal
	 */
	public String getPartial() {
		return partial;
	}

	/**
	 * For open reading frames predicted by Prodigal, set the partial indicator.
	 * This is a string of length 2 where a value of 1 means the orf is partial
	 * at that boundary: e.g. the string "01" indicates the orf is partial at
	 * the right boundary.
	 * 
	 * @param partial
	 *            the partial indicator set by Prodigal
	 */
	public void setPartial(String partial) {
		this.partial = partial;
	}
	
}
