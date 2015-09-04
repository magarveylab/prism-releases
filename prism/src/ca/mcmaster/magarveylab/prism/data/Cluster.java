package ca.mcmaster.magarveylab.prism.data;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ca.mcmaster.magarveylab.enums.ClusterFamilies;
import ca.mcmaster.magarveylab.enums.DomainFamilies;
import ca.mcmaster.magarveylab.enums.Frames;
import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.clusters.ClusterType;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.prism.data.sugar.Sugar;
import ca.mcmaster.magarveylab.prism.homology.data.HomologousCluster;
import ca.mcmaster.magarveylab.prism.tanimoto.data.TanimotoScore;
import ca.mcmaster.magarveylab.prism.util.Sorter;
import ca.mcmaster.magarveylab.prism.web.html.graph.CircularGenomeGraph;

/**
 * A biosynthetic cluster encompassing multiple orfs.
 * 
 * @author skinnider
 *
 */
public class Cluster {

	protected int index;
	protected Frames frame;
	protected List<Orf> orfs = new ArrayList<Orf>();
	protected List<ClusterFamilies> families = new ArrayList<ClusterFamilies>();
	protected List<ClusterType> types = new ArrayList<ClusterType>();
	protected List<List<Sugar>> sugars = new ArrayList<List<Sugar>>();

	private Library library = new Library();
	private CombinatorialData combinatorialData = new CombinatorialData();
	private List<HomologousCluster> homologs = new ArrayList<HomologousCluster>();
	private List<TanimotoScore> scores = new ArrayList<TanimotoScore>();
	private CircularGenomeGraph graph;
	private Map<String, String> files = new HashMap<String, String>();

	/**
	 * Instantiate a new, empty cluster.
	 */
	public Cluster() {
	}

	/**
	 * Get all orfs associated with this cluster.
	 * 
	 * @return all cluster orfs
	 */
	public List<Orf> orfs() {
		return orfs;
	}

	/**
	 * Get the name of this cluster, defined as the name of the first orf, or
	 * "empty cluster" if no orfs are present.
	 * 
	 * @return name of this cluster
	 */
	public String name() {
		if (orfs.size() > 0) {
			return orfs.get(0).toString();
		} else {
			return "empty cluster";
		}
	}

	/**
	 * Get the start index of this cluster, defined as the start index of its
	 * first orf.
	 * 
	 * @return start index of this cluster
	 */
	public int start() {
		Integer start = null;
		for (Orf orf : orfs) {
			if (start == null || orf.start() < start)
				start = orf.start();
		}
		return start;
	}
	
	/**
	 * Get the end index of this cluster, defined as the end index of its last
	 * orf.
	 * 
	 * @return end index of this cluster
	 */
	public int end() {
		Integer end = null;
		for (Orf orf : orfs) {
			if (end == null || orf.end() > end)
				end = orf.end();
		}
		return end;
	}

	/**
	 * Get the size of this cluster, defined as its end index minus its start
	 * index.
	 * 
	 * @return size of this cluster
	 */
	public int size() {
		return end() - start();
	}

	/**
	 * Get the index of this cluster.
	 * 
	 * @return index of this cluster
	 */
	public int index() {
		return index;
	}

	/**
	 * Set the index of this cluster in the genome.
	 * 
	 * @param index
	 *            index of this cluster
	 */
	public void setIndex(int index) {
		this.index = index;
	}

	/**
	 * Get the frame of this cluster.
	 * @return the cluster frame
	 */
	public Frames frame() {
		return frame;
	}

	/**
	 * Set the frame of this cluster.
	 * @param frame		the cluster frame
	 */
	public void setFrame(Frames frame) {
		this.frame = frame;
	}
	
	/**
	 * Check whether this cluster is colinear, i.e., all scaffold biosynthetic
	 * orfs are in the same frame.
	 * 
	 * @return true if all scaffold biosynthetic orfs are in the same frame,
	 *         i.e. cluster frame is either + or -
	 */
	public boolean isColinear() {
		return frame == Frames.NEGATIVE || frame == Frames.POSITIVE;
	}

	/**
	 * Get the family or families of this cluster (for circular genome graph).
	 * 
	 * @return cluster family/families
	 */
	public List<ClusterFamilies> families() {
		return families;
	}

	/**
	 * Get the family-specific subtype or subtypes of this cluster. E.g., for
	 * thiotemplated clusters, nonribosomal peptide and/or polyketide.
	 * 
	 * @return cluster type(s)
	 */
	public List<ClusterType> types() {
		return types;
	}

	/**
	 * Add a family-specific subtype to this cluster.
	 * 
	 * @param type
	 *            a biosynthetic family subtype to which this cluster belongs
	 */
	public void addType(ClusterType type) {
		types.add(type);
	}

	/**
	 * Add a list of family-specific subtypes to this cluster.
	 * 
	 * @param types
	 *            a list of biosynthetic family subtypes to which this cluster
	 *            belongs
	 */
	public void addTypes(List<ClusterType> types) {
		this.types.addAll(types);
	}
	
	/**
	 * Add a family to this cluster.
	 * 
	 * @param family
	 *            cluster family
	 */
	public void addFamily(ClusterFamilies family) {
		families.add(family);
	}

	/**
	 * Get the genome graph associated with this cluster (this cluster
	 * highlighted).
	 * 
	 * @return this cluster's genome graph
	 */
	public CircularGenomeGraph graph() {
		return graph;
	}

	/**
	 * Set the genome graph associated with this cluster (this cluster
	 * highlighted).
	 * 
	 * @param graph
	 *            this cluster's genome graph
	 */
	public void setGraph(CircularGenomeGraph graph) {
		this.graph = graph;
	}

	/**
	 * Get all known biosynthetic clusters with genetic homology to this
	 * cluster.
	 * 
	 * @return all homologous clusters
	 */
	public List<HomologousCluster> homologs() {
		return homologs;
	}
	
	/**
	 * Set the list of known biosynthetic clusters with genetic homology to this
	 * cluster
	 * 
	 * @param homologs
	 *            list of homologous clusters
	 */
	public void setHomologs(List<HomologousCluster> homologs) {
		this.homologs = homologs;
	}

	/**
	 * Get the combinatorial scaffold library associated with this cluster.
	 * 
	 * @return the cluster's scaffold library
	 */
	public Library library() {
		return library;
	}

	/**
	 * Set the combinatorial scaffold library associated with this cluster.
	 * 
	 * @param library
	 *            the cluster's scaffold library
	 */
	public void setLibrary(Library library) {
		this.library = library;
	}

	/**
	 * Get a file associated with this cluster.
	 * 
	 * @param name
	 *            key of the desired file
	 * @return file path associated with that key
	 */
	public String file(String name) {
		return files.get(name);
	}

	/**
	 * Add a file to this cluster.
	 * 
	 * @param name
	 *            key of the new file
	 * @param path
	 *            path of the new file
	 */
	public void addFile(String name, String path) {
		files.put(name, path);
	}

	/**
	 * Test whether this cluster contains at least one domain of a given type.
	 * 
	 * @param type
	 *            the type of domain
	 * @return true if this cluster contains at least one such domain
	 */
	public boolean contains(DomainType type) {
		boolean flag = false;
		for (Orf orf : orfs)
			if (orf.domains(type).size() > 0)
				flag = true;
		return flag;
	}

	/**
	 * Get all modules of a given type within this cluster.
	 * 
	 * @param type
	 *            type of module to find
	 * @return all modules of that type
	 */
	public List<Module> modules(ModuleTypes type) {
		List<Module> modules = new ArrayList<Module>();
		for (Orf orf : orfs) {
			for (Module module : orf.modules())
				if (module.type() == type)
					modules.add(module);
		}
		return modules;
	}

	/**
	 * Get all active modules of a given type within this cluster.
	 * 
	 * @param type
	 *            type of module to find
	 * @return all active modules of that type
	 */
	public List<Module> activeModules(ModuleTypes type) {
		List<Module> modules = new ArrayList<Module>();
		for (Orf orf : orfs) {
			for (Module module : orf.modules())
				if (module.type() == type && module.isActive())
					modules.add(module);
		}
		return modules;
	}

	/**
	 * Get all biosynthetic modules within a cluster.
	 * 
	 * @return all biosynthetic modules within this cluster
	 */
	public List<Module> modules() {
		List<Module> modules = new ArrayList<Module>();
		Sorter.sortOrfs(orfs);
		for (Orf orf : orfs)
			for (Module module : orf.modules())
				modules.add(module);
		return modules;
	}

	/**
	 * Get all orfs containing a biosynthetic module within this cluster.
	 * Note that trans-acyltransferase orfs are not considered to be scaffold orfs. 
	 * 
	 * @return all module-containing orfs
	 */
	public List<Orf> moduleOrfs() {
		List<Orf> moduleOrfs = new ArrayList<Orf>();
		for (Orf orf : orfs)
			if (orf.modules().size() > 0 && orf.hasActiveModule()
					&& !orf.contains(ModuleTypes.TRANS_AT))
				moduleOrfs.add(orf);
		return moduleOrfs;
	}

	/**
	 * Get all orfs containing a domain within this cluster.
	 * 
	 * @return all domain-containing orfs
	 */
	public List<Orf> domainOrfs() {
		List<Orf> domainOrfs = new ArrayList<Orf>();
		for (Orf orf : orfs)
			if (orf.domains().size() > 0)
				domainOrfs.add(orf);
		return domainOrfs;
	}

	/**
	 * Get the first orf in this cluster.
	 * 
	 * @return the cluster's first orf
	 */
	public Orf first() {
		Orf first = null;
		Integer idx = null;
		for (Orf orf : orfs)
			if (idx == null || orf.start() < idx) {
				first = orf;
				idx = orf.start();
			}
		return first;
	}

	/**
	 * Get the last orf in this cluster.
	 * 
	 * @return the cluster's last orf
	 */
	public Orf last() {
		Orf last = null;
		Integer idx = null;
		for (Orf orf : orfs)
			if (idx == null || orf.end() > idx) {
				last = orf;
				idx = orf.end();
			}
		return last;
	}

	/**
	 * Get all domains associated with this cluster.
	 * 
	 * @return all orf domains
	 */
	public List<Domain> domains() {
		List<Domain> domains = new ArrayList<Domain>();
		for (Orf orf : orfs)
			domains.addAll(orf.domains());
		return domains;
	}

	/**
	 * Get all domains of a given type in this cluster.
	 * 
	 * @param type
	 *            the type of domain to return
	 * @return all domains of this type
	 */
	public List<Domain> domains(DomainType type) {
		List<Domain> typeDomains = new ArrayList<Domain>();
		for (Domain domain : domains())
			if (domain.type() == type)
				typeDomains.add(domain);
		return typeDomains;
	}

	/**
	 * Get all domains of a given family in this cluster.
	 * 
	 * @param type
	 *            the type of domain to return
	 * @return all domains of this type
	 */
	public List<Domain> domains(DomainFamilies family) {
		List<Domain> familyDomains = new ArrayList<Domain>();
		for (Domain domain : domains())
			if (domain.family() == family)
				familyDomains.add(domain);
		return familyDomains;
	}

	/**
	 * Get all Tanimoto scores associated with this cluster's scaffold library.
	 * 
	 * @return cluster Tanimoto scores
	 */
	public List<TanimotoScore> scores() {
		return scores;
	}
	
	/**
	 * Get the combinations of sugars that can be envisioned based on this
	 * cluster's sugar biosynthesis genes. Returns a list of lists, where each
	 * list in the list represents a possible combination of sugars.
	 * 
	 * @return all biosynthetically plausible sugar combinations
	 */
	public List<List<Sugar>> sugars() {
		return sugars;
	}

	/**
	 * Set the combinations of sugars that might be produced by this cluster,
	 * based on an analysis of this cluster's sugar biosynthesis genes.
	 * 
	 * @param sugars
	 *            all biosynthetically plausible sugar combinations
	 */
	public void setSugars(List<List<Sugar>> sugars) {
		this.sugars = sugars;
	}
	
	/**
	 * Get information about the combinatorial data evaluated in the
	 * construction of this cluster's scaffold library.
	 * 
	 * @return the package of combinatorial data evaluated
	 */
	public CombinatorialData combinatorialData() {
		return combinatorialData;
	}

}
