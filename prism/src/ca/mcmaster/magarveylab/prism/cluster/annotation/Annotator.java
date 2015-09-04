package ca.mcmaster.magarveylab.prism.cluster.annotation;

import java.io.IOException;
import java.util.List;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;

/**
 * Annotates potential sites for biosynthetic reactions.
 * @author skinnider
 *
 */
public interface Annotator {
	
	/**
	 * Get the types of domains which use this annotator to determine potential substrates.
	 * @return	all domain types which use this annotator 
	 */
	public DomainType[] domains();
 
	/**
	 * Find all potential reaction substrates for a particular domain within a given module permutation. 
	 * @param domain		the domain to analyze
	 * @param permutation	the module permutation for which substrates should be found
	 * @param cluster		the cluster in which this domain and module permutation is found 
	 * @return				all potential substrate sets 
	 * @throws InvalidSmilesException
	 * @throws IOException
	 * @throws CDKException 
	 */
	public List<SubstrateSet> findSubstrates(Domain domain, List<Module> permutation, Cluster cluster) 
			throws IOException, CDKException;
	
}
