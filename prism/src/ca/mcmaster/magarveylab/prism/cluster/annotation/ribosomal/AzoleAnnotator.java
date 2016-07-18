package ca.mcmaster.magarveylab.prism.cluster.annotation.ribosomal;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.clusters.RibosomalClusterTypes;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;
import ca.mcmaster.magarveylab.prism.cluster.analysis.RibosomalClusterAnalyzer;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.cluster.reactions.RibosomalUtil;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;

/**
 * Get potential points of reaction for the cyanobactin heterocyclase PatD,
 * which forms oxazolines and thiazolines (but never bis-thiazoline/oxazoline
 * systems!), the linear azole-containing peptide dehydrogenase B and
 * cyclodehydratases C/D, the thiopeptide cyclodehydratase and oxidase LazE and
 * LazF, the bottromycin cyclodehydratase BotCD, the N- and C-terminal domains
 * of the YM-216391 family fused cyclodehydratase/oxidase YmBC, and the unique
 * goadsporin enymes GodG and GodF.
 * 
 * @author skinnider
 *
 */
public class AzoleAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.PatD, RibosomalDomains.McbB,
				RibosomalDomains.McbC, RibosomalDomains.McbD,
				RibosomalDomains.LazE, RibosomalDomains.LazF,
				RibosomalDomains.YmBC_a, RibosomalDomains.YmBC_b,
				RibosomalDomains.BotCD, RibosomalDomains.GodD,
				RibosomalDomains.GodE };
	}

	public List<SubstrateSet> findSubstrates(Domain domain,
			List<Module> permutation, Cluster cluster)
					throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		List<Module> serine = RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.SERINE, permutation);
		List<Module> threonine = RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.THREONINE, permutation);
		List<Module> cysteine = RibosomalClusterAnalyzer.getModules(
				ProteinogenicAminoAcids.CYSTEINE, permutation);
		
		// remove N-terminal residues
		Iterator<Module> serItr = serine.iterator();
		while (serItr.hasNext()) 
			if (permutation.indexOf(serItr.next()) == 0)
				serItr.remove();
		Iterator<Module> thrItr = threonine.iterator();
		while (thrItr.hasNext()) 
			if (permutation.indexOf(thrItr.next()) == 0)
				thrItr.remove();
		Iterator<Module> cysItr = cysteine.iterator();
		while (cysItr.hasNext()) 
			if (permutation.indexOf(cysItr.next()) == 0)
				cysItr.remove();

		List<Module> input = new ArrayList<Module>();
		if (domain.type() != RibosomalDomains.McbC) { 
			// McbC acts on Cys only
			input.addAll(serine);
			input.addAll(threonine);
		}
		if (domain.type() != RibosomalDomains.McbD) { 
			// McbD acts on Ser/Thr only
			input.addAll(cysteine);
		}
		
		System.out.println("Got " + input.size() + " modules for domain " + domain.type());

		// get subsets of modules for combinatorialization
		List<List<Module>> subsets = new ArrayList<List<Module>>();
		if (domain.type() == RibosomalDomains.YmBC_a) {
			// YmBC_a acts at all
			List<Module> subset = new ArrayList<Module>();
			subset.addAll(serine);
			subset.addAll(threonine);
			subset.addAll(cysteine);
			subsets.add(subset);
		} else if (domain.type() == RibosomalDomains.McbC) {
			subsets = RibosomalUtil.getSubsets(input, input.size());
		} else {
			if (cluster.types().contains(
					RibosomalClusterTypes.LINEAR_AZOLE_CONTAINING_PEPTIDE)) {
				int min = input.size() - 4 > 0 ? input.size() - 4 : 0;
				subsets = RibosomalUtil.getSubsetsWithMinSize(input, min);
			} else {
				subsets = RibosomalUtil.getSubsets(input);
			}
		}

		Iterator<List<Module>> itr = subsets.iterator();
		subsetLoop: while (itr.hasNext()) {
			List<Module> next = itr.next();

			// check there are no bis-thiazoles/oxazoles for PatD
			if (domain.type() == RibosomalDomains.PatD)
				for (int i = 0; i < next.size() - 1; i++) {
					Module m1 = next.get(i);
					Module m2 = next.get(i + 1);
					if (Math.abs(permutation.indexOf(m1)
							- permutation.indexOf(m2)) == 1) {
						itr.remove();
						continue subsetLoop;
					}
				}

			// create SubstrateSet
			SubstrateSet substrate = new SubstrateSet();
			for (int i = 0; i < next.size(); i++) {
				Module module = next.get(i);
				Module last = permutation.get(permutation.indexOf(module) - 1);
				substrate.add(module);
				substrate.add(last);
			}
			substrates.add(substrate);
		}

		return substrates;
	}

}
