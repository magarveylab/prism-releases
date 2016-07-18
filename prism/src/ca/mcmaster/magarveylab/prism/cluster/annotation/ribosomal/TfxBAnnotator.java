package ca.mcmaster.magarveylab.prism.cluster.annotation.ribosomal;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.enums.ModuleTypes;
import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.enums.interfaces.SubstrateType;
import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.SubstrateSet;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;

/**
 * Find substrates for the tetrapeptide substrate of the TfxB (and/or TfxC)
 * enzyme from the trifolitoxin cluster, which catalyzes the formation of the
 * trifolitoxin chromophore. Since the enzymology of this reaction is poorly
 * understood, the presence of TfxB and TfxC is considered sufficient to
 * catalyze the reaction in PRISM.
 * 
 * @author skinnider
 *
 */
public class TfxBAnnotator implements Annotator {

	public DomainType[] domains() {
		return new DomainType[] { RibosomalDomains.TfxB };
	}

	public List<SubstrateSet> findSubstrates(Domain domain,
			List<Module> permutation, Cluster cluster)
			throws InvalidSmilesException, IOException {
		List<SubstrateSet> substrates = new ArrayList<SubstrateSet>();

		if (permutation.size() == 0) 
			return substrates;
		if (!cluster.contains(RibosomalDomains.TfxC))
			return substrates;

		for (int i = 0; i < permutation.size() - 3; i++) {
			Module m1 = permutation.get(i);
			Module m2 = permutation.get(i + 1);
			Module m3 = permutation.get(i + 2);
			Module m4 = permutation.get(i + 3);

			if (m1.type() != ModuleTypes.RIBOSOMAL || m1.scaffold() == null
					|| m1.scaffold().topSubstrate() == null)
				continue;
			if (m2.type() != ModuleTypes.RIBOSOMAL || m2.scaffold() == null
					|| m2.scaffold().topSubstrate() == null)
				continue;
			if (m3.type() != ModuleTypes.RIBOSOMAL || m3.scaffold() == null
					|| m3.scaffold().topSubstrate() == null)
				continue;
			if (m4.type() != ModuleTypes.RIBOSOMAL || m4.scaffold() == null
					|| m4.scaffold().topSubstrate() == null)
				continue;

			SubstrateType aa2 = m2.scaffold().topSubstrate().type();
			SubstrateType aa3 = m3.scaffold().topSubstrate().type();
			SubstrateType aa4 = m4.scaffold().topSubstrate().type();
			if (aa2 == ProteinogenicAminoAcids.GLUTAMINE
					&& aa3 == ProteinogenicAminoAcids.GLYCINE
					&& aa4 == ProteinogenicAminoAcids.CYSTEINE) {
				SubstrateSet substrate = new SubstrateSet();
				substrate.add(m1);
				substrate.add(m2);
				substrate.add(m3);
				substrate.add(m4);
				substrates.add(substrate);
			}
		}

		return substrates;
	}

}
