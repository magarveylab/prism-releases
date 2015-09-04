package ca.mcmaster.magarveylab.prism.cluster.analysis;

import java.util.List;

import ca.mcmaster.magarveylab.enums.OrfTypes;
import ca.mcmaster.magarveylab.prism.cluster.module.ModuleFinder;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.Orf;

/**
 * Detects insertion modules in trans-acyltransferase polyketide synthases which
 * are split across multiple open reading frames by creating 'artificial' open
 * reading frames, which are used to generate the polyketide scaffold.
 * 
 * @author skinnider
 *
 */
public class ArtificialOrfCreator {

	int counter;
	List<Orf> orfs;

	/**
	 * Instantiate a new artificial open reading frame creator. 
	 * @param orfs	copy of the open reading frames in the trans-AT cluster 
	 */
	public ArtificialOrfCreator(List<Orf> orfs) {
		this.counter = 0;
		this.orfs = orfs;
	}

	/**
	 * Construct artificial open reading frames.
	 */
	public void initiate() {
		Orf artificial = new Orf("Artificial orf", "");
		Orf orf = orfs.get(counter);

		if (orf.hasActiveModule()) {
			List<Domain> domains = orf.domains();
			Module m1 = orf.getLastModule();
			Domain last = m1.last();
			int lastModuleDomainIndex = domains.indexOf(last);
			if (domains.size() - lastModuleDomainIndex > 1) {
				System.out.println("[ArtificialOrfCreator] Expanding from initial orf "
								+ orf.name());
				for (int i = lastModuleDomainIndex + 1; i < domains.size(); i++)
					artificial.add(domains.get(i));
			} else {
				System.out.println("[ArtificialOrfCreator] Skipping possible initial orf "
								+ orf.name());
				counter++;
				if (!abort())
					initiate();
			}
		} else {
			System.out.println("[ArtificialOrfCreator] Expanding from initial orf without domains "
							+ orf.name());
			artificial.add(orf.domains());
		}

		if (!abort()) {
			counter++;
			expand(artificial);
		}
	}

	/**
	 * Recursively expand construction of an artificial open reading frame.
	 * @param artificial	the artificial open reading frame to expand 
	 */
	public void expand(Orf artificial) {
		Orf orf = orfs.get(counter);

		if (orf.hasActiveModule()) {
			System.out.println("[ArtificialOrfCreator] Terminating with module-containing orf "
							+ orf.name());

			List<Domain> domains = orf.domains();
			Module m2 = orf.modules().get(0);
			Domain first = m2.first();
			if (domains.indexOf(first) > 0)
				for (int i = 0; i < domains.indexOf(first); i++)
					artificial.add(domains.get(i));

			StringBuffer sb = new StringBuffer();
			for (Domain domain : artificial.domains())
				sb.append(domain.type() + ", ");
			System.out.println("\tDomains: " + sb.toString());

			ModuleFinder.detectModules(artificial);
			if (artificial.modules().size() > 0) {
				orf.setType(OrfTypes.PKS);
				orfs.get(counter-1).setType(OrfTypes.PKS);
			}
			orfs.add(counter, artificial);
			counter++;

			if (!abort())
				initiate();
		} else {
			System.out.println("[ArtificialOrfCreator] Expanding with module-containing orf "
							+ orf.name());

			artificial.add(orf.domains());
			counter++;
			if (!abort()) {
				expand(artificial);
			} else {
				ModuleFinder.detectModules(artificial);
				orfs.add(counter - 1, artificial);
				return;
			}
		}
	}

	/**
	 * Determine whether or not to stop the recursive process of artificial orf creation.
	 * @return	true if the current orf is the last 
	 */
	public boolean abort() {
		return (orfs.size() - counter <= 1);
	}

}
