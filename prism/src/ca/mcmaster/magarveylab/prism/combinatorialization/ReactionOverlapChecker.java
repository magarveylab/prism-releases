package ca.mcmaster.magarveylab.prism.combinatorialization;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.enums.domains.TailoringDomains;
import ca.mcmaster.magarveylab.enums.interfaces.SubstrateType;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;
import ca.mcmaster.magarveylab.prism.data.Module;
import ca.mcmaster.magarveylab.prism.data.reactions.ReactionPlan;
import ca.mcmaster.magarveylab.prism.util.exception.TailoringSubstrateException;

/**
 * Some tailoring reactions are highly likely to overlap with other tailoring
 * reactions. This can create an excessive memory demand. This method helps
 * remove overlapping reaction plans with a high performance impact.
 * 
 * @author skinnider
 *
 */
public class ReactionOverlapChecker {

	private List<ReactionPlan> combination; 
	private List<Module> permutation;

	public ReactionOverlapChecker(List<ReactionPlan> combination,
			List<Module> permutation) {
		this.combination = combination;
		this.permutation = permutation;
	}

	public boolean check() {
		if (!checkAzoles())
			return false;
		if (!checkAzoleAzolines())
			return false;
		if (!checkLanC())
			return false;
		if (!checkCannotOverlap())
			return false;
		if (!checkMustOverlapOnce())
			return false;
		if (!checkMustOverlapTwice())
			return false;
		if (!checkThiopeptideEsterification())
			return false;
		if (!checkPyridineSubstituents())
			return false;
		return true;
	}

	/**
	 * Check reaction plan combinations for overlaps between domains which must
	 * overlap at at least one site.
	 * 
	 * @return true if the combination is allowable; false if it is not
	 */
	public boolean checkMustOverlapOnce() {
		boolean flag = true;

		List<DomainType[]> mustOverlap = new ArrayList<DomainType[]>();
		// LanB-LazB and NocQ MUST overlap (Dhb methoxylation) 
		mustOverlap.add(new DomainType[] { RibosomalDomains.LazB, RibosomalDomains.NocQ });
		// if pyridine is hydroxylated, must be at a pyridine site
		mustOverlap.add(new DomainType[] { RibosomalDomains.LazC, RibosomalDomains.NosC });
		// CltM acts at Dha
		mustOverlap.add(new DomainType[] { RibosomalDomains.LazB, RibosomalDomains.CltM });
		// LanB must act before LanD (LanD acts at Dha)  
		mustOverlap.add(new DomainType[] { RibosomalDomains.LanB, RibosomalDomains.LanD });

		for (DomainType[] overlap : mustOverlap) {
			boolean[] usedDomain = new boolean[overlap.length];
			boolean[][] used = new boolean[permutation.size()][overlap.length];
			for (ReactionPlan plan : combination) {
				DomainType type = plan.type();
				int typeIdx = Arrays.asList(overlap).indexOf(type);
				if (typeIdx == -1) {
					continue;
				} else {
					usedDomain[typeIdx] = true;
				}

				List<Module> modules = plan.modules().getAllModules();
				for (Module module : modules) {
					int moduleIdx = permutation.indexOf(module);
					used[moduleIdx][typeIdx] = true;
				}
			}

			// if all domains are used...
			if (!containsFalse(usedDomain)) {
				// ... and this reaction plan combination doesn't have 1 or
				// more sites of overlap, reject it
				int overlaps = 0;
				for (boolean[] u : used) 
					if (!containsFalse(u)) 
						overlaps++;
				if (overlaps < 1) {
					return false;
				}
			}

		}
		return flag;
	}

	/**
	 * Check reaction plan combinations for overlaps between domains which must
	 * overlap at at least two sites, such as the pyridine-forming cycloaddition
	 * enzyme LazC and the thiopeptide dehydratase LazB, or LazC and the
	 * pyridine hydroxylase NosC.
	 * 
	 * @return true if the combination is allowable; false if it is not
	 */
	public boolean checkMustOverlapTwice() {
		boolean flag = true;

		List<DomainType[]> mustOverlap = new ArrayList<DomainType[]>();
		// to form pyridine, need 2 dehydrated serines
		mustOverlap.add(new DomainType[] { RibosomalDomains.LazC, RibosomalDomains.LazB });
		mustOverlap.add(new DomainType[] { RibosomalDomains.LazC_b, RibosomalDomains.LazB });
		mustOverlap.add(new DomainType[] { RibosomalDomains.LazC, RibosomalDomains.NosC });

		for (DomainType[] overlap : mustOverlap) {
			boolean[] usedDomain = new boolean[overlap.length];
			boolean[][] used = new boolean[permutation.size()][overlap.length];
			for (ReactionPlan plan : combination) {
				DomainType type = plan.type();
				int typeIdx = Arrays.asList(overlap).indexOf(type);
				if (typeIdx == -1) {
					continue;
				} else {
					usedDomain[typeIdx] = true;
				}

				List<Module> modules = plan.modules().getAllModules();
				for (Module module : modules) {
					int moduleIdx = permutation.indexOf(module);
					used[moduleIdx][typeIdx] = true;
				}
			}

			// if all domains are used...
			if (!containsFalse(usedDomain)) {
				// ... and this reaction plan combination doesn't have 2 or
				// more sites of overlap, reject it
				int overlaps = 0;
				for (boolean[] u : used) 
					if (!containsFalse(u))
						overlaps++;
				if (overlaps < 2) 
					return false;
			}
		}
		return flag;
	}

	/**
	 * Check reaction plan combinations for overlaps between domains which
	 * cannot overlap at any sites.
	 * 
	 * @return true if the combination is allowable; false if it is not
	 */
	public boolean checkCannotOverlap() {
		boolean flag = true;

		List<DomainType[]> cannotOverlap = new ArrayList<DomainType[]>();
		cannotOverlap.add(new DomainType[] { RibosomalDomains.TsrI, TailoringDomains.GLYCOSYLTRANSFERASE });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.NosI, TailoringDomains.GLYCOSYLTRANSFERASE });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.NocQ, TailoringDomains.GLYCOSYLTRANSFERASE });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.LazB, TailoringDomains.GLYCOSYLTRANSFERASE });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.LanB, TailoringDomains.GLYCOSYLTRANSFERASE });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.LanC, TailoringDomains.GLYCOSYLTRANSFERASE });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.LanM, TailoringDomains.GLYCOSYLTRANSFERASE });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.LanKC, TailoringDomains.GLYCOSYLTRANSFERASE });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.SunA, RibosomalDomains.SunS });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.BdbB, RibosomalDomains.SunS });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.TpdI, RibosomalDomains.TpdL });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.TpaJ, RibosomalDomains.LazB });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.TpaJ, RibosomalDomains.LazE });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.LanC, RibosomalDomains.LanD });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.NosA, RibosomalDomains.LazB });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.NosA, RibosomalDomains.LazC });

		for (DomainType[] overlap : cannotOverlap) {
			boolean[] usedDomain = new boolean[overlap.length];
			boolean[][] used = new boolean[permutation.size()][overlap.length];
			for (ReactionPlan plan : combination) {
				DomainType type = plan.type();
				int typeIdx = Arrays.asList(overlap).indexOf(type);
				if (typeIdx == -1) {
					continue;
				} else {
					usedDomain[typeIdx] = true;
				}

				List<Module> modules = plan.modules().getAllModules();
				for (Module module : modules) {
					if (module == null) // LanKC annotator 
						continue;
					int moduleIdx = permutation.indexOf(module);
					used[moduleIdx][typeIdx] = true;
				}
			}

			// if all domains are used...
			if (!containsFalse(usedDomain)) {
				// ... and this reaction plan combination doesn't have 1 or
				// more sites of overlap, reject it
				int overlaps = 0;
				for (boolean[] u : used) 
					if (!containsFalse(u)) {
						overlaps++;
					}
				if (overlaps > 0) 
					return false;
			}

		}
		return flag;
	}

	/**
	 * Check reaction plan combinations for overlaps between thiopeptide
	 * esterification enzymes and thiopeptide dehydratases.
	 * 
	 * @return true if the combination is allowable; false if it is not
	 */
	public boolean checkThiopeptideEsterification() {
		boolean flag = true;

		List<DomainType[]> cannotOverlap = new ArrayList<DomainType[]>();
		// NosI and LazE cannot overlap 
		cannotOverlap.add(new DomainType[] { RibosomalDomains.LazE, RibosomalDomains.NosI });
		// NosI and LazB cannot overlap 
		cannotOverlap.add(new DomainType[] { RibosomalDomains.LazB, RibosomalDomains.NosI });
		// TsrI and LazB cannot overlap 
		cannotOverlap.add(new DomainType[] { RibosomalDomains.LazB, RibosomalDomains.TsrI });

		for (DomainType[] overlap : cannotOverlap) {
			boolean[] usedDomain = new boolean[overlap.length];
			boolean[][] used = new boolean[permutation.size()][overlap.length];
			for (ReactionPlan plan : combination) {
				DomainType type = plan.type();
				int typeIdx = Arrays.asList(overlap).indexOf(type);
				if (typeIdx == -1) {
					continue;
				} else {
					usedDomain[typeIdx] = true;
				}

				List<Module> modules = plan.modules().getAllModules();
				if (type == RibosomalDomains.LazB) {
					for (Module module : modules) {
						int moduleIdx = permutation.indexOf(module);
						used[moduleIdx][typeIdx] = true;
					}
				} else if (type == RibosomalDomains.LazE) {
					for (int i = 0; i < modules.size() - 1; i += 2) {
						try {
							Module module = plan.get(i);
							int moduleIdx = permutation.indexOf(module);
							used[moduleIdx][typeIdx] = true;
						} catch (TailoringSubstrateException e) {
							continue;
						}
					}
				} else if (type == RibosomalDomains.NosI) {
					try {
						Module module = plan.get(0);
						int moduleIdx = permutation.indexOf(module);
						used[moduleIdx][typeIdx] = true;
					} catch (TailoringSubstrateException e) {
						continue;
					}
				} else if (type == RibosomalDomains.TsrI) {
					try {
						Module module = plan.get(1);
						int moduleIdx = permutation.indexOf(module);
						used[moduleIdx][typeIdx] = true;
					} catch (TailoringSubstrateException e) {
						continue;
					}
				}
			}

			// if all domains are used...
			if (!containsFalse(usedDomain)) {
				// ... and this reaction plan combination has a site of
				// overlap
				int overlaps = 0;
				for (boolean[] u : used) 
					if (!containsFalse(u))
						overlaps++;
				if (overlaps > 0)
					return false;
			}
		}
		return flag;
	}

	/**
	 * Check that LanB acts before LanC at all potential Dha/Dhb residues.
	 * 
	 * @return true if the combination is allowable; false if it is not
	 */
	public boolean checkLanC() {
		boolean flag = true;

		List<DomainType[]> mustOverlap = new ArrayList<DomainType[]>();
		mustOverlap.add(new DomainType[] { RibosomalDomains.LanB, RibosomalDomains.LanC });

		for (DomainType[] overlap : mustOverlap) {
			boolean[] usedDomain = new boolean[overlap.length];
			boolean[][] used = new boolean[permutation.size()][overlap.length];
			for (ReactionPlan plan : combination) {
				DomainType type = plan.type();
				int typeIdx = Arrays.asList(overlap).indexOf(type);
				if (typeIdx == -1) {
					continue;
				} else {
					usedDomain[typeIdx] = true;
				}

				List<Module> modules = plan.modules().getAllModules();
				if (type == RibosomalDomains.LanC) {
					for (int i = 0; i < modules.size() - 1; i += 2)
						try {
							Module module = plan.get(i);
							int moduleIdx = permutation.indexOf(module);
							used[moduleIdx][typeIdx] = true;
						} catch (TailoringSubstrateException e) {
							continue;
						}
				} else {
					for (Module module : modules) {
						int moduleIdx = permutation.indexOf(module);
						used[moduleIdx][typeIdx] = true;
					}
				}
			}

			if (!containsFalse(usedDomain)) 
				for (boolean[] u : used)
					if (u[0] == false && u[1] == true) { // LanB and LanC don't overlap 
						return false;
					}
		}

		return flag;
	}

	/**
	 * Check that reaction plans contain required overlaps, and do not contain
	 * prohibited overlaps, with azol(in)e-forming domains. Azol(in)es require a
	 * unique method because reaction plans contain two domains, and only the
	 * first typically requires a specific overlap checker.
	 * 
	 * @return true if the combination is allowable; false if it is not
	 */
	public boolean checkAzoles() {
		boolean flag = true;

		List<DomainType[]> mustOverlap = new ArrayList<DomainType[]>();
		// azole transferases must overlap with azole+azoline 
		mustOverlap.add(new DomainType[] { RibosomalDomains.LazE, RibosomalDomains.TpdI });
		mustOverlap.add(new DomainType[] { RibosomalDomains.LazE, RibosomalDomains.TpdM });
		mustOverlap.add(new DomainType[] { RibosomalDomains.LazE, RibosomalDomains.TpdL });
		mustOverlap.add(new DomainType[] { RibosomalDomains.LazF, RibosomalDomains.TpdI });
		mustOverlap.add(new DomainType[] { RibosomalDomains.LazF, RibosomalDomains.TpdM });
		mustOverlap.add(new DomainType[] { RibosomalDomains.LazF, RibosomalDomains.TpdL });

		List<DomainType[]> cannotOverlap = new ArrayList<DomainType[]>();
		// azoles or azolines can't overlap with LazB (LanB) or LazC 
		cannotOverlap.add(new DomainType[] { RibosomalDomains.LazB, RibosomalDomains.LazE });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.LazB, RibosomalDomains.LazF });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.LanB, RibosomalDomains.LazE });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.LanB, RibosomalDomains.LazF });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.LazC, RibosomalDomains.LazE });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.LazC, RibosomalDomains.LazF });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.LazC_b, RibosomalDomains.LazE });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.LazC_b, RibosomalDomains.LazF });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.GodF, RibosomalDomains.GodD });
		cannotOverlap.add(new DomainType[] { RibosomalDomains.GodF, RibosomalDomains.GodE });
		// oxazoles or oxazolines cannot overlap with GTrs 
		cannotOverlap.add(new DomainType[] { TailoringDomains.GLYCOSYLTRANSFERASE, RibosomalDomains.LazE });

		// now merge the two
		List<DomainType[]> o = new ArrayList<DomainType[]>();
		o.addAll(mustOverlap);
		o.addAll(cannotOverlap);

		for (DomainType[] overlap : o) {
			boolean[] usedDomain = new boolean[overlap.length];
			boolean[][] used = new boolean[permutation.size()][overlap.length];
			for (ReactionPlan plan : combination) {
				DomainType type = plan.type();
				int typeIdx = Arrays.asList(overlap).indexOf(type);
				if (typeIdx == -1) {
					continue;
				} else {
					usedDomain[typeIdx] = true;
				}

				List<Module> modules = plan.modules().getAllModules();
				if (type == RibosomalDomains.LazE
						|| type == RibosomalDomains.LazF
						|| type == RibosomalDomains.GodD
						|| type == RibosomalDomains.GodE) {
					for (int i = 0; i < modules.size() - 1; i += 2)
						try {
							Module module = plan.get(i);
							int moduleIdx = permutation.indexOf(module);
							used[moduleIdx][typeIdx] = true;
						} catch (TailoringSubstrateException e) {
							continue;
						}
				} else if (type == RibosomalDomains.LazC
						|| type == RibosomalDomains.LazC_b) {
					for (int i = 0; i < modules.size() - 1; i ++)
						try {
							Module module = plan.get(i);
							int moduleIdx = permutation.indexOf(module);
							used[moduleIdx][typeIdx] = true;
						} catch (TailoringSubstrateException e) {
							continue;
						}
				} else {
					for (Module module : modules) {
						int moduleIdx = permutation.indexOf(module);
						used[moduleIdx][typeIdx] = true;
					}
				}
			}

			if (!containsFalse(usedDomain)) {
				int overlaps = 0;
				for (boolean[] u : used) 
					if (!containsFalse(u))
						overlaps++;
				if (mustOverlap.contains(overlap)) {
					// if the domains must overlap, reject if they don't
					if (overlaps < 1) 
						return false;
				} else if (cannotOverlap.contains(overlap)) {
					// if the domains cannot overlap, reject if they do
					if (overlaps > 0) 
						return false;
				}
			}
		}

		return flag;
	}

	/**
	 * Check that LazF never acts at a site that LazE does not also act on.
	 * 
	 * @return true if the combination is allowable; false if it is not
	 */
	public boolean checkAzoleAzolines() {
		boolean flag = true; 

		// LazE and LazF must overlap everywhere 
		List<DomainType[]> mustOverlap = new ArrayList<DomainType[]>();
		mustOverlap.add(new DomainType[] { RibosomalDomains.LazE, RibosomalDomains.LazF });
		mustOverlap.add(new DomainType[] { RibosomalDomains.GodD, RibosomalDomains.GodE });
		mustOverlap.add(new DomainType[] { RibosomalDomains.YmBC_a, RibosomalDomains.YmBC_b });

		for (DomainType[] overlap : mustOverlap) {
			boolean[] usedDomain = new boolean[overlap.length];
			boolean[][] used = new boolean[permutation.size()][overlap.length];
			for (ReactionPlan plan : combination) {
				DomainType type = plan.type();
				int typeIdx = Arrays.asList(overlap).indexOf(type);
				if (typeIdx == -1) {
					continue;
				} else {
					usedDomain[typeIdx] = true;
				}

				List<Module> modules = plan.modules().getAllModules();
				for (int i = 0; i < modules.size() - 1; i += 2)
					try {
						Module module = plan.get(i);
						int moduleIdx = permutation.indexOf(module);
						used[moduleIdx][typeIdx] = true;
					} catch (TailoringSubstrateException e) {
						continue;
					}
			}

			if (!containsFalse(usedDomain)) {
				// need to see LazE at all LazF sites
				for (boolean[] u : used) 
					if (u[0] == false && u[1] == true) 
						return false;
			}
		}

		return flag; 
	}

	/**
	 * Substituents at the 2 and 3 positions of the pyridine ring are always
	 * heterocyclized and fully oxidized.
	 * 
	 * @return true if the combination is allowable; false if it is not
	 */
	public boolean checkPyridineSubstituents() {
		// both LazE and LazF must overlap with the residues at pyridine positions 2 and 3, if possible
		List<DomainType[]> mustOverlap = new ArrayList<DomainType[]>();
		mustOverlap.add(new DomainType[] { RibosomalDomains.LazC, RibosomalDomains.LazE });
		mustOverlap.add(new DomainType[] { RibosomalDomains.LazC, RibosomalDomains.LazF });
		mustOverlap.add(new DomainType[] { RibosomalDomains.LazC_b, RibosomalDomains.LazE });
		mustOverlap.add(new DomainType[] { RibosomalDomains.LazC_b, RibosomalDomains.LazF });

		for (DomainType[] overlap : mustOverlap) {
			// the cluster must contain all domains to execute this check
			ReactionPlan plan1 = null;
			ReactionPlan plan2 = null;
			for (ReactionPlan plan : combination) {
				if (plan.type() == overlap[0])
					plan1 = plan;
				if (plan.type() == overlap[1])
					plan2 = plan;
			}
			if (plan1 == null || plan2 == null)
				continue;

			// get LazC sites
			int site1, site2;
			try {
				int dha1 = permutation.indexOf(plan1.get(0));
				site1 = dha1 + 1;
				site2 = permutation.indexOf(plan1.get(2));
			} catch (TailoringSubstrateException e) {
				continue;
			}

			// now, test if these residues must be heterocyclized/oxidized
			Module m1 = permutation.get(site1);
			Module m2 = permutation.get(site2);
			if (m1.scaffold() == null || m1.scaffold().topSubstrate() == null)
				continue;
			if (m2.scaffold() == null || m2.scaffold().topSubstrate() == null)
				continue;
			SubstrateType s1 = m1.scaffold().topSubstrate().type();
			SubstrateType s2 = m2.scaffold().topSubstrate().type();

			if (s1 != ProteinogenicAminoAcids.CYSTEINE
					&& s1 != ProteinogenicAminoAcids.SERINE
					&& s1 != ProteinogenicAminoAcids.THREONINE)
				continue;
			if (s2 != ProteinogenicAminoAcids.CYSTEINE
					&& s2 != ProteinogenicAminoAcids.SERINE
					&& s2 != ProteinogenicAminoAcids.THREONINE)
				continue;

			// finally, check for overlap between C/S/T residue and LazC 
			List<Module> modules = plan2.modules().getAllModules();
			if (overlap[1] == RibosomalDomains.LazE) {
				boolean containsM1 = false, containsM2 = false;
				for (int i = 0; i < modules.size() - 1; i += 2) {
					Module module = modules.get(i);
					if (module == m1)
						containsM1 = true;
					if (module == m2)
						containsM2 = true;
				}
				if (!containsM1 || !containsM2)
					return false;
			} else {
				if (!modules.contains(m1) || !modules.contains(m2))
					return false;
			}
		}

		return true;
	}

	private boolean containsFalse(boolean[] array) {
		boolean containsFalse = false;
		for (boolean b : array)
			if (b == false)
				containsFalse = true;
		return containsFalse;
	}

}
