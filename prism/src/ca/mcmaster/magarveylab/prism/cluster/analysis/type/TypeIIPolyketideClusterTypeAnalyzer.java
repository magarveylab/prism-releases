package ca.mcmaster.magarveylab.prism.cluster.analysis.type;

import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.enums.clusters.ClusterType;
import ca.mcmaster.magarveylab.enums.clusters.TypeIIPolyketideClusterTypes;
import ca.mcmaster.magarveylab.enums.domains.TypeIIPolyketideDomains;
import ca.mcmaster.magarveylab.prism.cluster.analysis.TypeIIPolyketideAnalyzer;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;

public class TypeIIPolyketideClusterTypeAnalyzer {
	
	/**
	 * Get the subtype(s) of this type II polyketide cluster.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return type II polyketide-specific cluster subtype(s)
	 */
	public static List<ClusterType> getTypes(Cluster cluster) {
		 List<ClusterType> types = new ArrayList<ClusterType>();
		if (TypeIIPolyketideAnalyzer.isTypeIIPolyketideCluster(cluster))
			for (Domain domain : cluster.domains(TypeIIPolyketideDomains.CLF)) {
				if (domain.blastResults().size() == 0) {
					continue;
				}
				String name = domain.blastResults().get(0).subject();
				if (name.contains("enterocin")) {
					types.add(TypeIIPolyketideClusterTypes.ENTEROCIN);
				} else if (name.contains("resistomycin")) {
					types.add(TypeIIPolyketideClusterTypes.RESISTOMYCIN);
				} else if (name.contains("oxytetracycline")
						|| name.contains("chlortetracycline")
						|| name.contains("dactylocycline")
						|| name.contains("SF2575")
						|| name.contains("chelocardin")
						) {
					types.add(TypeIIPolyketideClusterTypes.TETRACYCLINE);
				} else if (name.contains("alnumycin")
						|| name.contains("actinorhodin")
						|| name.contains("frenolicin")
						|| name.contains("griseusin")
						|| name.contains("granaticin")
						|| name.contains("medermycin")
						|| name.contains("naphthocyclinone")
						) {
					types.add(TypeIIPolyketideClusterTypes.BENZOISOCHROMANEQUINONE);
				} else if (name.contains("daunorubicin")
						|| name.contains("doxorubicin")
						|| name.contains("aclacinomycin")
						|| name.contains("aranciamycin")
						|| name.contains("chartreusin")
						|| name.contains("cosmomycin")
						|| name.contains("nogalamycin")
						|| name.contains("steffimycin")
						) {
					types.add(TypeIIPolyketideClusterTypes.ANTHRACYCLINE);
				} else if (name.contains("auricin")
						|| name.contains("azicemicin")
						|| name.contains("BE-7585A")
						|| name.contains("chattamycin")
						|| name.contains("chrysomycin")
						|| name.contains("gilvocarcin")
						|| name.contains("grincamycin")
						|| name.contains("hatomarubigin")
						|| name.contains("jadomycin")
						|| name.contains("kinamycin")
						|| name.contains("landomycin")
						|| name.contains("lomaiviticin")
						|| name.contains("oviedomycin")
						|| name.contains("PD-116740")
						|| name.contains("R1128")
						|| name.contains("radivomycin")
						|| name.contains("saquayamycin")
						|| name.contains("Sch-47554")
						|| name.contains("simocyclinone")
						|| name.contains("urdamycin")
						) {
					types.add(TypeIIPolyketideClusterTypes.ANGUCYCLINE);
				} else if (name.contains("chromomycin")
						|| name.contains("mithramycin")
						|| name.contains("X26")
						) {
					types.add(TypeIIPolyketideClusterTypes.AUREOLIC_ACID);
				} else if (name.contains("tetracenomycin")
						|| name.contains("elloramycin")
						|| name.contains("tetarimycin")
						) {
					types.add(TypeIIPolyketideClusterTypes.TETRACENOMYCIN);
				} else if (name.contains("A-74528")
						|| name.contains("AZ154")
						|| name.contains("benastatin")
						|| name.contains("FD-594")
						|| name.contains("fredericamycin")
						|| name.contains("griseorhodin")
						|| name.contains("lysolipin")
						|| name.contains("pradimicin")
						|| name.contains("rubromycin")
						|| name.contains("TLN-05220")
						|| name.contains("WhiE")
						|| name.contains("xantholipin")
						) {
					types.add(TypeIIPolyketideClusterTypes.PENTANGULAR_POLYPHENOL);
				} else if (name.contains("hedamycin")) { 
					types.add(TypeIIPolyketideClusterTypes.PLURAMYCIN);
				} else {
					types.add(TypeIIPolyketideClusterTypes.OTHER);
				}
			}
		return types;
	}

}
