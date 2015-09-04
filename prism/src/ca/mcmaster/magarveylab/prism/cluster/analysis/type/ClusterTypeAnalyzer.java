package ca.mcmaster.magarveylab.prism.cluster.analysis.type;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ca.mcmaster.magarveylab.enums.clusters.*;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.util.Strings;

public class ClusterTypeAnalyzer {

	public static List<ClusterType> getAllClusterTypes() {
		List<ClusterType> all = new ArrayList<ClusterType>();
		all.addAll(Arrays.asList(AminoglycosideClusterTypes.values()));
		all.addAll(Arrays.asList(ThiotemplatedClusterTypes.values()));
		all.addAll(Arrays.asList(TypeIIPolyketideClusterTypes.values()));
		all.addAll(Arrays.asList(NucleosideClusterTypes.values()));
		all.addAll(Arrays.asList(RibosomalClusterTypes.values()));
		all.addAll(Arrays.asList(BetaLactamClusterTypes.values()));
		return all;
	}
	
	public static String getTypeString(Cluster cluster) {
		List<ClusterType> types = cluster.types();
		StringBuffer sb = new StringBuffer();
		
		if (types.size() == 0) 
			return "Could not determine type";
			
		ClusterType first = types.get(0);
		String capitalized = Strings.capitalizeFirstLetter(first.fullName());
		sb.append(capitalized);
		
		if (types.size() > 1)
			sb.append("/");
		
		for (int i = 1; i < types.size() - 1; i++) {
			ClusterType t = types.get(i);
			sb.append(t.fullName() + "/");
		}
		
		if (types.size() > 1) {
			ClusterType last = types.get(types.size() - 1);
			sb.append(last.fullName());
		}
		
		return sb.toString();
	}
	
}
