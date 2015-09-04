package ca.mcmaster.magarveylab.prism.cluster.analysis;

import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;
import ca.mcmaster.magarveylab.prism.data.Cluster;
import ca.mcmaster.magarveylab.prism.data.Domain;
import ca.mcmaster.magarveylab.prism.fasta.FastaUtil;

public class ReductiveLoopAnalyzer {
	
	public static boolean isActiveKetoreductase(Domain domain, Cluster cluster) {
		boolean flag = false;
		if (domain.type() == ThiotemplatedDomains.KETOREDUCTASE) {
			String sequence = FastaUtil.getDomainSequence(domain, cluster);
			for (int i = 0; i < sequence.length() - 5; i++) {
				char g1 = sequence.charAt(i);
				char g2 = sequence.charAt(i+2);
				char g3 = sequence.charAt(i+5);
				if (g1 != 'G' && g1 != 'g')
					continue;
				if (g2 != 'G' && g2 != 'g')
					continue;
				if (g3 != 'G' && g3 != 'A'
						&& g3 != 'g' && g3 != 'a')
					continue;
				flag = true;
			}
		}
		return flag;
	}

}
