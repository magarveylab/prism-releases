package ca.mcmaster.magarveylab.prism.enums.hmms;

import ca.mcmaster.magarveylab.enums.interfaces.SubstrateType;
import ca.mcmaster.magarveylab.enums.substrates.PolyketideMonomers;

/**
 * HMMs used to predict acyltransferase domain substrate specificity.
 * 
 * @author skinnider
 */
public enum AcyltransferaseHmms implements SubstrateHmm {

	_2_METHYLBUTERYL_COA("2MBuC_1.hmm", PolyketideMonomers._2_METHYLBUTERYL_COA),
	BENZOYL_COA("BzC_1.hmm", PolyketideMonomers.BENZOYL_COA),
	ETHYLMALONYL_COA_1("EMC_1.hmm", PolyketideMonomers.ETHYLMALONYL_COA),
	ETHYLMALONYL_COA_2("EMC_2.hmm", PolyketideMonomers.ETHYLMALONYL_COA),
	ETHYLMALONYL_COA_3("EMC_3.hmm", PolyketideMonomers.ETHYLMALONYL_COA),
	ISOBUTERYL_COA("IBuC_1.hmm", PolyketideMonomers.ISOBUTERYL_COA),
	MALONYL_COA_1("MC_1.hmm", PolyketideMonomers.MALONYL_COA),
	MALONYL_COA_2("MC_2.hmm", PolyketideMonomers.MALONYL_COA),
	METHYLMALONYL_COA_1("MMC_1.hmm", PolyketideMonomers.METHYLMALONYL_COA),
	METHYLMALONYL_COA_2("MMC_2.hmm", PolyketideMonomers.METHYLMALONYL_COA),
	METHYLMALONYL_COA_3("MMC_3.hmm", PolyketideMonomers.METHYLMALONYL_COA),
	METHYLMALONYL_COA_4("MMC_4.hmm", PolyketideMonomers.METHYLMALONYL_COA),
	METHOXYLMALONYL_COA_1("MOMC_1.hmm", PolyketideMonomers.METHOXYLMALONYL_COA),
	METHOXYLMALONYL_COA_2("MOMC_2.hmm", PolyketideMonomers.METHOXYLMALONYL_COA),
	PROPIONYL_COA_1("PC_1.hmm", PolyketideMonomers.PROPIONYL_COA);  

	private final String hmm;
	private final SubstrateType substrate;
	
	private AcyltransferaseHmms(final String hmm, final SubstrateType substrate) {
		this.hmm = hmm;
		this.substrate = substrate;
	}

	public String hmm() { 
		return hmm; 
	}
	
	public String fullName() { 
		return substrate.fullName(); 
	}
	
	public String abbreviation() { 
		return substrate.abbreviation(); 
	}
	
	public String smiles() { 
		return substrate.smiles(); 
	}
	
}