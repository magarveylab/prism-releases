package ca.mcmaster.magarveylab.prism.enums.hmms;

import ca.mcmaster.magarveylab.enums.interfaces.SubstrateType;
import ca.mcmaster.magarveylab.enums.substrates.AlphaKetoAcids;
import ca.mcmaster.magarveylab.enums.substrates.FattyAcids;
import ca.mcmaster.magarveylab.enums.substrates.StarterUnits;

/**
 * Substrates of a distinct clade of acyl-adenylating enzymes revealed by
 * phylogenetic analysis (see doi:10.1093/nar/gkv1012).
 * 
 * @author skinnider
 *
 */
public enum AcylAdenylatingHmms implements SubstrateHmm {
	
	// fatty acids
	MYRISTATE("myristate.hmm", FattyAcids.MYRISTATE),
	LONG_CHAIN_FATTY_ACID("long_chain_fatty_acid.hmm", FattyAcids.LONG_CHAIN_FATTY_ACID), 
	SHORT_CHAIN_FATTY_ACID("short_chain_fatty_acid.hmm", FattyAcids.SHORT_CHAIN_FATTY_ACID),
	_3_AMINONON_5_ENOIC_ACID("3-aminonon-5-enoic_acid.hmm", FattyAcids._3_AMINONON_5_ENOIC_ACID),

	// aromatic starters
	_2_3_DIHYDROXYBENZOIC_ACID("2_3_dihydroxybenzoic_acid.hmm", StarterUnits._2_3_DIHYDROXYBENZOIC_ACID),
	SALICYLIC_ACID("salicylic_acid.hmm", StarterUnits.SALICYLIC_ACID),
	CINNAMIC_ACID("cinnamic_acid.hmm", StarterUnits.CINNAMIC_ACID),
	_3_AMINO_5_HYDROXYBENZOIC_ACID("AHBA.hmm", StarterUnits._3_AMINO_5_HYDROXYBENZOIC_ACID),
	_3_FORMAMIDO_5_HYDROXYBENZOIC_ACID("FHBA.hmm", StarterUnits._3_FORMAMIDO_5_HYDROXYBENZOIC_ACID),
	_3_HYDROXYPICOLINIC_ACID("3-HPA.hmm", StarterUnits._3_HYDROXYPICOLINIC_ACID),
	_3_HYDROXYQUINALDIC_ACID("3-HQA.hmm", StarterUnits._3_HYDROXYQUINALDIC_ACID),
	QUINOXALINE("quinoxaline.hmm", StarterUnits.QUINOXALINE),
	PHENYLACETATE("phenylacetate.hmm", StarterUnits.PHENYLACETATE),
	PARA_HYDROXYBENZOIC_ACID("PHBA.hmm", StarterUnits.PARA_HYDROXYBENZOIC_ACID),
	PARA_AMINOBENZOIC_ACID("PABA.hmm", StarterUnits.PARA_AMINOBENZOIC_ACID),
	
	// small starters
	BETA_AMINOALANINAMIDE("beta_aminoalaninamide.hmm", StarterUnits.BETA_AMINOALANINAMIDE),
	CYCLOHEXANECARBOXYLATE("CHC.hmm", StarterUnits.CYCLOHEXANECARBOXYLATE),
	DIHYDROXYCYCLOHEXANECARBOXYLIC_ACID("DHCHC.hmm", StarterUnits.DIHYDROXYCYCLOHEXANECARBOXYLIC_ACID),
	ALKYLMALONYL_COA("alkylmalonyl_CoA.hmm", StarterUnits.ALKYLMALONYL_COA),
	_3_HYDROXYBUTANOIC_ACID("3-hydroxybutanoic_acid.hmm", StarterUnits._3_HYDROXYBUTANOIC_ACID),
	AMINOLEVULINIC_ACID("aminolevulinic_acid.hmm", StarterUnits.AMINOLEVULINIC_ACID),
	
	// alpha keto/alpha hydroxy acids
	PYRUVATE("pyruvate.hmm", AlphaKetoAcids.PYRUVATE),
	ALPHA_KETOISOVALERATE("alpha-ketoisovalerate.hmm", AlphaKetoAcids.ALPHA_KETOISOVALERATE),
	ALPHA_KETOISOCAPROATE("alpha-ketoisocaproate.hmm", AlphaKetoAcids.ALPHA_KETOISOCAPROATE),
	_3_METHYL_2_OXOPENTANOIC_ACID("3-methyl-2-oxopentanoate.hmm", AlphaKetoAcids._3_METHYL_2_OXOPENTANOIC_ACID),
	PHENYLPYRUVATE("phenylpyruvate.hmm", AlphaKetoAcids.PHENYLPYRUVATE),
	
	//amni - fungal additions
	//TODO need to add SMILES
//	_4_HYDROXYPHENYLPYRUVIC_ACID("4-Hydroxyphenylpyruvicacid", "4-Hydroxyphenylpyruvicacid.hmm", "4Hpp", ""),
//	FATTY_ACID_DERIVED_AMINO_ACID("Fatty-acid-derived_Amino_acid.hmm", "Fatty-acid-derived_Amino_acid.hmm", "Faa", ""),
//	HYDROXYISOCAPROIC("Hydroxyisocaproic acid", "Hydroxyisocaproic_acid.hmm", "Hic", ""),
//	HYDROXYISOVALERIC("Hydroxyisovaleric acid", "Hydroxyisovaleric_acid.hmm", "Hiv", ""),
//	_2_HYDROXY_3_METHYLPENTANOIC_ACID("2-hydroxy-3-methylpentanoic acid", "2-hydroxy-3-methylpentanoic_acid.hmm", "Hmp", ""),
//	INDOLE_PYRUVIC_ACID("Indole-pyruvic acid", "Indole-pyruvic_acid.hmm", "Ipa", ""),
//	LINOLEIC_ACID("Linoleic acid", "Linoleic_acid.hmm", "Lin", ""),
//	METHYLATED_MYRISTIC_ACID("Methylated-myristic acid", "Methylated-myristic_acid.hmm", "Mma", ""),
	;
	
	private final String hmm;
	private final SubstrateType substrate;
	
	private AcylAdenylatingHmms(final String hmm, final SubstrateType substrate) {
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
