package ca.mcmaster.magarveylab.prism.enums.hmms;

import ca.mcmaster.magarveylab.enums.interfaces.SubstrateType;
import ca.mcmaster.magarveylab.enums.substrates.AlphaKetoAcids;
import ca.mcmaster.magarveylab.enums.substrates.DecoySubstrates;
import ca.mcmaster.magarveylab.enums.substrates.FungalAminoAcids;
import ca.mcmaster.magarveylab.enums.substrates.NonproteinogenicAminoAcids;
import ca.mcmaster.magarveylab.enums.substrates.ProteinogenicAminoAcids;

/**
 * Substrates, mostly alpha- or beta-amino acids, of adenylation domains.
 * 
 * @author skinnider
 *
 */
public enum AdenylationHmms implements SubstrateHmm {

	_2_AMINO_ADIPIC_ACID("Khayatt_ref_aad_1.hmm", NonproteinogenicAminoAcids._2_AMINO_ADIPIC_ACID),
	ALANINE_1("Alanine_1.hmm", ProteinogenicAminoAcids.ALANINE),
	ALANINE_2("Khayatt_ref_ala_2.hmm", ProteinogenicAminoAcids.ALANINE),
	ALANINE_3("Khayatt_ref_ala_3.hmm", ProteinogenicAminoAcids.ALANINE),
	BETA_ALANINE("Khayatt_ref_alabeta_1.hmm", NonproteinogenicAminoAcids.BETA_ALANINE),
	_3_HYDROXYANTHRANILATE("3-hydroxyanthranilate.hmm", NonproteinogenicAminoAcids._3_HYDROXYANTHRANILATE),
	AMINOCAPROLACTAM("aminocaprolactam.hmm", NonproteinogenicAminoAcids.AMINOCAPROLACTAM),
	ARGININE("Khayatt_ref_arg_1.hmm", ProteinogenicAminoAcids.ARGININE),
	ASPARAGINE_1("Khayatt_ref_asn_1.hmm", ProteinogenicAminoAcids.ASPARAGINE),
	ASPARAGINE_2("Khayatt_ref_asn_2.hmm", ProteinogenicAminoAcids.ASPARAGINE),
	BETA_HYDROXY_ASPARAGINE("beta_hydroxy_asparagine.hmm", NonproteinogenicAminoAcids.BETA_HYDROXY_ASPARAGINE),
	BETA_HYDROXY_ASPARTATE("beta_hydroxy_aspartate.hmm", NonproteinogenicAminoAcids.BETA_HYDROXY_ASPARTATE),
	ASPARTIC_ACID("Khayatt_ref_asp_1.hmm", ProteinogenicAminoAcids.ASPARTATE),
	METHYL_ASPARTIC_ACID("Khayatt_ref_aspme_1.hmm", NonproteinogenicAminoAcids.METHYL_ASPARTIC_ACID),
	BETA_HYDROXYTYROSINE("Khayatt_ref_bh_1.hmm", NonproteinogenicAminoAcids.BETA_HYDROXYTYROSINE),
	_4R_4E_2_BUTENYL_4_METHYL_L_THREONINE("Khayatt_ref_bm_1.hmm", NonproteinogenicAminoAcids._4R_4E_2_BUTENYL_4_METHYL_L_THREONINE),
	CAPREOMYCIDINE("capreomycidine.hmm", NonproteinogenicAminoAcids.CAPREOMYCIDINE),
	CYSTEINE_1("Khayatt_ref_cys_1.hmm", ProteinogenicAminoAcids.CYSTEINE),
	CYSTEINE_2("Khayatt_ref_cys_2.hmm", ProteinogenicAminoAcids.CYSTEINE),
	DIAMINOPROPIONATE("diaminopropionate.hmm", NonproteinogenicAminoAcids.DIAMINOPROPIONATE),
	DIAMINOPROPIONATE_2("Diaminopropanoic_acid_2.hmm", NonproteinogenicAminoAcids.DIAMINOPROPIONATE),
	_2_4_DIAMINO_BUTYRIC_ACID("Khayatt_ref_dab_1.hmm", NonproteinogenicAminoAcids._2_4_DIAMINO_BUTYRIC_ACID),
	DHAB_AKA_DEHYDROTHREONINE("Khayatt_ref_dhab_dh_1.hmm", NonproteinogenicAminoAcids.DEHYDROAMINOBUTYRIC_ACID), 
	DEHYDROAMINOBUTYRIC_ACID("dehydroaminobutyric_acid.hmm", NonproteinogenicAminoAcids.DEHYDROAMINOBUTYRIC_ACID),
	_3_5_DIHYDROXYPHENYLGLYCINE("Khayatt_ref_dhpg_dpg_1.hmm",  NonproteinogenicAminoAcids._3_5_DIHYDROXYPHENYLGLYCINE),
	GLUTAMINE("Khayatt_ref_gln_1.hmm", ProteinogenicAminoAcids.GLUTAMINE), 
	GLUTAMATE("Khayatt_ref_glu_1.hmm", ProteinogenicAminoAcids.GLUTAMATE),
	GLUTAMATE_2("Glutamic_acid_2.hmm", ProteinogenicAminoAcids.GLUTAMATE),
	METHYL_GLUTAMATE("methyl-glutamate.hmm", NonproteinogenicAminoAcids.METHYL_GLUTAMATE),
	GLYCINE_1("Khayatt_ref_gly_1.hmm", ProteinogenicAminoAcids.GLYCINE),
	GLYCINE_2("Khayatt_ref_gly_2.hmm", ProteinogenicAminoAcids.GLYCINE),
	HISTIDINE_1("Khayatt_ref_his_1.hmm", ProteinogenicAminoAcids.HISTIDINE),
	HISTIDINE_2("histidine.hmm", ProteinogenicAminoAcids.HISTIDINE),
	HISTIDINE_3("histidine_3.hmm", ProteinogenicAminoAcids.HISTIDINE),
	_4_HYDROXY_PHENYLGLYCINE("Khayatt_ref_hpg_hpg2Cl_1.hmm", NonproteinogenicAminoAcids._4_HYDROXY_PHENYLGLYCINE),
	HYDROXYVALINE("Khayatt_ref_hyv_1.hmm", NonproteinogenicAminoAcids.HYDROXYVALINE),
	ISOLEUCINE("Khayatt_ref_ile_1.hmm", ProteinogenicAminoAcids.ISOLEUCINE),
	AMINOBUTYRIC_ACID("Khayatt_ref_iva_abu_1.hmm", NonproteinogenicAminoAcids.AMINOBUTYRIC_ACID),
	BETA_HYDROXY_LEUCINE("beta_hydroxy_leucine.hmm", NonproteinogenicAminoAcids.BETA_HYDROXY_LEUCINE),
	LEUCINE_1("Khayatt_ref_leu_1.hmm", ProteinogenicAminoAcids.LEUCINE),
	LEUCINE_2("Khayatt_ref_leu_2.hmm", ProteinogenicAminoAcids.LEUCINE),
	LEUCINE_3("Khayatt_ref_leu_3.hmm", ProteinogenicAminoAcids.LEUCINE),
	LYSINE("Khayatt_ref_lys_1.hmm", ProteinogenicAminoAcids.LYSINE),
	BETA_LYSINE("Khayatt_ref_lysbeta_1.hmm", NonproteinogenicAminoAcids.BETA_LYSINE),
	METHIONINE("methionine.hmm", ProteinogenicAminoAcids.METHIONINE),
	ORNITHINE("Khayatt_ref_orn_1.hmm", NonproteinogenicAminoAcids.ORNITHINE),
	N5_HYDROXYORNITHINE("Khayatt_ref_ornfn5h_1.hmm", NonproteinogenicAminoAcids.N5_HYDROXYORNITHINE),
	HYDROXY_ACETYL_ORNITHINE("Khayatt_ref_ornha_1.hmm", NonproteinogenicAminoAcids.HYDROXY_ACETYL_ORNITHINE),
	BETA_PHENYLALANINE("beta_phenylalanine.hmm", NonproteinogenicAminoAcids.BETA_PHENYLALANINE),
	BETA_METHYL_PHENYLALANINE("beta_methyl_phenylalanine.hmm", NonproteinogenicAminoAcids.BETA_METHYL_PHENYLALANINE),
	BETA_HYDROXY_PHENYLALANINE("beta_hydroxy_phenylalanine.hmm", NonproteinogenicAminoAcids.BETA_HYDROXY_PHENYLALANINE),
	PHENYLALANINE("Khayatt_ref_phe_1.hmm", ProteinogenicAminoAcids.PHENYLALANINE),
	PIPECOLIC_ACID("Khayatt_ref_pip_1.hmm", NonproteinogenicAminoAcids.PIPECOLIC_ACID),
	HYDROXYPIPECOLIC_ACID("pipecolic_acid_2.hmm", NonproteinogenicAminoAcids.HYDROXYPIPECOLIC_ACID),
	PIPERAZIC_ACID("piperazic_acid.hmm", NonproteinogenicAminoAcids.PIPERAZIC_ACID),
	PROLINE_1("Khayatt_ref_pro_1.hmm", ProteinogenicAminoAcids.PROLINE),
	PROLINE_2("Khayatt_ref_pro_2.hmm", ProteinogenicAminoAcids.PROLINE),
	PROLINE_3("proline_3.hmm", ProteinogenicAminoAcids.PROLINE),
	METHYL_PROLINE("Khayatt_ref_prom_1.hmm", NonproteinogenicAminoAcids.METHYL_PROLINE),
	SERINE("Khayatt_ref_ser_1.hmm", ProteinogenicAminoAcids.SERINE),
	THREONINE_1("Khayatt_ref_thr_1.hmm", ProteinogenicAminoAcids.THREONINE),
	THREONINE_2("Khayatt_ref_thr_2.hmm", ProteinogenicAminoAcids.THREONINE),
	ALLO_THREONINE("Khayatt_ref_thrallo_1.hmm", NonproteinogenicAminoAcids.ALLO_THREONINE),
	TRYPTOPHAN("Khayatt_ref_trp_1.hmm", ProteinogenicAminoAcids.TRYPTOPHAN),
	TYROSINE_1("Khayatt_ref_tyr_1.hmm", ProteinogenicAminoAcids.TYROSINE),
	TYROSINE_2("Khayatt_ref_tyr_2.hmm", ProteinogenicAminoAcids.TYROSINE),
	VALINE_1("Khayatt_ref_val_1.hmm", ProteinogenicAminoAcids.VALINE),
	VALINE_2("Khayatt_ref_val_2.hmm", ProteinogenicAminoAcids.VALINE),
	VALINE_3("Khayatt_ref_val_3.hmm", ProteinogenicAminoAcids.VALINE),
	
	// decoy hmms 
	QUINOMYCIN_STARTER_UNIT_BIOSYNTHESIS("hydroxyquinaldic_acid_biosynthesis.hmm", 
			DecoySubstrates.QUINOMYCIN_STARTER_UNIT_BIOSYNTHESIS),
	GLYCOPEPTIDE_STARTER_UNIT_BIOSYNTHESIS("beta_hydroxytyrosine_biosynthesis.hmm", 
			DecoySubstrates.GLYCOPEPTIDE_STARTER_UNIT_BIOSYNTHESIS),
			
	// fungal  
	FUNGAL_4_HYDROXYPHENYLPYRUVATE("fungal_4-Hydroxyphenylpyruvicacid.hmm", FungalAminoAcids.HYDROXYPHENYLPYRUVATE),
	FUNGAL_ALANINE("fungal_alanine.hmm", ProteinogenicAminoAcids.ALANINE),
	// FUNGAL_AMINOBUTYRIC_ACID("fungal_Abu.hmm", FungalAminoAcids.AMINOBUTYRIC_ACID), 
	FUNGAL_AMINOISOBUTYRIC_ACID("fungal_Aib.hmm", FungalAminoAcids.AMINOISOBUTYRIC_ACID),
	FUNGAL_ANTHRANILIC_ACID("fungal_Anthrinilic_acid.hmm", FungalAminoAcids.ANTHRANILATE),
	FUNGAL_BETA_ALANINE("fungal_Beta-alanine.hmm", NonproteinogenicAminoAcids.BETA_ALANINE),
	FUNGAL_BETA_AMINOISOBUTYRIC_ACID("fungal_Beta-Aib.hmm", FungalAminoAcids.BETA_AMINOISOBUTYRIC_ACID),
	FUNGAL_DEHYDROALANINE("fungal_Dehydro-ala.hmm", FungalAminoAcids.DEHYDROALANINE),
	FUNGAL_GLUTAMINE("fungal_Glutamine.hmm", ProteinogenicAminoAcids.GLUTAMINE),
	FUNGAL_HYDROXYGLUTAMINE("fungal_Hydroxyglutamine.hmm", FungalAminoAcids.HYDROXYGLUTAMINE),
	FUNGAL_GLYCINE("fungal_Glycine.hmm", ProteinogenicAminoAcids.GLYCINE),
	FUNGAL_HYDROXYISOCAPROIC_ACID("fungal_Hic.hmm", AlphaKetoAcids.ALPHA_KETOISOCAPROATE),
	FUNGAL_HYDROXYISOVALERIC_ACID("fungal_Hiv.hmm", AlphaKetoAcids.ALPHA_KETOISOVALERATE),
	FUNGAL_HYDROXYMETHYLPENTANOIC_ACID("fungal_HMP.hmm", AlphaKetoAcids._3_METHYL_2_OXOPENTANOIC_ACID),
	FUNGAL_HOMOSERINE("fungal_Homo-serine.hmm", FungalAminoAcids.HOMOSERINE),
	FUNGAL_HOMOTYROSINE("fungal_Homo-tyrosine.hmm", FungalAminoAcids.HOMOTYROSINE),
	FUNGAL_HYDROXYHOMOTYROSINE_SULFATE("fungal_Hydroxy-homotyrosine-sulfate.hmm", FungalAminoAcids.HYDROXYHOMOTYROSINE_SULFATE),
	FUNGAL_INDOLE_PYRUVIC_ACID("fungal_IPA.hmm", FungalAminoAcids.INDOLE_PYRUVIC_ACID),
	FUNGAL_ISOLEUCINE("fungal_Ile.hmm", ProteinogenicAminoAcids.ISOLEUCINE),
	FUNGAL_LEUCINE("fungal_Leu.hmm", ProteinogenicAminoAcids.LEUCINE),
	FUNGAL_N_METHOXY_TRYPTOPHAN("fungal_N-Methoxy-tryptophan.hmm", FungalAminoAcids.N_METHOXY_TRYPTOPHAN),
	FUNGAL_ORNITHINE("fungal_Orn.hmm", NonproteinogenicAminoAcids.ORNITHINE),
	FUNGAL_PHENYLALANINE("fungal_Phe.hmm", ProteinogenicAminoAcids.PHENYLALANINE),
	FUNGAL_PIPECOLIC_ACID("fungal_Pipecolic_acid.hmm", NonproteinogenicAminoAcids.PIPECOLIC_ACID),
	FUNGAL_PROLINE("fungal_Pro.hmm", ProteinogenicAminoAcids.PROLINE),
	FUNGAL_HYDROXYPROLINE("fungal_OHPro.hmm", FungalAminoAcids.HYDROXYPROLINE),
	FUNGAL_METHYLPROLINE("fungal_MePro.hmm", NonproteinogenicAminoAcids.METHYL_PROLINE),
	FUNGAL_SERINE("fungal_Ser.hmm", ProteinogenicAminoAcids.SERINE),
	FUNGAL_THREONINE("fungal_Thr.hmm", ProteinogenicAminoAcids.THREONINE),
	FUNGAL_TRYPTOPHAN("fungal_Trp.hmm", ProteinogenicAminoAcids.TRYPTOPHAN),
	FUNGAL_TYROSINE("fungal_Tyr.hmm", ProteinogenicAminoAcids.TYROSINE),
	FUNGAL_VALINE("fungal_Val.hmm", ProteinogenicAminoAcids.VALINE),
	;
	
	private final String hmm;
	private final SubstrateType substrate;
	
	private AdenylationHmms(final String hmm, final SubstrateType substrate) {
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
