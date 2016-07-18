package ca.mcmaster.magarveylab.prism.enums.hmms;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.NucleosideDomains;
import ca.mcmaster.magarveylab.enums.domains.PrerequisiteDomains;
import ca.mcmaster.magarveylab.enums.domains.TypeIIPolyketideDomains;
import ca.mcmaster.magarveylab.enums.interfaces.SubstrateType;

/**
 * Substrates which require the detection of one or more additional biosynthetic
 * enzymes within the cluster to be conclusively identified.
 * 
 * @author skinnider
 *
 */
public enum SubstratePrerequisites {

	METHOXYMALONATE_1(AcyltransferaseHmms.METHOXYLMALONYL_COA_1, new DomainType[] { 
			PrerequisiteDomains._3_HYDROXYACYL_COA_DEHYDROGENASE,
			PrerequisiteDomains.ACYL_COA_DEHYDROGENASE }),
	METHOXYMALONATE_2(AcyltransferaseHmms.METHOXYLMALONYL_COA_2, new DomainType[] { 
			PrerequisiteDomains._3_HYDROXYACYL_COA_DEHYDROGENASE,
			PrerequisiteDomains.ACYL_COA_DEHYDROGENASE }),
/*	ENDURACIDIDINE(AdenylationSubstrates.ENDURACIDIDINE, new DomainType[] { 
			PrerequisiteDomains.PLP_DEPENDENT_AMINOTRANSFERASE_1, 
			PrerequisiteDomains.PLP_DEPENDENT_AMINOTRANSFERASE_2 }), */
	DIAMINOPROPIONATE(AdenylationHmms.DIAMINOPROPIONATE, new DomainType[] { 
			PrerequisiteDomains.DIAMINOPROPIONATE_SYNTHASE }),
	CAPREOMYCIDINE(AdenylationHmms.CAPREOMYCIDINE, new DomainType[] { 
			PrerequisiteDomains.CAPREOMYCIDINE_SYNTHASE_1,
			PrerequisiteDomains.CAPREOMYCIDINE_SYNTHASE_2 }),
	PIPECOLIC_ACID(AdenylationHmms.HYDROXYPIPECOLIC_ACID, new DomainType[] {
			PrerequisiteDomains.PIPECOLIC_ACID_SYNTHASE_1,
			PrerequisiteDomains.PIPECOLIC_ACID_SYNTHASE_2 }),
	BETA_PHENYLALANINE(AdenylationHmms.BETA_PHENYLALANINE, new DomainType[] { 
			TypeIIPolyketideDomains.PAL }),
	_3_HYDROXYANTHRANILATE(AdenylationHmms._3_HYDROXYANTHRANILATE, new DomainType[] { 
			PrerequisiteDomains.ARYL_FORMAMIDASE,
			PrerequisiteDomains.TRYPTOPHAN_DIOXYGENASE }),
	AMINOCAPROLACTAM(AdenylationHmms.AMINOCAPROLACTAM, new DomainType[] { 
			NucleosideDomains.CapW }),
	HISTIDINE_3(AdenylationHmms.HISTIDINE_3, new DomainType[] { 
			NucleosideDomains.NikQ }),
//	HOMOTYROSINE(AdenylationSubstrates.FUNGAL_HOMOTYROSINE, new DomainType[] { 
//			PrerequisiteDomains.ISOPROPYL_MALATE_SYNTHASE, PrerequisiteDomains.TRANSAMINASE }),
	;
	
	private SubstrateType substrate;
	private DomainType[] domains;
	
	private SubstratePrerequisites(SubstrateType substrate, DomainType[] domains) {
		this.substrate = substrate;
		this.domains = domains;
	}
	
	public SubstrateType substrate() {
		return substrate;
	}
	
	public DomainType[] domains() {
		return domains;
	}
	
}
