package ca.mcmaster.magarveylab.prism.enums.hmms;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.ThiotemplatedDomains;

public enum SubstrateDomainSearches {

	ADENYLATION(ThiotemplatedDomains.ADENYLATION, AdenylationHmms.values()),
	ACYL_ADENYLATING(ThiotemplatedDomains.ACYL_ADENYLATING, AcylAdenylatingHmms.values()),
	ACYLTRANSFERASE(ThiotemplatedDomains.ACYLTRANSFERASE, AcyltransferaseHmms.values()),
	;
	
	private DomainType type;
	private SubstrateHmm[] substrates;
	
	private SubstrateDomainSearches(DomainType type, SubstrateHmm[] substrates) {
		this.type = type;
		this.substrates = substrates;
	}
	
	public DomainType type() {
		return type;
	}
	
	public SubstrateHmm[] substrates() {
		return substrates;
	}
	
}
