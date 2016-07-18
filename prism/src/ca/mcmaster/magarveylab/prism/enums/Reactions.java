package ca.mcmaster.magarveylab.prism.enums;

import ca.mcmaster.magarveylab.prism.cluster.annotation.Annotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.impl.*;
import ca.mcmaster.magarveylab.prism.cluster.annotation.typeII.AllTypeIIModulesAnnotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.typeII.AntibioticMonooxygenaseAnnotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.typeII.BVMOAnnotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.typeII.C11C12Annotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.typeII.C15C16Annotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.typeII.C17C18Annotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.typeII.C19C20Annotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.typeII.C1C2Annotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.typeII.C1DRingCyclaseAnnotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.typeII.C2AnthracyclineDRingCyclaseAnnotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.typeII.C2CDRingCyclaseAnnotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.typeII.C2DERingCyclaseAnnotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.typeII.C2TetracenomycinDRingCyclaseAnnotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.typeII.C3CRingCyclaseAnnotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.typeII.C5C6Annotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.typeII.C5CRingCyclaseAnnotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.typeII.C7ABCRingCyclaseAnnotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.typeII.C7ABRingCyclaseAnnotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.typeII.C7C8Annotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.typeII.C9C10Annotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.typeII.C9C14FirstRingCyclaseAnnotator;
import ca.mcmaster.magarveylab.prism.cluster.annotation.typeII.CGlycosyltransferaseAnnotator;
import ca.mcmaster.magarveylab.prism.cluster.reactions.GenericReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.impl.*;
import ca.mcmaster.magarveylab.prism.cluster.reactions.typeII.ABCRingCyclaseReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.typeII.ABRingCyclaseReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.typeII.AminotransferaseReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.typeII.AngucyclineCyclaseReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.typeII.AnthracyclineFourthRingCyclaseReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.typeII.AntibioticMonooxygenaseReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.typeII.BVMOReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.typeII.C2MethyltransferaseReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.typeII.CGlycosyltransferaseReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.typeII.FavorskiiaseReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.typeII.PentangularPolyphenolCyclaseReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.typeII.PyranCyclaseReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.typeII.TetracenomycinFourthRingCyclaseReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.typeII.TetracyclineFourthRingCyclaseReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.typeII.ThirdRingCyclaseReaction;
import ca.mcmaster.magarveylab.prism.cluster.reactions.typeII.TypeIIKetoreductaseReaction;
import ca.mcmaster.magarveylab.prism.util.exception.ClassInstantiationException;
import ca.mcmaster.magarveylab.prism.cluster.reactions.ribosomal.*;
import ca.mcmaster.magarveylab.prism.cluster.annotation.ribosomal.*;

public enum Reactions {
	
	// type II polyketide cyclases and beta-lactam formation
	ISOPENICILLIN_N_SYNTHASE(IsopenicillinNSynthaseReaction.class, 
			new Annotator[] { new IsopenicillinNSynthaseAnnotator() }),
	DEACETOXYCEPHALOSPORIN_C_SYNTHASE(DeacetoxycephalosporinCSynthaseReaction.class, 
			new Annotator[] { new DeacetoxycephalosporinCSynthaseAnnotator() }),
	ABC_RING_CYCLASE(ABCRingCyclaseReaction.class,
			new Annotator[] { new C7ABCRingCyclaseAnnotator(), new C9C14FirstRingCyclaseAnnotator() }),
	AB_RING_CYCLASE(ABRingCyclaseReaction.class, new Annotator[] { new C7ABRingCyclaseAnnotator() }),
	THIRD_RING_CYCLASE(ThirdRingCyclaseReaction.class, new Annotator[] { new C3CRingCyclaseAnnotator(),
			new C5CRingCyclaseAnnotator() }),
	PYRAN_CYCLASE(PyranCyclaseReaction.class, new Annotator[] { new C3CRingCyclaseAnnotator() }),
	ANGCYCLINE_CYCLASE(AngucyclineCyclaseReaction.class, new Annotator[] { new C2CDRingCyclaseAnnotator() }),
	TETRACENOMYCIN_FOURTH_RING_CYCLASE(TetracenomycinFourthRingCyclaseReaction.class, 
			new Annotator[] { new C2TetracenomycinDRingCyclaseAnnotator() }),
	ANTHRACYCLINE_FOURTH_RING_CYCLASE(AnthracyclineFourthRingCyclaseReaction.class, 
			new Annotator[] { new C2AnthracyclineDRingCyclaseAnnotator() }),
	TETRACYCLINE_FOURTH_RING_CYCLASE(TetracyclineFourthRingCyclaseReaction.class, 
			new Annotator[] { new C1DRingCyclaseAnnotator() }),
	PENTANGULAR_POLYPHENOL_CYCLASE(PentangularPolyphenolCyclaseReaction.class, 
			new Annotator[] { new C2DERingCyclaseAnnotator() }),
	FAVORSKIIASE(FavorskiiaseReaction.class, new Annotator[] { new AllTypeIIModulesAnnotator() }),
	BVMO(BVMOReaction.class, new Annotator[] { new BVMOAnnotator() }),
			
	// other multi-modular tailoring enzymes
	P450A(P450AReaction.class, new Annotator[] { new P450AAnnotator() }),
	P450B(P450BReaction.class, new Annotator[] { new P450BAnnotator() }),
	P450C(P450CReaction.class, new Annotator[] { new P450CAnnotator() }),
	P450D(P450DReaction.class, new Annotator[] { new P450DAnnotator() }),
	ABM(AntibioticMonooxygenaseReaction.class, new Annotator[] { new AntibioticMonooxygenaseAnnotator() }),

	// oxazole/thiazole formation
	CONDENSATION_CYCLIZATION(HeterocyclizationReaction.class, 
			new Annotator[] { new HeterocyclizationAnnotator() }),
	NITROREDUCTASE(NitroreductaseReaction.class, new Annotator[] { new ModularAnnotator() }),
			
	// oxidoreductases
	PROLINE_DEHYDROGENASE(ProlineDehydrogenaseReaction.class, new Annotator[] { new ProlineDehydrogenaseAnnotator() }),
	TRYPTOPHAN_DIOXYGENASE(TryptophanDioxygenaseReaction.class, new Annotator[] { new TryptophanDioxygenaseAnnotator() }),
	KETOREDUCTASE(KetoreductaseReaction.class, new Annotator[] { new ReductiveLoopAnnotator() }),
	DEHYDRATASE(DehydrataseReaction.class, new Annotator[] { new ReductiveLoopAnnotator() }),
	ENOLREDUCTASE(EnolreductaseReaction.class, new Annotator[] { new ReductiveLoopAnnotator() }),
	C9_KETOREDUCTASE(TypeIIKetoreductaseReaction.class, new Annotator[] { new C9C10Annotator() }),
	C11_KETOREDUCTASE(TypeIIKetoreductaseReaction.class, new Annotator[] { new C11C12Annotator() }),
	C15_KETOREDUCTASE(TypeIIKetoreductaseReaction.class, new Annotator[] { new C15C16Annotator() }),
	C17_KETOREDUCTASE(TypeIIKetoreductaseReaction.class, new Annotator[] { new C17C18Annotator() }),
	C19_KETOREDUCTASE(TypeIIKetoreductaseReaction.class, new Annotator[] { new C19C20Annotator() }),
			
	// modular tailoring reactions
	O_METHYLTRANSFERASE(OMethyltransferaseReaction.class, new Annotator[] { new ModularAnnotator() }),
	C_METHYLTRANSFERASE(CMethyltransferaseReaction.class, new Annotator[] { new ModularAnnotator() }),
	N_METHYLTRANSFERASE(NMethyltransferaseReaction.class, new Annotator[] { new ModularAnnotator() }),

	// simple addition reactions
	CARBAMOYLTRANSFERASE(CarbamoyltransferaseReaction.class, new Annotator[] { new HydroxylAnnotator() }),
	PHOSPHOTRANSFERASE(PhosphotransferaseReaction.class, new Annotator[] { new HydroxylAnnotator() }),
	FORMYLTRANSFERASE(FormylationReaction.class, new Annotator[] { new FormyltransferaseAnnotator() }),
	SULFOTRANSFERASE(SulfotransferaseReaction.class, new Annotator[] { new HydroxylAnnotator() }),
	ISOPENICILLIN_N_ACYLTRANSFERASE(IsopenicillinNAcyltransferaseReaction.class, 
			new Annotator[] { new IsopenicillinNAcyltransferaseAnnotator() }),
	CHLORINATION(ChlorinaseReaction.class, new Annotator[] { new ChlorinaseAnnotator() }),
	GLYCOSYLTRANSFERASE(GlycosyltransferaseReaction.class, new Annotator[] { new HydroxylAnnotator() }),
	C_GLYCOSYLTRANSFERASE(CGlycosyltransferaseReaction.class, new Annotator[] { new CGlycosyltransferaseAnnotator() }),
	ACYL_ADENYLATE_LIGASE(AcylAdenylateLigaseReaction.class, new Annotator[] { new AcylAdenylateLigaseAnnotator() }),
	C2_AMINOTRANSFERASE(AminotransferaseReaction.class, new Annotator[] { new C1C2Annotator() }),
	CARBOXY_METHYLTRANSFERASE(OMethyltransferaseReaction.class, new Annotator[] { new C1C2Annotator() }),
	C11_OMT(OMethyltransferaseReaction.class, new Annotator[] { new C11C12Annotator() }),
	C2MT(C2MethyltransferaseReaction.class, new Annotator[] { new C1C2Annotator() }),
	C6_CMT(CMethyltransferaseReaction.class, new Annotator[] { new C5C6Annotator() }),
	C8_CMT(CMethyltransferaseReaction.class, new Annotator[] { new C7C8Annotator() }),
	C10_CMT(CMethyltransferaseReaction.class, new Annotator[] { new C9C10Annotator() }),
	
	//ribosomal reactions
	AGEF1(AgeF1Reaction.class, new Annotator[] { new MacrocyclizationAnnotator() } ),
	AGRB(AgrBReaction.class, new Annotator[] { new AgrBAnnotator() } ),
	AMINOVINYLCYSTEINE(AminovinylcysteineReaction.class, new Annotator[] { new AminovinylcysteineAnnotator() } ),
	AZOLE(AzoleReaction.class, new Annotator[] { new AzoleAnnotator(), new YmB1Annotator() } ),
	AZOLINE(AzolineReaction.class, new Annotator[] { new AzoleAnnotator(), new YmB1Annotator() } ),
	BETAHYDROXYLASE(BetaHydroxylaseReaction.class, new Annotator[] { new AspartateAnnotator(), new PhenylalanineAnnotator(), new HistidineAnnotator(), new ValineAnnotator() } ),
	BETAMETHYLTRANSFERASE(BetaMethyltransferaseReaction.class, new Annotator[] { new PhenylalanineAnnotator(), new BotRMT2Annotator(), new ProlineAnnotator() } ),
	BOTOMT(BotOMTReaction.class, new Annotator[] { new AspartateAnnotator() } ),
	BOTP(BotPReaction.class, new Annotator[] { new BotPAnnotator() } ),
	CINORF7(Cinorf7Reaction.class, new Annotator[] { new Cinorf7Annotator() } ),
	CLTM(CltMReaction.class, new Annotator[] { new CltMAnnotator() } ),
	COMQ(ComQReaction.class, new Annotator[] { new TryptophanAnnotator() } ),
	DIHYDROXYLASE(DihydroxylaseReaction.class, new Annotator[] { new IsoleucineAnnotator() } ),
	DISULFIDE(DisulfideReaction.class, new Annotator[] { new LassoPeptideDisulfideAnnotator(), new SactipeptideDisulfideAnnotator()} ),
	ELXO(ElxOReaction.class, new Annotator[] { new NTerminusAnnotator() } ),
	GAMMAHYDROXYLASE(GammaHydroxylaseReaction.class, new Annotator[] { new GlutamateAnnotator()} ),
	GARO(GarOReaction.class, new Annotator[] { new LanthionineAnnotator() } ),
	HYDROXYMETHYLPROLINE(HydroxymethylprolineReaction.class, new Annotator[] { new IsoleucineAnnotator() } ),
	LANB(LanBReaction.class, new Annotator[] { new LanBAnnotator(), new ThreonineDehydrataseAnnotator(), new SerineDehydrataseAnnotator() } ),
	LANC(LanCReaction.class, new Annotator[] { new LanthionineAnnotator() } ),
	LANJ(LanJReaction.class, new Annotator[] { new LanthionineAnnotator() } ),
	LANK(LanKCReaction.class, new Annotator[] { new LabioninAnnotator() } ),
	LANM(LanMReaction.class, new Annotator[] { new LanMAnnotator() } ),
	LASSOPEPTIDE(LassoPeptideReaction.class, new Annotator[] { new LassoPeptideAnnotator() } ),
	MACROCYCLIZATION(MacrocyclizationReaction.class, new Annotator[] { new MacrocyclizationAnnotator() } ),
	MACROCYCLOHYDRATION(MacrocyclodehydrationReaction.class, new Annotator[] { new BotCAnnotator() } ),
	MDNB(MdnBReaction.class, new Annotator[] { new MdnBAnnotator() } ),
	MDNC(MdnCReaction.class, new Annotator[] { new MdnCAnnotator() } ),
	MIBH(MibHReaction.class, new Annotator[] { new TryptophanAnnotator() } ),
	MIBO(MibOReaction.class, new Annotator[] { new ProlineAnnotator() } ),
	NACETYLATION(NAcetylationReaction.class, new Annotator[] { new NTerminusAnnotator() } ),
	NNDIMETHYLTRANSFERASE(NNDimethyltransferaseReaction.class, new Annotator[] { new NTerminusAnnotator() } ),
	NOCQ(NocQReaction.class, new Annotator[] { new SerineAnnotator() } ),
	NOSA(NosAReaction.class, new Annotator[] { new CTerminusAnnotator() } ),
	NOSI(NosIReaction.class, new Annotator[] { new NosIAnnotator(), new TsrIAnnotator() } ),
	PATGOX(PatGOxReaction.class, new Annotator[] { new PatGOxAnnotator() } ),
	POYBC(BetaMethyltransferaseReaction.class, new Annotator[] { new PoyBCAnnotator() }),
	POYE(PoyEReaction.class, new Annotator[] { new PoyEAnnotator() }),
	POYF(PoyFReaction.class, new Annotator[] { new PoyFAnnotator() }),
	POYI(BetaHydroxylaseReaction.class, new Annotator[] { new PoyIAnnotator() }),
	PRENYLATION(PrenylationReaction.class, new Annotator[] { new PatFAnnotator() } ),
	PYRIDINEHYDROXYLATION(PyridineHydroxylationReaction.class, new Annotator[] { new PyridineAnnotator() } ),
	PYRIDINE(PyridineReaction.class, new Annotator[] { new PyridineAnnotator() } ),
	SACTIPEPTIDE(SactipeptideReaction.class, new Annotator[] { new SactipeptideAnnotator() } ), 
	STRB(StrBReaction.class, new Annotator[] { new StrBAnnotator() } ),
	SUNA(DisulfideReaction.class, new Annotator[] { new SunAAnnotator() } ),
	SUNS(SunSReaction.class, new Annotator[] { new SunSAnnotator() } ),
	TCLO(TclOReaction.class, new Annotator[] { new ThreonineAnnotator() } ),
	TFXB(TfxBReaction.class, new Annotator[] { new TfxBAnnotator() } ),
	THIOAMIDE(ThioamideReaction.class, new Annotator[] { new TvaHAnnotator() } ),
	TPAJ(TpaJReaction.class, new Annotator[] { new CTerminusAnnotator() }),
	TPDI(TpdIReaction.class, new Annotator[] { new CysteineAnnotator() } ),
	TPDJ12(TpdJ12Reaction.class, new Annotator[] { new TpdJ12Annotator() } ),
	TPDM(TpdMReaction.class, new Annotator[] { new CysteineAnnotator() } ),
	TPDT(TpdTReaction.class, new Annotator[] { new AsparagineAnnotator() } ),
	TSRB(TsrBReaction.class, new Annotator[] { new CTerminusAnnotator() } ),
	TSRC(TsrCReaction.class, new Annotator[] { new CTerminusAnnotator() } ),
	TSRI(TsrIReaction.class, new Annotator[] { new TsrIAnnotator() } ),
	TVAFAMINOVINYLCYSTEINE(TvaFAminovinylcysteineReaction.class, new Annotator[] { new AminovinylcysteineAnnotator() } ),
	TVAG(TvaGReaction.class, new Annotator[] { new HistidineAnnotator() } ),
	TVAD(TvaDReaction.class, new Annotator[] { new TvaDAnnotator() } ),
	;
	
	private Class<? extends GenericReaction> reactionClass;
	private Annotator[] annotators;
	
	private Reactions(Class<? extends GenericReaction> reaction, Annotator[] annotators) {
		this.reactionClass = reaction;
		this.annotators = annotators;
	}
	
	public Annotator[] annotators() throws ClassInstantiationException {
		return annotators;
	}
	
	public Class<? extends GenericReaction> reactionClass() {
		return reactionClass;
	}
	
}
