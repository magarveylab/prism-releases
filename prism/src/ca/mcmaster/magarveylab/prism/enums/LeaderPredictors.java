package ca.mcmaster.magarveylab.prism.enums;

import ca.mcmaster.magarveylab.enums.domains.DomainType;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.motif.LeaderPredictor;
import ca.mcmaster.magarveylab.prism.motif.leadercleavage.AutoInducingPeptidesPredictor;
import ca.mcmaster.magarveylab.prism.motif.leadercleavage.BacterialHeadToTailCyclizedPeptidePredictor;
import ca.mcmaster.magarveylab.prism.motif.leadercleavage.BottromycinsPredictor;
import ca.mcmaster.magarveylab.prism.motif.leadercleavage.ComXPredictor;
import ca.mcmaster.magarveylab.prism.motif.leadercleavage.CyanobactinsPredictor;
import ca.mcmaster.magarveylab.prism.motif.leadercleavage.GlycocinsPredictor;
import ca.mcmaster.magarveylab.prism.motif.leadercleavage.LantipeptidesPredictor;
import ca.mcmaster.magarveylab.prism.motif.leadercleavage.LassoPeptidesPredictor;
import ca.mcmaster.magarveylab.prism.motif.leadercleavage.LinaridinsPredictor;
import ca.mcmaster.magarveylab.prism.motif.leadercleavage.LinearAzoleContainingPeptidesPredictor;
import ca.mcmaster.magarveylab.prism.motif.leadercleavage.MicroviridinsPredictor;
import ca.mcmaster.magarveylab.prism.motif.leadercleavage.ProchlorosinsPredictor;
import ca.mcmaster.magarveylab.prism.motif.leadercleavage.ProteusinsPredictor;
import ca.mcmaster.magarveylab.prism.motif.leadercleavage.SactipeptidesPredictor;
import ca.mcmaster.magarveylab.prism.motif.leadercleavage.StreptidesPredictor;
import ca.mcmaster.magarveylab.prism.motif.leadercleavage.ThiopeptidesPredictor;
import ca.mcmaster.magarveylab.prism.motif.leadercleavage.ThioviridamidePredictor;
import ca.mcmaster.magarveylab.prism.motif.leadercleavage.TrifolitoxinsPredictor;
import ca.mcmaster.magarveylab.prism.motif.leadercleavage.YM216391Predictor;

/**
 * Associates leader peptide cleavage predictors with the domain types they
 * predict cleavage patterns for.
 * 
 * @author skinnider
 *
 */
public enum LeaderPredictors {
	
	AUTO_INDUCING_PEPTIDE(RibosomalDomains.AgrD, new AutoInducingPeptidesPredictor()),
	BOTTROMYCIN(RibosomalDomains.BotA, new BottromycinsPredictor()),
	COMX(RibosomalDomains.ComX, new ComXPredictor()),
	CYANOBACTIN(RibosomalDomains.PatE, new CyanobactinsPredictor()),
	GLYCOCINS(RibosomalDomains.SunA, new GlycocinsPredictor()),
	HEAD_TO_TAIL_CYCLIZED(RibosomalDomains.Head_to_tail_precursor, new BacterialHeadToTailCyclizedPeptidePredictor()),
	LANTIPEPTIDE_A(RibosomalDomains.LanA_1, new LantipeptidesPredictor()),
	LANTIPEPTIDE_B(RibosomalDomains.LanA_2, new LantipeptidesPredictor()),
	LANTIPEPTIDE_C(RibosomalDomains.LanA_3, new LantipeptidesPredictor()),
	LANTIPEPTIDE_D1(RibosomalDomains.LanA_4A, new LantipeptidesPredictor()),
	LANTIPEPTIDE_D2(RibosomalDomains.LanA_4B, new LantipeptidesPredictor()),
	LANTIPEPTIDE_E(RibosomalDomains.LanA_5, new LantipeptidesPredictor()),
	LANTIPEPTIDE_F(RibosomalDomains.LanA_6, new LantipeptidesPredictor()),
	LANTIPEPTIDE_G1(RibosomalDomains.LanA_7A, new LantipeptidesPredictor()),
	LANTIPEPTIDE_G2(RibosomalDomains.LanA_7B, new LantipeptidesPredictor()),
	LANTIPEPTIDE_H1(RibosomalDomains.LanA_8A, new LantipeptidesPredictor()),
	LANTIPEPTIDE_H2(RibosomalDomains.LanA_8B, new LantipeptidesPredictor()),
	LANTIPEPTIDE_I1(RibosomalDomains.LanA_9A, new LantipeptidesPredictor()),
	LANTIPEPTIDE_I2(RibosomalDomains.LanA_9B, new LantipeptidesPredictor()),
	PUTATIVE_LANTIPEPTIDE(RibosomalDomains.Putative_lantipeptide, new LantipeptidesPredictor()),
	LASSO_PEPTIDE(RibosomalDomains.Lasso_precursors, new LassoPeptidesPredictor()),
	PUTATIVE_LASSO_PEPTIDE(RibosomalDomains.Putative_lasso_peptide, new LassoPeptidesPredictor()),
	LINARDIN(RibosomalDomains.CypA, new LinaridinsPredictor()),
	LINARDIN_2(RibosomalDomains.LegA, new LinaridinsPredictor()),
	LAP_1(RibosomalDomains.LAP_1, new LinearAzoleContainingPeptidesPredictor()),
	LAP_2(RibosomalDomains.LAP_2, new LinearAzoleContainingPeptidesPredictor()),
	LAP_3(RibosomalDomains.LAP_3, new LinearAzoleContainingPeptidesPredictor()),
	LAP_4(RibosomalDomains.LAP_4, new LinearAzoleContainingPeptidesPredictor()),
	LAP_5(RibosomalDomains.LAP_5, new LinearAzoleContainingPeptidesPredictor()),
	MICROVIRIDIN(RibosomalDomains.MdnA, new MicroviridinsPredictor()),
	PROCHLOROSIN(RibosomalDomains.ProcA, new ProchlorosinsPredictor()),
	PROTEUSIN(RibosomalDomains.PoyA, new ProteusinsPredictor()),
	SACTIPEPTIDE(RibosomalDomains.SboA, new SactipeptidesPredictor()),
	STREPTIDE(RibosomalDomains.StrA, new StreptidesPredictor()),
	TRIFOLITOXINS(RibosomalDomains.TfxA, new TrifolitoxinsPredictor()),
	THIOPEPTIDE(RibosomalDomains.LazA, new ThiopeptidesPredictor()),
	THIOVIRIDAMIDE(RibosomalDomains.TvaA, new ThioviridamidePredictor()),
	YM(RibosomalDomains.YmA, new YM216391Predictor()),
	;
	
	private DomainType domain;
	private LeaderPredictor predictor;
	
	private LeaderPredictors(DomainType domain, LeaderPredictor predictor) {
		this.domain = domain;
		this.predictor = predictor;
	}

	/**
	 * Get the type of domain associated with this leader cleavage site(s)
	 * predictor.
	 * 
	 * @return the domain associated with this leader cleavage site(s)
	 *         predictor.
	 */
	public DomainType domain() {
		return domain;
	}

	/**
	 * Get the leader cleavage site(s) predictor for this class of ribosomal
	 * compounds.
	 * 
	 * @return the leader cleavage site(s) predictor for this class of ribosomal
	 *         compounds
	 */
	public LeaderPredictor predictor() {
		return predictor;
	}

}
