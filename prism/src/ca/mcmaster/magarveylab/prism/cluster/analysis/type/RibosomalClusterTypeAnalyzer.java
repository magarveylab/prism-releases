package ca.mcmaster.magarveylab.prism.cluster.analysis.type;

import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.enums.clusters.ClusterType;
import ca.mcmaster.magarveylab.enums.clusters.RibosomalClusterTypes;
import ca.mcmaster.magarveylab.enums.domains.RibosomalDomains;
import ca.mcmaster.magarveylab.prism.data.Cluster;

/**
 * Determine whether clusters in question are ribosomal clusters, and determine
 * their subtypes.
 * 
 * @author skinnider
 *
 */
public class RibosomalClusterTypeAnalyzer {

	/**
	 * Get the subtype(s) of this ribosomal cluster.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return ribosomal cluster subtype(s)
	 */
	public static List<ClusterType> getTypes(Cluster cluster) {
		List<ClusterType> types = new ArrayList<ClusterType>();
		if (isAutoinducingPeptideCluster(cluster))
			types.add(RibosomalClusterTypes.AUTO_INDUCING_PEPTIDE);
		if (isBacterialHeadToTailCyclizedCluster(cluster))
			types.add(RibosomalClusterTypes.BACTERIAL_HEAD_TO_TAIL_CYCLIZED);
		if (isBottromycinCluster(cluster))
			types.add(RibosomalClusterTypes.BOTTROMYCIN);
		if (isComXCluster(cluster))
			types.add(RibosomalClusterTypes.COMX);
		if (isCyanobactinCluster(cluster))
			types.add(RibosomalClusterTypes.CYANOBACTIN);
		if (isGlyocinCluster(cluster))
			types.add(RibosomalClusterTypes.GLYCOCIN);
		if (isClassILantipeptideCluster(cluster))
			types.add(RibosomalClusterTypes.CLASS_I_LANTIPEPTIDE);
		if (isClassIILantipeptideCluster(cluster))
			types.add(RibosomalClusterTypes.CLASS_II_LANTIPEPTIDE);
		if (isClassIIIorIVLantipeptideCluster(cluster))
			types.add(RibosomalClusterTypes.CLASS_III_IV_LANTIPEPTIDE);
		if (isProchlorosinLantipeptideCluster(cluster))
			types.add(RibosomalClusterTypes.PROCHLOROSIN);
		if (isLassoPeptideCluster(cluster))
			types.add(RibosomalClusterTypes.LASSO_PEPTIDE);
		if (isLinardinCluster(cluster))
			types.add(RibosomalClusterTypes.LINARIDIN);
		if (isLinearAzoleCluster(cluster))
			types.add(RibosomalClusterTypes.LINEAR_AZOLE_CONTAINING_PEPTIDE);
		if (isMicroviridinCluster(cluster))
			types.add(RibosomalClusterTypes.MICROVIRIDIN);
		if (isProteusinCluster(cluster))
			types.add(RibosomalClusterTypes.PROTEUSIN);
		if (isSactipeptideCluster(cluster))
			types.add(RibosomalClusterTypes.SACTIPEPTIDE);
		if (isStreptideCluster(cluster))
			types.add(RibosomalClusterTypes.STREPTIDE);
		if (isThiopeptideCluster(cluster))
			types.add(RibosomalClusterTypes.THIOPEPTIDE);
		if (isTrifolitoxinCluster(cluster))
			types.add(RibosomalClusterTypes.TRIFOLITOXIN);
		if (isThioviridamideCluster(cluster))
			types.add(RibosomalClusterTypes.THIOVIRIDAMIDE);
		if (isYM216391Cluster(cluster))
			types.add(RibosomalClusterTypes.YM216391);
		return types;
	}

	/**
	 * Determine whether this cluster is a ribosomal cluster of the autoinducing
	 * peptide type. Requires the precursor (AgrD) and the endopeptidase (AgrB).
	 * <br><br>
	 * Auto-inducing peptides (AIPs) are macrocyclized peptides used for
	 * Staphylococcus (and Clostridia) quorum sensing. AgrB cleaves the
	 * precursor and apparently catalyzes the intramolecular ester or thioester
	 * formation between the C-terminus and a cysteine or serine, generally
	 * located four amino acids from the C-terminus (SKACFMFV).
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if autoinducing peptide-type cluster
	 */
	public static boolean isAutoinducingPeptideCluster(Cluster cluster) {
		return (cluster.contains(RibosomalDomains.AgrD)
				&& cluster.contains(RibosomalDomains.AgrB));
	}

	/**
	 * Determine whether this cluster is a ribosomal cluster of the bacterial
	 * head-to-tail cyclized type. Requires the precursor peptide and DUF95
	 * family protein (for cyclization, seemingly).
	 * <br><br>
	 * These ribosomal peptides have a highly variable leader peptide (between 2
	 * and 30+ residues at the N-terminus) which is removed, followed by
	 * head-to-tail macrolactamization, theoretically facilitated by the thus
	 * uncharacterized DUF95 family enzyme.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if bacterial head-to-tail cyclized-type cluster
	 */
	public static boolean isBacterialHeadToTailCyclizedCluster(Cluster cluster) {
		return (cluster.contains(RibosomalDomains.DUF95)
				&& cluster.contains(RibosomalDomains.Head_to_tail_precursor));
	}

	/**
	 * Determine whether this cluster is a ribosomal cluster of the bottromycin
	 * type. Requires the precursor (botA) and a series of specific radical SAM
	 * enzymes and Ycao's, including BotC.
	 * <br><br>
	 * At the moment, very little is known about bottromycin biosynthesis, and
	 * BLAST results indicate that only the canonical bottromycin cluster exists
	 * in a handful of Streptomyces. I would recommend simply using these two
	 * genes to define the presence or absence of bottromycin. The primary
	 * sequence tells whether the molecule is the standard bottromycin (A,B,C
	 * series) or is the Alanine variant (D).
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if bottromycin-type cluster
	 */
	public static boolean isBottromycinCluster(Cluster cluster) {
		return (cluster.contains(RibosomalDomains.BotA) && cluster
				.contains(RibosomalDomains.BotC));
	}

	/**
	 * Determine whether this cluster is a ribosomal cluster of the ComX type.
	 * Requires the ComX precursor and the isoprenyltransferase ComQ.
	 * <br><br>
	 * ComX is a small modified peptide used for Bacillus (or Gram Positive
	 * related) quorum sensing. The structure is a linear peptide with a
	 * modified tryptophan, including an intramolecular cyclization which is
	 * also prenylated. There are very few fully elucidated structures, but
	 * they're generally less than 10 residues long and contain only the one
	 * modification.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if ComX-type cluster
	 */
	public static boolean isComXCluster(Cluster cluster) {
		return (cluster.contains(RibosomalDomains.ComQ) && cluster
				.contains(RibosomalDomains.ComX));
	}

	/**
	 * Determine whether this cluster is a ribosomal cluster of the cyanobactin
	 * type. Requires the precursor (patE) and two proteases (patA and patG).
	 * <br><br>
	 * Cyanobactins are a fairly poorly understood class of cyclic ribosomal
	 * peptides. Multiple cyanobactins may be produced from a single precursor
	 * peptide, flanked by N- and C-terminal cleavage sites. The C-terminal
	 * cyclase and protease patG often also has an oxidative domain that can
	 * convert thiazolines and oxazolines into thiazoles and oxazoles (adding
	 * another double bond in the heterocycle). Thiazolines and oxazolines are
	 * formed by a dehydratase domain (patD). Cyanobactins can also be linear if
	 * a hybrid prenyltransferase/methyltransferase is present, which generates
	 * an N-prenylated and terminal COOH O-methylated structure.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if cyanobactin-type cluster
	 */
	public static boolean isCyanobactinCluster(Cluster cluster) {
		return (cluster.contains(RibosomalDomains.PatA) 
				&& (cluster.contains(RibosomalDomains.PatG) 
				|| cluster.contains(RibosomalDomains.PatG_ox)));
	}

	/**
	 * Determine whether this cluster is a ribosomal cluster of the glyocin
	 * type. Requires a precursor (sunA) as well as the S-glycosyltransferase
	 * SunS.
	 * <br><br>
	 * Glycocins are extremely poorly-studied family of ribosomal peptides,
	 * classified by a hairpin structure held together by two disulfide bonds,
	 * along with an S-glycosylation. Typically the turn or the hairpin is
	 * glycosylated at either a Cys or Ser. If the hairpin is O-glycosylated,
	 * the tail (which should be a cysteine) will be S-glycosylated.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if glyocin-type cluster
	 */
	public static boolean isGlyocinCluster(Cluster cluster) {
		return (cluster.contains(RibosomalDomains.SunA) 
				&& cluster.contains(RibosomalDomains.SunS));
	}

	/**
	 * Determine whether this cluster is a class I, II, or III lantipeptide.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if this cluster meets criteria for class I, II, or III
	 *         lantipeptides
	 */
	public static boolean isLantipeptideCluster(Cluster cluster) {
		return isClassILantipeptideCluster(cluster)
				|| isClassIILantipeptideCluster(cluster)
				|| isClassIIIorIVLantipeptideCluster(cluster);
	}

	/**
	 * Determine whether this cluster is a class I lantipeptide. Requires a LanB
	 * (dehydratase) and LanC (cyclase).<br>
	 * <br>
	 * Lantibiotics are a very diverse class of ribosomal peptides that can be
	 * characterized by their lantithione bonds between serines and threonines
	 * converted into Dha and Dhb that have been attacked by cysteines. In class
	 * 1 lantibiotics, LanB converts Ser/Thr to Dha/Dhb before LanC catalyzes
	 * lantithione bond formation.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if class I lantipeptide cluster
	 */
	public static boolean isClassILantipeptideCluster(Cluster cluster) {
		return cluster.contains(RibosomalDomains.LanB)
				&& cluster.contains(RibosomalDomains.LanC);
	}

	/**
	 * Determine whether this cluster is a class II lantipeptide. Requires a
	 * LanM (fused dehydratase/cyclase).<br>
	 * <br>
	 * In class 2 lantibiotics, LanM catalyzes Dha/Dhb formation and lantithione
	 * bond formation.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if class II lantipeptide cluster
	 */
	public static boolean isClassIILantipeptideCluster(Cluster cluster) {
		return cluster.contains(RibosomalDomains.LanM);
	}

	/**
	 * Determine whether this cluster is a class III/IV lantipeptide. Requires a
	 * LanKC/LanL (Lyase/Kinase/Cyclase). LanL is distinguished from LanKC by
	 * point mutations.<br>
	 * <br>
	 * In class 3 lantibiotics, the LanKC enzyme catalyzes both reactions as
	 * well as labionin bond formation. Labionin bonds are lantithione bonds
	 * with a second Dhb/Dha R-group attacking the lantithione Dha/Dhb residue
	 * at the alpha-carbon.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if class III lantipeptide cluster
	 */
	public static boolean isClassIIIorIVLantipeptideCluster(Cluster cluster) {
		return cluster.contains(RibosomalDomains.LanKC);
	}
	
	/**
	 * Prochlorosins are a unique class of type II lantipeptides produced by
	 * cyanobacteria, in which a single LanM-type enzyme encoded anywhere within
	 * the genome catalyzes lantithione bond formation in trans within a diverse
	 * range of precursor peptides. A putative prochlorosin cluster therefore
	 * only requires a homolog of the ProcA precursor peptide.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if prochlorosin cluster
	 */
	public static boolean isProchlorosinLantipeptideCluster(Cluster cluster) {
		return (cluster.contains(RibosomalDomains.ProcA));
	}

	/**
	 * Determine whether this cluster is a ribosomal cluster of the lassopeptide
	 * type. Requires a transglutaminase and asparagine synthase. <br>
	 * <br>
	 * Lasso peptides are thoroughly predictable, and are not known to be
	 * modified. After cleavage of the leader peptide by a conserved
	 * transglutaminase-like enzyme, the amine of the N-terminal residue
	 * (typically glycine) is used to attack a aspartate or glutamate located 7
	 * or 8 residues downstream, forming a macrolactam loop with the tail
	 * passing through. Disulfide bonds occasionally form to hold the tail in
	 * place.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if lasso peptide cluster
	 */
	public static boolean isLassoPeptideCluster(Cluster cluster) {
		return (cluster.contains(RibosomalDomains.Transglutaminase)
				&& cluster.contains(RibosomalDomains.Asparagine_synthase));
	}

	/**
	 * Determine whether this cluster is a ribosomal cluster of the linardin
	 * type. Requires: cypA/legA (precursor), cypH/legH, and cypL.
	 * <br><br>
	 * Linaridins are a small class of ribosomal peptides that are very closely
	 * related to lantibiotics and can be distinguished by the lack of
	 * lantithione bonds. Cypemycin is the prototypical linaridin, possessing a
	 * number of Dha and Dhb residues, along with an AviCys macrocycle tail
	 * (also seen in lantibiotics like epidermin / gallidermin). After leader
	 * peptide cleavage the lantibiotic dehydratases CypH and CypL convert
	 * threonines and serines into Dhb and Dha. The C-terminal cysteine is
	 * decarboxylated before the thiol cycles back to form the AviCys moiety.
	 * Other modifying enzymes then step in the finish the molecule.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if linardin-type cluster
	 */
	public static boolean isLinardinCluster(Cluster cluster) {
		boolean precursor = cluster.contains(RibosomalDomains.CypA)
				|| cluster.contains(RibosomalDomains.LegA);
		boolean cypH = cluster.contains(RibosomalDomains.CypH)
				|| cluster.contains(RibosomalDomains.LegH);
		return (precursor && cypH && cluster.contains(RibosomalDomains.CypL));
	}

	/**
	 * Determine whether this cluster is a ribosomal cluster of the linear
	 * azole-containing peptide type. Requires the precursor peptide (A), a
	 * dehydrogenase (B), and a potentially fused cyclodehydratase (C/D).
	 * <br><br>
	 * LAPs possess an N-terminal leader peptide that must be cleaved prior to
	 * molecule construction. All cysteines, serines, and threonines are
	 * converted to thiazolines and oxazolines respectively by the C and D
	 * dehyratases, then to thiazoles and oxazoles by the dehydrogenase (B).
	 * There are only three fully characterized examples, and except for the
	 * large number of thiazoles and oxazoles, other modifications are sparse,
	 * including Dha formation (godG and H), N-terminal acetylation (godH) and
	 * N-terminal dimethylation (pznL). Dha formation is basically identical to
	 * the system used in thiazolyl peptide biosynthesis.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if linear azole-containing peptide-type cluster
	 */
	public static boolean isLinearAzoleCluster(Cluster cluster) {
		return (cluster.contains(RibosomalDomains.McbB) 
				&& (cluster.contains(RibosomalDomains.McbC) 
				|| cluster.contains(RibosomalDomains.McbD)) 
				|| cluster.contains(RibosomalDomains.GodG));
	}

	/**
	 * Determine whether this cluster is a ribosomal cluster of the microviridin
	 * type. Requires MdnA and one or two ATP-grasp enzymes (MdnB and/or MdnC).
	 * <br><br>
	 * After cleaving the leader peptide, the ATP grasp enzymes MdnC and/or B
	 * will form intramolecular ester and amide bonds (respectively). The two
	 * known molecule types include 3 (microviridins) or 2 (marinostatins)
	 * crosslinks, with marinostatin possessing only two ester bonds, and thus,
	 * only one ATP grasp enzyme.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if microviridin-type cluster
	 */
	public static boolean isMicroviridinCluster(Cluster cluster) {
		return (cluster.contains(RibosomalDomains.MdnA) 
				&& (cluster.contains(RibosomalDomains.MdnB) 
				|| cluster.contains(RibosomalDomains.MdnC)));
	}

	/**
	 * Determine whether this cluster is a ribosomal cluster of the proteusin
	 * type. Requires poyA (precursor) and poyD (epimerase). Usually found with
	 * at least one other radical SAM
	 * <br><br>
	 * Proteusins are a very unusual class of ribosomal peptides with one
	 * partially characterized gene cluster and three preliminarily described
	 * molecules. The leader peptide each molecule is quite large and typically
	 * resembles a protein (often nitrile hydratase). At the moment the
	 * precursor and the poyD epimerase seem to be fairly diagnostic, but that
	 * may change with time. These molecules are essentially unpredictable since
	 * the modifying enzymes appear to work several times with no obvious
	 * 'code', instead relying on the secondary structure of the small molecule
	 * to place residues in the correct position for modification. Modifying
	 * enzymes found the described gene clusters include Asparagine
	 * N-methyltransferases, beta-carbon hydroxylating enzymes, and two closely
	 * related but apparently distinct C-methyltransferases (radical SAMs).
	 * Predictability is limited to leader peptide cleavage, which seems to
	 * follow the general bacteriocin GG protease signal.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if proteusin-type cluster
	 */
	public static boolean isProteusinCluster(Cluster cluster) {
		return (cluster.contains(RibosomalDomains.PoyA) 
				&& cluster.contains(RibosomalDomains.PoyD));
	}

	/**
	 * Determine whether this cluster is a ribosomal cluster of the sactipeptide
	 * type. Requires precursor (SboA) and Fe-S protein (AlbA).
	 * <br><br>
	 * Sactipeptides are a relatively simple class of ribosomal peptides that
	 * form hairpins held together by S-C bonds between cysteine R-groups and
	 * alpha-carbons. After leader peptide cleavage, these bonds are installed
	 * by the distinct sactipeptide enzyme AlbA.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if sactipeptide-type cluster
	 */
	public static boolean isSactipeptideCluster(Cluster cluster) {
		return (cluster.contains(RibosomalDomains.SboA) 
				&& cluster.contains(RibosomalDomains.AlbA));
	}
	
	/**
	 * Determine whether this cluster is a ribosomal cluster of the streptide
	 * type. Requires precursor (StrA), SPASM domain-containing radical SAM
	 * enzyme (StrB), and transporter (StrC).<br>
	 * <br>
	 * Streptide is a quorum-sensing peptide from Streptococcus and Lactococcus
	 * (though not described yet) with a single described structure. The
	 * prepeptide (StrA) is is cut on the N-terminal side, and, in the one known
	 * structure, on the C-terminal side as well. Streptide possesses a unique
	 * modification: a C-C bond between a conserved lysine and a conserved
	 * tryptophan, catalyzed by a conserved Fe-S cluster protein (StrB). The
	 * final small molecule is exported by a conserved pump (StrC).
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if streptide cluster
	 */
	public static boolean isStreptideCluster(Cluster cluster) {
		return (cluster.contains(RibosomalDomains.StrA) 
				&& cluster.contains(RibosomalDomains.StrB) 
				&& cluster.contains(RibosomalDomains.StrC));
	}

	/**
	 * Determine whether this cluster is a ribosomal cluster of the thiopeptide
	 * type. Requires a LazB and LazC. 
	 * <br><br>
	 * Thiopeptides (or thiazolyl peptides) are far and away the most
	 * complicated ribosomal natural products. In their most basic form, these
	 * molecules are azole containing peptides with a unique 4+2 cycloaddition
	 * between two Dha residues that forms their trademark macrocycle and
	 * pyrizine core.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if thiopeptide-type cluster
	 */
	public static boolean isThiopeptideCluster(Cluster cluster) {
		return (cluster.contains(RibosomalDomains.LazB) 
				&& cluster.contains(RibosomalDomains.LazC));
	}

	/**
	 * Determine whether this cluster is a ribosomal cluster of the trifolitoxin
	 * type. Requires tfxA (precursor), tfxB (nitroreductase and sagB), tfxC
	 * (nitroreductase and FMN)
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if thiopeptide-type cluster
	 */
	public static boolean isTrifolitoxinCluster(Cluster cluster) {
		return (cluster.contains(RibosomalDomains.TfxA)
				&& cluster.contains(RibosomalDomains.TfxB) 
				&& cluster.contains(RibosomalDomains.TfxC));
	}

	/**
	 * Determine whether this cluster is a ribosomal cluster of the
	 * thioviridamide type. Requires the TvaA precursor and the
	 * thioamide-forming protein TvaH. <br>
	 * <br>
	 * Thioviridamide is the only described molecule from an apparent class of
	 * thioamide containing ribosomal peptides. Given that its biosynthesis is
	 * largely based on hypothetical protein functions, I'd hesitate to say what
	 * any other thioamide containing molecule may look like, so at the moment
	 * we'll only describe thioviridamide and related structures / gene
	 * clusters. <br>
	 * <br>
	 * As a conserved first step, the thioviridamide leader should be removed,
	 * followed by C-terminal cysteine decarboxylation (tvaF) and cyclization,
	 * creating the C-terminal avicys macrocycle. All amide ketones in the
	 * remaining linear portion of the molecule are converted to thioamides by
	 * the unique tvaH enzyme. The penultimate histidine is then
	 * beta-hydroxylated (tvaJ) and methylated on both side chain nitrogens
	 * (tvaG). <br>
	 * <br>
	 * A remaining mystery is the conversion of the N-terminal serine to the
	 * unusual monomer shown below. Later studies will hopefully clarify exactly
	 * which enzymes are responsible.
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if thioviridamide-type cluster
	 */
	public static boolean isThioviridamideCluster(Cluster cluster) {
		return cluster.contains(RibosomalDomains.TvaA)
				&& cluster.contains(RibosomalDomains.TvaH);
	}
	
	/**
	 * Determine whether this cluster is a ribosomal cluster of the YM-216391
	 * type. Requires the YmA precursor and unique macrocyclase, which
	 * distinguishes this cluster from the cyanobactins (see van der Donk
	 * review). <br>
	 * <br>
	 * This ribosomal peptide is one of a family of three molecules, which are
	 * macrocyclized thiazole / oxazole containing compounds that incorporate a
	 * beta-hydroxylated phenylalanine into one of the oxazoles. <br>
	 * <br>
	 * Like other thiazole and oxazole containing molecules, the leader (and
	 * follower, in this case) peptide is removed, and heterocycles are
	 * installed by a cyclodehydratase (ymD [not listed here] and ymBC
	 * N-terminal domain), and then oxidized to azoles (by the ymBC C-terminal
	 * domain). Macrocyclization likely follows the formation of the azoles. In
	 * this unique case, a phenylalanine is hydroxylated at the beta position
	 * (ymE), and heterocyclized on an adjacent amide ketone (ymB1), and then
	 * oxidized to an oxazole (ymC1).
	 * 
	 * @param cluster
	 *            cluster to analyze
	 * @return true if YM-216391-type cluster
	 */
	public static boolean isYM216391Cluster(Cluster cluster) {
		return cluster.contains(RibosomalDomains.YmA)
				&& cluster.contains(RibosomalDomains.YmF);
	}

}
