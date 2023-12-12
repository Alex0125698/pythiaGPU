// Settings.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the settings database.

#ifndef Pythia8_Settings_H
#define Pythia8_Settings_H

// #include "Pythia8/Info.h"
#include "Pythia8/PythiaStdlib.h"

namespace Pythia8 {

// setting names

enum class Flag
{
  Alpgen_setHeavyMasses,
  Alpgen_setLightMasses,
  Alpgen_setMLM,
  Alpgen_setNjet,
  BeamRemnants_allowBeamJunction,
  BeamRemnants_allowJunction,
  BeamRemnants_beamJunction,
  BeamRemnants_primordialKT,
  BeamRemnants_rescatterRestoreY,
  Beams_allowMomentumSpread,
  Beams_allowVertexSpread,
  Beams_newLHEFsameInit,
  Beams_readLHEFheaders,
  Beams_setProductionScalesFromLHEF,
  Beams_strictLHEFscale,
  BoseEinstein_Eta,
  BoseEinstein_Kaon,
  BoseEinstein_Pion,
  Bottomonium_all,
  Charmonium_all,
  Check_abortIfVeto,
  Check_event,
  Check_history,
  Check_particleData,
  ColourReconnection_allowDoubleJunRem,
  ColourReconnection_allowJunctions,
  ColourReconnection_forceHadronLevelCR,
  ColourReconnection_forceResonance,
  ColourReconnection_lowerLambdaOnly,
  ColourReconnection_reconnect,
  ColourReconnection_sameNeighbourColours,
  ColourReconnection_singleReconnection,
  ContactInteractions_QCffbar2eebar,
  ContactInteractions_QCffbar2mumubar,
  ContactInteractions_QCffbar2tautaubar,
  ContactInteractions_QCqq2qq,
  ContactInteractions_QCqqbar2qqbar,
  Diffraction_doHard,
  ExcitedFermion_all,
  ExcitedFermion_bg2bStar,
  ExcitedFermion_cg2cStar,
  ExcitedFermion_dg2dStar,
  ExcitedFermion_egm2eStar,
  ExcitedFermion_mugm2muStar,
  ExcitedFermion_qq2bStarq,
  ExcitedFermion_qq2cStarq,
  ExcitedFermion_qq2dStarq,
  ExcitedFermion_qq2sStarq,
  ExcitedFermion_qq2uStarq,
  ExcitedFermion_qqbar2eStare,
  ExcitedFermion_qqbar2eStareStar,
  ExcitedFermion_qqbar2muStarmu,
  ExcitedFermion_qqbar2muStarmuStar,
  ExcitedFermion_qqbar2nueStarnue,
  ExcitedFermion_qqbar2nueStarnueStar,
  ExcitedFermion_qqbar2numuStarnumu,
  ExcitedFermion_qqbar2numuStarnumuStar,
  ExcitedFermion_qqbar2nutauStarnutau,
  ExcitedFermion_qqbar2nutauStarnutauStar,
  ExcitedFermion_qqbar2tauStartau,
  ExcitedFermion_qqbar2tauStartauStar,
  ExcitedFermion_sg2sStar,
  ExcitedFermion_taugm2tauStar,
  ExcitedFermion_ug2uStar,
  ExtraDimensionsGstar_all,
  ExtraDimensionsGstar_ffbar2Gstar,
  ExtraDimensionsGstar_gg2Gstar,
  ExtraDimensionsGstar_gg2Gstarg,
  ExtraDimensionsGstar_qg2Gstarq,
  ExtraDimensionsGstar_qqbar2Gstarg,
  ExtraDimensionsGstar_qqbar2KKgluonstar,
  ExtraDimensionsGstar_SMinBulk,
  ExtraDimensionsGstar_VLVL,
  ExtraDimensionsLED_dijets,
  ExtraDimensionsLED_ffbar2gammagamma,
  ExtraDimensionsLED_ffbar2Ggamma,
  ExtraDimensionsLED_ffbar2GZ,
  ExtraDimensionsLED_ffbar2llbar,
  ExtraDimensionsLED_gg2DJgg,
  ExtraDimensionsLED_gg2DJqqbar,
  ExtraDimensionsLED_gg2gammagamma,
  ExtraDimensionsLED_gg2Gg,
  ExtraDimensionsLED_gg2llbar,
  ExtraDimensionsLED_GravScalar,
  ExtraDimensionsLED_monojet,
  ExtraDimensionsLED_qg2DJqg,
  ExtraDimensionsLED_qg2Gq,
  ExtraDimensionsLED_qq2DJqq,
  ExtraDimensionsLED_qqbar2DJgg,
  ExtraDimensionsLED_qqbar2DJqqbarNew,
  ExtraDimensionsLED_qqbar2Gg,
  ExtraDimensionsTEV_ffbar2bbbar,
  ExtraDimensionsTEV_ffbar2ccbar,
  ExtraDimensionsTEV_ffbar2ddbar,
  ExtraDimensionsTEV_ffbar2epluseminus,
  ExtraDimensionsTEV_ffbar2muplusmuminus,
  ExtraDimensionsTEV_ffbar2nuenuebar,
  ExtraDimensionsTEV_ffbar2numunumubar,
  ExtraDimensionsTEV_ffbar2nutaunutaubar,
  ExtraDimensionsTEV_ffbar2ssbar,
  ExtraDimensionsTEV_ffbar2tauplustauminus,
  ExtraDimensionsTEV_ffbar2ttbar,
  ExtraDimensionsTEV_ffbar2uubar,
  ExtraDimensionsUnpart_ffbar2gammagamma,
  ExtraDimensionsUnpart_ffbar2llbar,
  ExtraDimensionsUnpart_ffbar2Ugamma,
  ExtraDimensionsUnpart_ffbar2UZ,
  ExtraDimensionsUnpart_gg2gammagamma,
  ExtraDimensionsUnpart_gg2llbar,
  ExtraDimensionsUnpart_gg2Ug,
  ExtraDimensionsUnpart_monojet,
  ExtraDimensionsUnpart_qg2Uq,
  ExtraDimensionsUnpart_qqbar2Ug,
  FourthBottom_all,
  FourthBottom_ffbar2bPrimebPrimebar_s_gmZ_,
  FourthBottom_ffbar2bPrimeqbar_s_W_,
  FourthBottom_ffbar2bPrimetbar_s_W_,
  FourthBottom_gg2bPrimebPrimebar,
  FourthBottom_qq2bPrimeq_t_W_,
  FourthBottom_qqbar2bPrimebPrimebar,
  FourthPair_ffbar2tauPrimenuPrimebar_s_W_,
  FourthPair_ffbar2tPrimebPrimebar_s_W_,
  FourthTop_all,
  FourthTop_ffbar2tPrimeqbar_s_W_,
  FourthTop_ffbar2tPrimetPrimebar_s_gmZ_,
  FourthTop_gg2tPrimetPrimebar,
  FourthTop_qq2tPrimeq_t_W_,
  FourthTop_qqbar2tPrimetPrimebar,
  HadronLevel_all,
  HadronLevel_BoseEinstein,
  HadronLevel_Decay,
  HadronLevel_Hadronize,
  HadronScatter_afterDecay,
  HadronScatter_allowDecayProd,
  HadronScatter_scatter,
  HadronScatter_scatterRepeat,
  HadronScatter_tile,
  HardQCD_3parton,
  HardQCD_all,
  HardQCD_gg2bbbar,
  HardQCD_gg2ccbar,
  HardQCD_gg2gg,
  HardQCD_gg2ggg,
  HardQCD_gg2qqbar,
  HardQCD_gg2qqbarg,
  HardQCD_hardbbbar,
  HardQCD_hardccbar,
  HardQCD_qg2qg,
  HardQCD_qg2qgg,
  HardQCD_qg2qqqbarDiff,
  HardQCD_qg2qqqbarSame,
  HardQCD_qq2qq,
  HardQCD_qq2qqgDiff,
  HardQCD_qq2qqgSame,
  HardQCD_qqbar2bbbar,
  HardQCD_qqbar2ccbar,
  HardQCD_qqbar2gg,
  HardQCD_qqbar2ggg,
  HardQCD_qqbar2qqbargDiff,
  HardQCD_qqbar2qqbargSame,
  HardQCD_qqbar2qqbarNew,
  HiddenValley_all,
  HiddenValley_doKinMix,
  HiddenValley_ffbar2BvBvbar,
  HiddenValley_ffbar2CvCvbar,
  HiddenValley_ffbar2DvDvbar,
  HiddenValley_ffbar2EvEvbar,
  HiddenValley_ffbar2MUvMUvbar,
  HiddenValley_ffbar2nuEvnuEvbar,
  HiddenValley_ffbar2nuMUvnuMUvbar,
  HiddenValley_ffbar2nuTAUvnuTAUvbar,
  HiddenValley_ffbar2SvSvbar,
  HiddenValley_ffbar2TAUvTAUvbar,
  HiddenValley_ffbar2TvTvbar,
  HiddenValley_ffbar2UvUvbar,
  HiddenValley_ffbar2Zv,
  HiddenValley_fragment,
  HiddenValley_FSR,
  HiddenValley_gg2BvBvbar,
  HiddenValley_gg2CvCvbar,
  HiddenValley_gg2DvDvbar,
  HiddenValley_gg2SvSvbar,
  HiddenValley_gg2TvTvbar,
  HiddenValley_gg2UvUvbar,
  HiddenValley_qqbar2BvBvbar,
  HiddenValley_qqbar2CvCvbar,
  HiddenValley_qqbar2DvDvbar,
  HiddenValley_qqbar2SvSvbar,
  HiddenValley_qqbar2TvTvbar,
  HiddenValley_qqbar2UvUvbar,
  Higgs_clipWings,
  Higgs_cubicWidth,
  Higgs_runningLoopMass,
  Higgs_useBSM,
  HiggsBSM_all,
  HiggsBSM_allA3,
  HiggsBSM_allHplusminus,
  HiggsBSM_allH1,
  HiggsBSM_allH2,
  HiggsBSM_allHpair,
  HiggsBSM_bg2Hplusminust,
  HiggsBSM_ff2A3ff_t_WW_,
  HiggsBSM_ff2A3ff_t_ZZ_,
  HiggsBSM_ff2H1ff_t_WW_,
  HiggsBSM_ff2H1ff_t_ZZ_,
  HiggsBSM_ff2H2ff_t_WW_,
  HiggsBSM_ff2H2ff_t_ZZ_,
  HiggsBSM_ffbar2A3,
  HiggsBSM_ffbar2A3H1,
  HiggsBSM_ffbar2A3H2,
  HiggsBSM_ffbar2A3W,
  HiggsBSM_ffbar2A3Z,
  HiggsBSM_ffbar2Hplusminus,
  HiggsBSM_ffbar2HplusminusH1,
  HiggsBSM_ffbar2HplusminusH2,
  HiggsBSM_ffbar2HplusHminus,
  HiggsBSM_ffbar2H1,
  HiggsBSM_ffbar2H1W,
  HiggsBSM_ffbar2H1Z,
  HiggsBSM_ffbar2H2,
  HiggsBSM_ffbar2H2W,
  HiggsBSM_ffbar2H2Z,
  HiggsBSM_gg2A3,
  HiggsBSM_gg2A3bbbar,
  HiggsBSM_gg2A3g_l_t_,
  HiggsBSM_gg2A3ttbar,
  HiggsBSM_gg2H1,
  HiggsBSM_gg2H1bbbar,
  HiggsBSM_gg2H1g_l_t_,
  HiggsBSM_gg2H1ttbar,
  HiggsBSM_gg2H2,
  HiggsBSM_gg2H2bbbar,
  HiggsBSM_gg2H2g_l_t_,
  HiggsBSM_gg2H2ttbar,
  HiggsBSM_gmgm2A3,
  HiggsBSM_gmgm2H1,
  HiggsBSM_gmgm2H2,
  HiggsBSM_qg2A3q,
  HiggsBSM_qg2A3q_l_t_,
  HiggsBSM_qg2H1q,
  HiggsBSM_qg2H1q_l_t_,
  HiggsBSM_qg2H2q,
  HiggsBSM_qg2H2q_l_t_,
  HiggsBSM_qqbar2A3bbbar,
  HiggsBSM_qqbar2A3g_l_t_,
  HiggsBSM_qqbar2A3ttbar,
  HiggsBSM_qqbar2H1bbbar,
  HiggsBSM_qqbar2H1g_l_t_,
  HiggsBSM_qqbar2H1ttbar,
  HiggsBSM_qqbar2H2bbbar,
  HiggsBSM_qqbar2H2g_l_t_,
  HiggsBSM_qqbar2H2ttbar,
  HiggsSM_all,
  HiggsSM_ff2Hff_t_WW_,
  HiggsSM_ff2Hff_t_ZZ_,
  HiggsSM_ffbar2H,
  HiggsSM_ffbar2HW,
  HiggsSM_ffbar2HZ,
  HiggsSM_gg2H,
  HiggsSM_gg2Hbbbar,
  HiggsSM_gg2Hg_l_t_,
  HiggsSM_gg2Httbar,
  HiggsSM_gmgm2H,
  HiggsSM_NLOWidths,
  HiggsSM_qg2Hq,
  HiggsSM_qg2Hq_l_t_,
  HiggsSM_qqbar2Hbbbar,
  HiggsSM_qqbar2Hg_l_t_,
  HiggsSM_qqbar2Httbar,
  Init_showAllParticleData,
  Init_showAllSettings,
  Init_showChangedParticleData,
  Init_showChangedResonanceData,
  Init_showChangedSettings,
  Init_showMultipartonInteractions,
  Init_showProcesses,
  JetMatching_doFxFx,
  JetMatching_doShowerKt,
  JetMatching_merge,
  JetMatching_setMad,
  LeftRightSymmmetry_all,
  LeftRightSymmmetry_ff2HLff,
  LeftRightSymmmetry_ff2HRff,
  LeftRightSymmmetry_ffbar2HLHL,
  LeftRightSymmmetry_ffbar2HRHR,
  LeftRightSymmmetry_ffbar2WR,
  LeftRightSymmmetry_ffbar2ZR,
  LeftRightSymmmetry_lgm2HLe,
  LeftRightSymmmetry_lgm2HLmu,
  LeftRightSymmmetry_lgm2HLtau,
  LeftRightSymmmetry_lgm2HRe,
  LeftRightSymmmetry_lgm2HRmu,
  LeftRightSymmmetry_lgm2HRtau,
  LeftRightSymmmetry_ll2HL,
  LeftRightSymmmetry_ll2HR,
  LeptoQuark_all,
  LeptoQuark_gg2LQLQbar,
  LeptoQuark_qg2LQl,
  LeptoQuark_ql2LQ,
  LeptoQuark_qqbar2LQLQbar,
  LesHouches_matchInOut,
  Main_LHEFskipInit,
  Main_spareFlag1,
  Main_spareFlag2,
  Main_spareFlag3,
  Merging_allowColourShuffling,
  Merging_allowIncompleteHistoriesInReal,
  Merging_allowSQCDClustering,
  Merging_allowWClustering,
  Merging_doCutBasedMerging,
  Merging_doKTMerging,
  Merging_doMGMerging,
  Merging_doNL3Loop,
  Merging_doNL3Subt,
  Merging_doNL3Tree,
  Merging_doPTLundMerging,
  Merging_doUMEPSSubt,
  Merging_doUMEPSTree,
  Merging_doUNLOPSLoop,
  Merging_doUNLOPSSubt,
  Merging_doUNLOPSSubtNLO,
  Merging_doUNLOPSTilde,
  Merging_doUNLOPSTree,
  Merging_doUserMerging,
  Merging_doXSectionEstimate,
  Merging_enforceCutOnLHE,
  Merging_enforceStrongOrdering,
  Merging_includeMassive,
  Merging_includeRedundant,
  Merging_mayRemoveDecayProducts,
  Merging_orderInRapidity,
  Merging_pickByFullP,
  Merging_pickByPoPT2,
  Merging_pickBySumPT,
  Merging_usePythiaQFacHard,
  Merging_usePythiaQRenHard,
  Merging_useShowerPlugin,
  MultipartonInteractions_allowDoubleRescatter,
  MultipartonInteractions_allowRescatter,
  NewGaugeBoson_ffbar2gmZZprime,
  NewGaugeBoson_ffbar2R0,
  NewGaugeBoson_ffbar2Wprime,
  Next_showMothersAndDaughters,
  Next_showScaleAndVertex,
  Onia_all,
  Onia_all_3DJ_,
  Onia_all_3PJ_,
  Onia_all_3S1_,
  Onia_forceMassSplit,
  ParticleDecays_allowPhotonRadiation,
  ParticleDecays_FSRinDecays,
  ParticleDecays_limitCylinder,
  ParticleDecays_limitRadius,
  ParticleDecays_limitTau,
  ParticleDecays_limitTau0,
  ParticleDecays_mixB,
  PartonLevel_all,
  PartonLevel_earlyResDec,
  PartonLevel_FSR,
  PartonLevel_FSRinProcess,
  PartonLevel_FSRinResonances,
  PartonLevel_ISR,
  PartonLevel_MPI,
  PartonLevel_Remnants,
  PDF_extrapolate,
  PDF_lepton,
  PDF_useHard,
  PhaseSpace_bias2Selection,
  PhaseSpace_increaseMaximum,
  PhaseSpace_sameForSecond,
  PhaseSpace_showSearch,
  PhaseSpace_showViolation,
  PhaseSpace_useBreitWigners,
  PhotonCollision_all,
  PhotonCollision_gmgm2bbbar,
  PhotonCollision_gmgm2ccbar,
  PhotonCollision_gmgm2ee,
  PhotonCollision_gmgm2mumu,
  PhotonCollision_gmgm2qqbar,
  PhotonCollision_gmgm2tautau,
  POWHEG_pythiaRandom,
  Print_quiet,
  ProcessLevel_all,
  ProcessLevel_resonanceDecays,
  PromptPhoton_all,
  PromptPhoton_ffbar2gammagamma,
  PromptPhoton_gg2gammagamma,
  PromptPhoton_gg2ggamma,
  PromptPhoton_qg2qgamma,
  PromptPhoton_qqbar2ggamma,
  Random_setSeed,
  RHadrons_allow,
  RHadrons_allowDecay,
  RHadrons_setMasses,
  SecondHard_Bottomonium,
  SecondHard_Charmonium,
  SecondHard_generate,
  SecondHard_GmZAndJet,
  SecondHard_PhotonAndJet,
  SecondHard_SingleGmZ,
  SecondHard_SingleTop,
  SecondHard_SingleW,
  SecondHard_TopPair,
  SecondHard_TwoBJets,
  SecondHard_TwoJets,
  SecondHard_TwoPhotons,
  SecondHard_WAndJet,
  SigmaDiffractive_dampen,
  SigmaElastic_setOwn,
  SigmaProcess_bMassiveME,
  SigmaProcess_cMassiveME,
  SigmaProcess_muMassiveME,
  SigmaProcess_tauMassiveME,
  SigmaTotal_setOwn,
  SigmaTotal_zeroAXB,
  SLHA_allowUserOverride,
  SLHA_keepSM,
  SLHA_NMSSM,
  SLHA_useDecayTable,
  SoftQCD_all,
  SoftQCD_centralDiffractive,
  SoftQCD_doubleDiffractive,
  SoftQCD_elastic,
  SoftQCD_inelastic,
  SoftQCD_nonDiffractive,
  SoftQCD_singleDiffractive,
  SpaceShower_alphaSuseCMW,
  SpaceShower_MEafterFirst,
  SpaceShower_MEcorrections,
  SpaceShower_phiIntAsym,
  SpaceShower_phiPolAsym,
  SpaceShower_QCDshower,
  SpaceShower_QEDshowerByL,
  SpaceShower_QEDshowerByQ,
  SpaceShower_rapidityOrder,
  SpaceShower_samePTasMPI,
  SpaceShower_useFixedFacScale,
  SpaceShower_weakShower,
  Stat_reset,
  Stat_showErrors,
  Stat_showPartonLevel,
  Stat_showProcessLevel,
  StringFlav_suppressLeadingB,
  StringZ_useNonstandardB,
  StringZ_useNonstandardC,
  StringZ_useNonstandardH,
  StringZ_usePetersonB,
  StringZ_usePetersonC,
  StringZ_usePetersonH,
  SUSY_all,
  SUSY_gg2gluinogluino,
  SUSY_gg2squarkantisquark,
  SUSY_qg2chiplusminussquark,
  SUSY_qg2chi0squark,
  SUSY_qg2squarkgluino,
  SUSY_qq2antisquark,
  SUSY_qq2squarksquark,
  SUSY_qq2squarksquark_onlyQCD,
  SUSY_qqbar2chiplusminuschi0,
  SUSY_qqbar2chiplusminusgluino,
  SUSY_qqbar2chipluschiminus,
  SUSY_qqbar2chi0chi0,
  SUSY_qqbar2chi0gluino,
  SUSY_qqbar2gluinogluino,
  SUSY_qqbar2sleptonantislepton,
  SUSY_qqbar2squarkantisquark,
  SUSY_qqbar2squarkantisquark_onlyQCD,
  SUSYResonance_3BodyMatrixElement,
  TimeShower_allowBeamRecoil,
  TimeShower_allowMPIdipole,
  TimeShower_alphaSuseCMW,
  TimeShower_dampenBeamRecoil,
  TimeShower_globalRecoil,
  TimeShower_interleave,
  TimeShower_limitPTmaxGlobal,
  TimeShower_MEafterFirst,
  TimeShower_MEcorrections,
  TimeShower_phiPolAsym,
  TimeShower_QCDshower,
  TimeShower_QEDshowerByGamma,
  TimeShower_QEDshowerByL,
  TimeShower_QEDshowerByQ,
  TimeShower_recoilToColoured,
  TimeShower_useFixedFacScale,
  TimeShower_weakShower,
  Top_all,
  Top_ffbar2tqbar_s_W_,
  Top_ffbar2ttbar_s_gmZ_,
  Top_gg2ttbar,
  Top_gmgm2ttbar,
  Top_qq2tq_t_W_,
  Top_qqbar2ttbar,
  WeakBosonAndParton_all,
  WeakBosonAndParton_ffbar2gmZgm,
  WeakBosonAndParton_ffbar2Wgm,
  WeakBosonAndParton_fgm2gmZf,
  WeakBosonAndParton_fgm2Wf,
  WeakBosonAndParton_qg2gmZq,
  WeakBosonAndParton_qg2Wq,
  WeakBosonAndParton_qqbar2gmZg,
  WeakBosonAndParton_qqbar2Wg,
  WeakBosonExchange_all,
  WeakBosonExchange_ff2ff_t_gmZ_,
  WeakBosonExchange_ff2ff_t_W_,
  WeakDoubleBoson_all,
  WeakDoubleBoson_ffbar2gmZgmZ,
  WeakDoubleBoson_ffbar2WW,
  WeakDoubleBoson_ffbar2ZW,
  WeakShower_singleEmission,
  WeakShower_vetoQCDjets,
  WeakShower_vetoWeakJets,
  WeakSingleBoson_all,
  WeakSingleBoson_ffbar2ffbar_s_gm_,
  WeakSingleBoson_ffbar2ffbar_s_gmZ_,
  WeakSingleBoson_ffbar2ffbar_s_W_,
  WeakSingleBoson_ffbar2gmZ,
  WeakSingleBoson_ffbar2W,
  Zprime_coup2gen4,
  Zprime_universality,
  END,
};

enum class Mode
{
  BeamRemnants_companionPower,
  BeamRemnants_maxValQuark,
  BeamRemnants_remnantMode,
  Beams_frameType,
  Beams_idA,
  Beams_idB,
  Beams_nSkipLHEFatInit,
  Check_levelParticleData,
  Check_nErrList,
  ColourReconnection_flipMode,
  ColourReconnection_lambdaForm,
  ColourReconnection_mode,
  ColourReconnection_nColours,
  ColourReconnection_timeDilationMode,
  ContactInteractions_etaLL,
  ContactInteractions_etaLR,
  ContactInteractions_etaRR,
  ContactInteractions_nQuarkNew,
  Diffraction_bProfile,
  Diffraction_PomFlux,
  Diffraction_sampleType,
  Event_startColTag,
  ExtraDimensionsGstar_KKintMode,
  ExtraDimensionsLED_CutOffMode,
  ExtraDimensionsLED_n,
  ExtraDimensionsLED_NegInt,
  ExtraDimensionsLED_nQuarkNew,
  ExtraDimensionsLED_opMode,
  ExtraDimensionsTEV_gmZmode,
  ExtraDimensionsTEV_nMax,
  ExtraDimensionsUnpart_CutOffMode,
  ExtraDimensionsUnpart_gXX,
  ExtraDimensionsUnpart_gXY,
  ExtraDimensionsUnpart_spinU,
  HadronScatter_hadronSelect,
  HadronScatter_scatterProb,
  HardQCD_nQuarkNew,
  HiddenValley_nFlav,
  HiddenValley_Ngauge,
  HiddenValley_spinFv,
  HiddenValley_spinqv,
  HiggsA3_parity,
  HiggsH1_parity,
  HiggsH2_parity,
  Init_showOneParticleData,
  JetMatching_exclusive,
  JetMatching_jetAlgorithm,
  JetMatching_jetAllow,
  JetMatching_jetMatch,
  JetMatching_nEta,
  JetMatching_nJet,
  JetMatching_nJetMax,
  JetMatching_nPartonsNow,
  JetMatching_nPhi,
  JetMatching_nQmatch,
  JetMatching_scheme,
  JetMatching_slowJetPower,
  LesHouches_idRenameBeams,
  LesHouches_setLeptonMass,
  LesHouches_setLifetime,
  Main_numberOfEvents,
  Main_numberOfSubruns,
  Main_spareMode1,
  Main_spareMode2,
  Main_spareMode3,
  Main_subrun,
  Main_timesAllowErrors,
  Merging_incompleteScalePrescrip,
  Merging_ktType,
  Merging_nJetMax,
  Merging_nJetMaxNLO,
  Merging_nQuarksMerge,
  Merging_nRecluster,
  Merging_nRequested,
  Merging_unorderedASscalePrescrip,
  Merging_unorderedPDFscalePrescrip,
  Merging_unorderedScalePrescrip,
  MiniStringFragmentation_nTry,
  MultipartonInteractions_alphaEMorder,
  MultipartonInteractions_alphaSorder,
  MultipartonInteractions_bProfile,
  MultipartonInteractions_bSelScale,
  MultipartonInteractions_enhanceScreening,
  MultipartonInteractions_nQuarkIn,
  MultipartonInteractions_nSample,
  MultipartonInteractions_processLevel,
  MultipartonInteractions_pTmaxMatch,
  MultipartonInteractions_rescatterMode,
  Next_numberCount,
  Next_numberShowEvent,
  Next_numberShowInfo,
  Next_numberShowLHA,
  Next_numberShowProcess,
  ParticleData_modeBreitWigner,
  PDF_PomSet,
  PDFinProcess_nQuarkIn,
  POWHEG_emitted,
  POWHEG_MPIveto,
  POWHEG_nFinal,
  POWHEG_pTdef,
  POWHEG_pTemt,
  POWHEG_pThard,
  POWHEG_QEDveto,
  POWHEG_veto,
  POWHEG_vetoCount,
  PromptPhoton_nQuarkLoop,
  Pythia_versionDate,
  Random_seed,
  RHadrons_idGluino,
  RHadrons_idSbottom,
  RHadrons_idStop,
  SigmaProcess_alphaEMorder,
  SigmaProcess_alphaSorder,
  SigmaProcess_factorScale1,
  SigmaProcess_factorScale2,
  SigmaProcess_factorScale3,
  SigmaProcess_factorScale3VV,
  SigmaProcess_renormScale1,
  SigmaProcess_renormScale2,
  SigmaProcess_renormScale3,
  SigmaProcess_renormScale3VV,
  SLHA_meMode,
  SLHA_readFrom,
  SLHA_verbose,
  SpaceShower_alphaEMorder,
  SpaceShower_alphaSorder,
  SpaceShower_nQuarkIn,
  SpaceShower_pTdampMatch,
  SpaceShower_pTmaxMatch,
  SpaceShower_weakShowerMode,
  StandardModel_alphaSnfmax,
  SUSY_idA,
  SUSY_idB,
  SUSY_sin2thetaWMode,
  TauDecays_externalMode,
  TauDecays_mode,
  TauDecays_tauMother,
  TimeShower_alphaEMorder,
  TimeShower_alphaSorder,
  TimeShower_globalRecoilMode,
  TimeShower_nGammaToLepton,
  TimeShower_nGammaToQuark,
  TimeShower_nGluonToQuark,
  TimeShower_nMaxGlobalBranch,
  TimeShower_nMaxGlobalRecoil,
  TimeShower_nPartonsInBorn,
  TimeShower_pTdampMatch,
  TimeShower_pTmaxMatch,
  TimeShower_weakShowerMode,
  TimeShower_weightGluonToQuark,
  Tune_ee,
  Tune_pp,
  Tune_preferLHAPDF,
  WeakZ0_gmZmode,
  Zprime_gmZmode,
  END,
};

enum class Param
{
  BeamRemnants_gluonPower,
  BeamRemnants_halfMassForKT,
  BeamRemnants_halfScaleForKT,
  BeamRemnants_primordialKThard,
  BeamRemnants_primordialKTremnant,
  BeamRemnants_primordialKTsoft,
  BeamRemnants_reducedKTatHighY,
  BeamRemnants_saturation,
  BeamRemnants_valenceDiqEnhance,
  BeamRemnants_valencePowerDinP,
  BeamRemnants_valencePowerMeson,
  BeamRemnants_valencePowerUinP,
  BeamRemnants_xGluonCutoff,
  Beams_eA,
  Beams_eB,
  Beams_eCM,
  Beams_maxDevA,
  Beams_maxDevB,
  Beams_maxDevTime,
  Beams_maxDevVertex,
  Beams_offsetTime,
  Beams_offsetVertexX,
  Beams_offsetVertexY,
  Beams_offsetVertexZ,
  Beams_pxA,
  Beams_pxB,
  Beams_pyA,
  Beams_pyB,
  Beams_pzA,
  Beams_pzB,
  Beams_sigmaPxA,
  Beams_sigmaPxB,
  Beams_sigmaPyA,
  Beams_sigmaPyB,
  Beams_sigmaPzA,
  Beams_sigmaPzB,
  Beams_sigmaTime,
  Beams_sigmaVertexX,
  Beams_sigmaVertexY,
  Beams_sigmaVertexZ,
  BoseEinstein_lambda,
  BoseEinstein_QRef,
  BoseEinstein_widthSep,
  Check_epTolErr,
  Check_epTolWarn,
  Check_mTolErr,
  Check_mTolWarn,
  ColourReconnection_blowR,
  ColourReconnection_blowT,
  ColourReconnection_dLambdaCut,
  ColourReconnection_fracGluon,
  ColourReconnection_fragmentationTime,
  ColourReconnection_junctionCorrection,
  ColourReconnection_kI,
  ColourReconnection_m0,
  ColourReconnection_m2Lambda,
  ColourReconnection_range,
  ColourReconnection_rHadron,
  ColourReconnection_timeDilationPar,
  ContactInteractions_Lambda,
  Diffraction_coreFraction,
  Diffraction_coreRadius,
  Diffraction_expPow,
  Diffraction_largeMassSuppress,
  Diffraction_MBRalpha,
  Diffraction_MBRbeta0,
  Diffraction_MBRdyminCD,
  Diffraction_MBRdyminCDflux,
  Diffraction_MBRdyminDD,
  Diffraction_MBRdyminDDflux,
  Diffraction_MBRdyminSD,
  Diffraction_MBRdyminSDflux,
  Diffraction_MBRdyminSigCD,
  Diffraction_MBRdyminSigDD,
  Diffraction_MBRdyminSigSD,
  Diffraction_MBRepsilon,
  Diffraction_MBRm2Min,
  Diffraction_MBRsigma0,
  Diffraction_mMinPert,
  Diffraction_mPowPomP,
  Diffraction_mRefPomP,
  Diffraction_mWidthPert,
  Diffraction_pickQuarkNorm,
  Diffraction_pickQuarkPower,
  Diffraction_PomFluxAlphaPrime,
  Diffraction_PomFluxEpsilon,
  Diffraction_primKTwidth,
  Diffraction_probMaxPert,
  Diffraction_sigmaRefPomP,
  ExcitedFermion_contactDec,
  ExcitedFermion_coupF,
  ExcitedFermion_coupFcol,
  ExcitedFermion_coupFprime,
  ExcitedFermion_Lambda,
  ExtraDimensionsGstar_Gbb,
  ExtraDimensionsGstar_Ggg,
  ExtraDimensionsGstar_Ggmgm,
  ExtraDimensionsGstar_Ghh,
  ExtraDimensionsGstar_Gll,
  ExtraDimensionsGstar_Gqq,
  ExtraDimensionsGstar_Gtt,
  ExtraDimensionsGstar_GWW,
  ExtraDimensionsGstar_GZZ,
  ExtraDimensionsGstar_kappaMG,
  ExtraDimensionsGstar_KKgbL,
  ExtraDimensionsGstar_KKgbR,
  ExtraDimensionsGstar_KKgqL,
  ExtraDimensionsGstar_KKgqR,
  ExtraDimensionsGstar_KKgtL,
  ExtraDimensionsGstar_KKgtR,
  ExtraDimensionsLED_c,
  ExtraDimensionsLED_g,
  ExtraDimensionsLED_LambdaT,
  ExtraDimensionsLED_MD,
  ExtraDimensionsLED_t,
  ExtraDimensionsTEV_mStar,
  ExtraDimensionsUnpart_dU,
  ExtraDimensionsUnpart_lambda,
  ExtraDimensionsUnpart_LambdaU,
  ExtraDimensionsUnpart_ratio,
  FourthGeneration_VcbPrime,
  FourthGeneration_VtbPrime,
  FourthGeneration_VtPrimeb,
  FourthGeneration_VtPrimebPrime,
  FourthGeneration_VtPrimed,
  FourthGeneration_VtPrimes,
  FourthGeneration_VubPrime,
  FragmentationSystems_mJoin,
  FragmentationSystems_mJoinJunction,
  HadronLevel_mStringMin,
  HadronScatter_j,
  HadronScatter_k,
  HadronScatter_N,
  HadronScatter_p,
  HadronScatter_rMax,
  HiddenValley_alphaFSR,
  HiddenValley_aLund,
  HiddenValley_bmqv2,
  HiddenValley_kappa,
  HiddenValley_kinMix,
  HiddenValley_probVector,
  HiddenValley_pTminFSR,
  HiddenValley_rFactqv,
  HiddenValley_sigmamqv,
  Higgs_wingsFac,
  HiggsA3_coup2d,
  HiggsA3_coup2H1H1,
  HiggsA3_coup2H1Z,
  HiggsA3_coup2H2Z,
  HiggsA3_coup2Hchg,
  HiggsA3_coup2HchgW,
  HiggsA3_coup2l,
  HiggsA3_coup2u,
  HiggsA3_coup2W,
  HiggsA3_coup2Z,
  HiggsA3_etaParity,
  HiggsA3_phiParity,
  HiggsH1_coup2d,
  HiggsH1_coup2Hchg,
  HiggsH1_coup2l,
  HiggsH1_coup2u,
  HiggsH1_coup2W,
  HiggsH1_coup2Z,
  HiggsH1_etaParity,
  HiggsH1_phiParity,
  HiggsH2_coup2A3A3,
  HiggsH2_coup2A3H1,
  HiggsH2_coup2d,
  HiggsH2_coup2H1H1,
  HiggsH2_coup2H1Z,
  HiggsH2_coup2Hchg,
  HiggsH2_coup2HchgW,
  HiggsH2_coup2l,
  HiggsH2_coup2u,
  HiggsH2_coup2W,
  HiggsH2_coup2Z,
  HiggsH2_etaParity,
  HiggsH2_phiParity,
  HiggsHchg_coup2H1W,
  HiggsHchg_coup2H2W,
  HiggsHchg_tanBeta,
  JetMatching_clFact,
  JetMatching_coneMatchHeavy,
  JetMatching_coneMatchLight,
  JetMatching_coneRadius,
  JetMatching_coneRadiusHeavy,
  JetMatching_etaJetMax,
  JetMatching_eTjetMin,
  JetMatching_eTseed,
  JetMatching_eTthreshold,
  JetMatching_qCut,
  JetMatching_qCutME,
  LeftRightSymmmetry_coupHee,
  LeftRightSymmmetry_coupHmue,
  LeftRightSymmmetry_coupHmumu,
  LeftRightSymmmetry_coupHtaue,
  LeftRightSymmmetry_coupHtaumu,
  LeftRightSymmmetry_coupHtautau,
  LeftRightSymmmetry_gL,
  LeftRightSymmmetry_gR,
  LeftRightSymmmetry_vL,
  LeptoQuark_kCoup,
  LesHouches_mRecalculate,
  Main_spareParm1,
  Main_spareParm2,
  Main_spareParm3,
  Merging_aCollFSR,
  Merging_aCollISR,
  Merging_Dparameter,
  Merging_dRijMS,
  Merging_fsrInRecNorm,
  Merging_kFactor0j,
  Merging_kFactor1j,
  Merging_kFactor2j,
  Merging_muFac,
  Merging_muFacInME,
  Merging_muRen,
  Merging_muRenInME,
  Merging_nonJoinedNorm,
  Merging_pTiMS,
  Merging_QijMS,
  Merging_scaleSeparationFactor,
  Merging_TMS,
  MultipartonInteractions_a1,
  MultipartonInteractions_alphaSvalue,
  MultipartonInteractions_coreFraction,
  MultipartonInteractions_coreRadius,
  MultipartonInteractions_deltaYRescatter,
  MultipartonInteractions_ecmPow,
  MultipartonInteractions_ecmRef,
  MultipartonInteractions_expPow,
  MultipartonInteractions_Kfactor,
  MultipartonInteractions_pT0Ref,
  MultipartonInteractions_pTmin,
  MultipartonInteractions_ySepRescatter,
  Onia_massSplit,
  ParticleData_alphaSvalueMRun,
  ParticleData_maxEnhanceBW,
  ParticleData_mbRun,
  ParticleData_mcRun,
  ParticleData_mdRun,
  ParticleData_msRun,
  ParticleData_mtRun,
  ParticleData_muRun,
  ParticleDecays_colRearrange,
  ParticleDecays_mSafety,
  ParticleDecays_multGoffset,
  ParticleDecays_multIncrease,
  ParticleDecays_multIncreaseWeak,
  ParticleDecays_multRefMass,
  ParticleDecays_rMax,
  ParticleDecays_sigmaSoft,
  ParticleDecays_tau0Max,
  ParticleDecays_tauMax,
  ParticleDecays_xBdMix,
  ParticleDecays_xBsMix,
  ParticleDecays_xyMax,
  ParticleDecays_zMax,
  PDF_PomGluonA,
  PDF_PomGluonB,
  PDF_PomQuarkA,
  PDF_PomQuarkB,
  PDF_PomQuarkFrac,
  PDF_PomRescale,
  PDF_PomStrangeSupp,
  PhaseSpace_bias2SelectionPow,
  PhaseSpace_bias2SelectionRef,
  PhaseSpace_mHatMax,
  PhaseSpace_mHatMaxSecond,
  PhaseSpace_mHatMin,
  PhaseSpace_mHatMinSecond,
  PhaseSpace_minWidthBreitWigners,
  PhaseSpace_pTHat3Max,
  PhaseSpace_pTHat3Min,
  PhaseSpace_pTHat5Max,
  PhaseSpace_pTHat5Min,
  PhaseSpace_pTHatMax,
  PhaseSpace_pTHatMaxSecond,
  PhaseSpace_pTHatMin,
  PhaseSpace_pTHatMinDiverge,
  PhaseSpace_pTHatMinSecond,
  PhaseSpace_RsepMin,
  Pythia_versionNumber,
  ResonanceWidths_minThreshold,
  ResonanceWidths_minWidth,
  RHadrons_diquarkSpin1,
  RHadrons_maxWidth,
  RHadrons_mCollapse,
  RHadrons_mOffsetCloud,
  RHadrons_probGluinoball,
  SigmaDiffractive_maxAX,
  SigmaDiffractive_maxAXB,
  SigmaDiffractive_maxXB,
  SigmaDiffractive_maxXX,
  SigmaElastic_bSlope,
  SigmaElastic_lambda,
  SigmaElastic_phaseConst,
  SigmaElastic_rho,
  SigmaElastic_tAbsMin,
  SigmaProcess_alphaSvalue,
  SigmaProcess_factorFixScale,
  SigmaProcess_factorMultFac,
  SigmaProcess_Kfactor,
  SigmaProcess_renormFixScale,
  SigmaProcess_renormMultFac,
  SigmaTotal_sigmaAX,
  SigmaTotal_sigmaAXB,
  SigmaTotal_sigmaAXB2TeV,
  SigmaTotal_sigmaEl,
  SigmaTotal_sigmaTot,
  SigmaTotal_sigmaXB,
  SigmaTotal_sigmaXX,
  SLHA_minMassSM,
  SpaceShower_alphaSuseCMW,
  SpaceShower_alphaSvalue,
  SpaceShower_ecmPow,
  SpaceShower_ecmRef,
  SpaceShower_factorMultFac,
  SpaceShower_fixedFacScale,
  SpaceShower_pT0Ref,
  SpaceShower_pTdampFudge,
  SpaceShower_pTmaxFudge,
  SpaceShower_pTmaxFudgeMPI,
  SpaceShower_pTmin,
  SpaceShower_pTminChgL,
  SpaceShower_pTminChgQ,
  SpaceShower_pTminWeak,
  SpaceShower_renormMultFac,
  SpaceShower_strengthIntAsym,
  StandardModel_alphaEM0,
  StandardModel_alphaEMmZ,
  StandardModel_GF,
  StandardModel_sin2thetaW,
  StandardModel_sin2thetaWbar,
  StandardModel_Vcb,
  StandardModel_Vcd,
  StandardModel_Vcs,
  StandardModel_Vtb,
  StandardModel_Vtd,
  StandardModel_Vts,
  StandardModel_Vub,
  StandardModel_Vud,
  StandardModel_Vus,
  StringFlav_decupletSup,
  StringFlav_etaPrimeSup,
  StringFlav_etaSup,
  StringFlav_heavyLeadingBSup,
  StringFlav_lightLeadingBSup,
  StringFlav_mesonBL1S0J1,
  StringFlav_mesonBL1S1J0,
  StringFlav_mesonBL1S1J1,
  StringFlav_mesonBL1S1J2,
  StringFlav_mesonBvector,
  StringFlav_mesonCL1S0J1,
  StringFlav_mesonCL1S1J0,
  StringFlav_mesonCL1S1J1,
  StringFlav_mesonCL1S1J2,
  StringFlav_mesonCvector,
  StringFlav_mesonSL1S0J1,
  StringFlav_mesonSL1S1J0,
  StringFlav_mesonSL1S1J1,
  StringFlav_mesonSL1S1J2,
  StringFlav_mesonSvector,
  StringFlav_mesonUDL1S0J1,
  StringFlav_mesonUDL1S1J0,
  StringFlav_mesonUDL1S1J1,
  StringFlav_mesonUDL1S1J2,
  StringFlav_mesonUDvector,
  StringFlav_popcornRate,
  StringFlav_popcornSmeson,
  StringFlav_popcornSpair,
  StringFlav_probQQ1toQQ0,
  StringFlav_probQQtoQ,
  StringFlav_probSQtoQQ,
  StringFlav_probStoUD,
  StringFlav_thetaL1S0J1,
  StringFlav_thetaL1S1J0,
  StringFlav_thetaL1S1J1,
  StringFlav_thetaL1S1J2,
  StringFlav_thetaPS,
  StringFlav_thetaV,
  StringFragmentation_eBothLeftJunction,
  StringFragmentation_eMaxLeftJunction,
  StringFragmentation_eMinLeftJunction,
  StringFragmentation_eNormJunction,
  StringFragmentation_stopMass,
  StringFragmentation_stopNewFlav,
  StringFragmentation_stopSmear,
  StringPT_enhancedFraction,
  StringPT_enhancedWidth,
  StringPT_sigma,
  StringZ_aExtraDiquark,
  StringZ_aExtraSQuark,
  StringZ_aLund,
  StringZ_aNonstandardB,
  StringZ_aNonstandardC,
  StringZ_aNonstandardH,
  StringZ_bLund,
  StringZ_bNonstandardB,
  StringZ_bNonstandardC,
  StringZ_bNonstandardH,
  StringZ_epsilonB,
  StringZ_epsilonC,
  StringZ_epsilonH,
  StringZ_rFactB,
  StringZ_rFactC,
  StringZ_rFactH,
  TauDecays_tauPolarization,
  TimeShower_alphaSvalue,
  TimeShower_factorMultFac,
  TimeShower_fixedFacScale,
  TimeShower_mMaxGamma,
  TimeShower_octetOniumColFac,
  TimeShower_octetOniumFraction,
  TimeShower_pTdampFudge,
  TimeShower_pTmaxFudge,
  TimeShower_pTmaxFudgeMPI,
  TimeShower_pTmin,
  TimeShower_pTminChgL,
  TimeShower_pTminChgQ,
  TimeShower_pTminWeak,
  TimeShower_renormMultFac,
  TimeShower_scaleGluonToQuark,
  WeakShower_enhancement,
  WeakShower_vetoWeakDeltaR,
  Wprime_al,
  Wprime_anglesWZ,
  Wprime_aq,
  Wprime_coup2WZ,
  Wprime_vl,
  Wprime_vq,
  Zprime_ab,
  Zprime_abPrime,
  Zprime_ac,
  Zprime_ad,
  Zprime_ae,
  Zprime_amu,
  Zprime_anglesWW,
  Zprime_anue,
  Zprime_anumu,
  Zprime_anutau,
  Zprime_anutauPrime,
  Zprime_as,
  Zprime_at,
  Zprime_atau,
  Zprime_atauPrime,
  Zprime_atPrime,
  Zprime_au,
  Zprime_coup2WW,
  Zprime_vb,
  Zprime_vbPrime,
  Zprime_vc,
  Zprime_vd,
  Zprime_ve,
  Zprime_vmu,
  Zprime_vnue,
  Zprime_vnumu,
  Zprime_vnutau,
  Zprime_vnutauPrime,
  Zprime_vs,
  Zprime_vt,
  Zprime_vtau,
  Zprime_vtauPrime,
  Zprime_vtPrime,
  Zprime_vu,
  END,
};

enum class Word
{
  Alpgen_file,
  Beams_LHEF,
  Beams_LHEFheader,
  Main_spareWord1,
  Main_spareWord2,
  Main_spareWord3,
  Merging_Process,
  PDF_p,
  PDF_pHardSet,
  PDF_pHardSetB,
  PDF_piSet,
  PDF_piSetB,
  PDF_pSet,
  PDF_pSetB,
  POWHEG_dir,
  SLHA_file,
  END,
};

enum class FlagList
{
  Bottomonium_gg2bbbar_3DJ__3DJ_1__g,
  Bottomonium_gg2bbbar_3DJ__3PJ_8__g,
  Bottomonium_gg2bbbar_3PJ__3PJ_1__g,
  Bottomonium_gg2bbbar_3PJ__3S1_8__g,
  Bottomonium_gg2bbbar_3S1__1S0_8__g,
  Bottomonium_gg2bbbar_3S1__3PJ_8__g,
  Bottomonium_gg2bbbar_3S1__3S1_1__g,
  Bottomonium_gg2bbbar_3S1__3S1_8__g,
  Bottomonium_qg2bbbar_3DJ__3PJ_8__q,
  Bottomonium_qg2bbbar_3PJ__3PJ_1__q,
  Bottomonium_qg2bbbar_3PJ__3S1_8__q,
  Bottomonium_qg2bbbar_3S1__1S0_8__q,
  Bottomonium_qg2bbbar_3S1__3PJ_8__q,
  Bottomonium_qg2bbbar_3S1__3S1_8__q,
  Bottomonium_qqbar2bbbar_3DJ__3PJ_8__g,
  Bottomonium_qqbar2bbbar_3PJ__3PJ_1__g,
  Bottomonium_qqbar2bbbar_3PJ__3S1_8__g,
  Bottomonium_qqbar2bbbar_3S1__1S0_8__g,
  Bottomonium_qqbar2bbbar_3S1__3PJ_8__g,
  Bottomonium_qqbar2bbbar_3S1__3S1_8__g,
  Charmonium_gg2ccbar_3DJ__3DJ_1__g,
  Charmonium_gg2ccbar_3DJ__3PJ_8__g,
  Charmonium_gg2ccbar_3PJ__3PJ_1__g,
  Charmonium_gg2ccbar_3PJ__3S1_8__g,
  Charmonium_gg2ccbar_3S1__1S0_8__g,
  Charmonium_gg2ccbar_3S1__3PJ_8__g,
  Charmonium_gg2ccbar_3S1__3S1_1__g,
  Charmonium_gg2ccbar_3S1__3S1_8__g,
  Charmonium_qg2ccbar_3DJ__3PJ_8__q,
  Charmonium_qg2ccbar_3PJ__3PJ_1__q,
  Charmonium_qg2ccbar_3PJ__3S1_8__q,
  Charmonium_qg2ccbar_3S1__1S0_8__q,
  Charmonium_qg2ccbar_3S1__3PJ_8__q,
  Charmonium_qg2ccbar_3S1__3S1_8__q,
  Charmonium_qqbar2ccbar_3DJ__3PJ_8__g,
  Charmonium_qqbar2ccbar_3PJ__3PJ_1__g,
  Charmonium_qqbar2ccbar_3PJ__3S1_8__g,
  Charmonium_qqbar2ccbar_3S1__1S0_8__g,
  Charmonium_qqbar2ccbar_3S1__3PJ_8__g,
  Charmonium_qqbar2ccbar_3S1__3S1_8__g,
  END,
};

enum class ModeList
{
  Bottomonium_states_3DJ_,
  Bottomonium_states_3PJ_,
  Bottomonium_states_3S1_,
  Charmonium_states_3DJ_,
  Charmonium_states_3PJ_,
  Charmonium_states_3S1_,
  SUSY_idVecA,
  SUSY_idVecB,
  END,
};

enum class ParamList
{
  Bottomonium_O_3DJ__3D1_1__,
  Bottomonium_O_3DJ__3P0_8__,
  Bottomonium_O_3PJ__3P0_1__,
  Bottomonium_O_3PJ__3S1_8__,
  Bottomonium_O_3S1__1S0_8__,
  Bottomonium_O_3S1__3P0_8__,
  Bottomonium_O_3S1__3S1_1__,
  Bottomonium_O_3S1__3S1_8__,
  Charmonium_O_3DJ__3D1_1__,
  Charmonium_O_3DJ__3P0_8__,
  Charmonium_O_3PJ__3P0_1__,
  Charmonium_O_3PJ__3S1_8__,
  Charmonium_O_3S1__1S0_8__,
  Charmonium_O_3S1__3P0_8__,
  Charmonium_O_3S1__3S1_1__,
  Charmonium_O_3S1__3S1_8__,
  StringFlav_probQQ1toQQ0join,
  END,
};

// class to store and access all settings

class Settings
{
public: // private data

  vector<bool> flagValue;
  vector<bool> flagDefault;
  vector<bool> modeFailOutOfRange;
  vector<int> modeValue;
  vector<int> modeDefault;
  vector<int> modeMax;
  vector<int> modeMin;
  vector<double> paramValue;
  vector<double> paramDefault;
  vector<double> paramMax;
  vector<double> paramMin;
  vector<string> wordValue;
  vector<string> wordDefault;
  vector<vector<bool>> flagListValue;
  vector<vector<bool>> flagListDefault;
  vector<vector<int>> modeListValue;
  vector<vector<int>> modeListDefault;
  vector<int> modeListMax;
  vector<int> modeListMin;
  vector<vector<double>> paramListValue;
  vector<vector<double>> paramListDefault;
  vector<double> paramListMax;
  vector<double> paramListMin;

public: // public functions

  // direct access to setting. WARNING: will ignore limits
  
  std::vector<bool>::reference operator[](const Flag flag);
  int& operator[](const Mode mode);
  double& operator[](const Param param);
  string& operator[](const Word word);
  vector<bool>& operator[](const FlagList flagList);
  vector<int>& operator[](const ModeList modeList);
  vector<double>& operator[](const ParamList paramList);

  // standard setting getters (const access only)
  
  bool get(const Flag flag);
  int get(const Mode mode);
  double get(const Param param);
  const string& get(const Word word);
  const vector<bool>& get(const FlagList flagList);
  const vector<int>& get(const ModeList modeList);
  const vector<double>& get(const ParamList paramList);

  // functions to get the default value (read only)

  bool getDefault(const Flag flag);
  int getDefault(const Mode mode);
  double getDefault(const Param param);
  const string& getDefault(const Word word);
  const vector<bool>& getDefault(const FlagList flagList);
  const vector<int>& getDefault(const ModeList modeList);
  const vector<double>& getDefault(const ParamList paramList);

  // modify any setting but respect the limits

  void set(Flag flag, bool value);
  void set(Mode mode, int value);
  void set(Param param, double value);
  void set(Word word, stringref value);
  void set(FlagList flagList, const vector<bool>& value);
  void set(ModeList modeList, const vector<int>& value);
  void set(ParamList paramList, const vector<double>& value);

  // restore default for a particular setting

  void restoreDefault(Flag flag);
  void restoreDefault(Mode mode);
  void restoreDefault(Param param);
  void restoreDefault(Word word);
  void restoreDefault(FlagList flagList);
  void restoreDefault(ModeList modeList);
  void restoreDefault(ParamList paramList);

  // restore default values for all settings

  void restoreDefault();

public: // private functions

  // adds the appropriate setting using string (default,min,max)

  void addFlag(stringref name, bool defaultVal);
  void addMode(stringref name, int defaultVal, bool hasMin, bool hasMax, int min, int max, bool boundsError);
  void addParam(stringref name, double defaultVal, bool hasMin, bool hasMax, double min, double max);
  void addWord(stringref name, stringref defaultVal);
  void addFlagList(stringref name, const vector<bool>& defaultVal);
  void addModeList(stringref name, const vector<int>& defaultVal, bool hasMin, bool hasMax, int min, int max);
  void addParamList(stringref name, const vector<double>& defaultVal, bool hasMin, bool hasMax, double min, double max);

public: // more functions

  // load all settings from XML database
  bool init(stringref startFile = "../xmldoc/Index.xml", bool append = false, bool reinit = false, ostream& os = cout);

  // update setting using a string (Format: name = value)
  bool readString(stringref line, bool warn = true, ostream& os = cout);

  // regulate level of printout by change of various settings
  void printQuiet(bool quiet);

  // write changed/all seting values to a file
  bool writeFile(stringref toFile, bool writeAll = false) ;

  // write changed/all seting values to stream
  bool writeFile(ostream& os=cout, bool writeAll = false);

  // wrtie table of settings (value|default|min|max) in lexigraphical order.
  void list(bool onlyChanged=true, string filter="", ostream& os=cout);

  // Restore all e+e- settings to their original values.
  void resetTuneEE();

  // Restore all pp settings to their original values.
  void resetTunePP();

  // Set the values related to a tune of e+e- data,
  // i.e. mainly for final-state radiation and hadronization.
  void initTuneEE();

  // Set the values related to a tune of pp/ppbar data,
  // i.e. mainly for initial-state radiation and multiparton interactions.
  void initTunePP();

};

} // end namespace Pythia8

#endif // Pythia8_Settings_H
