// Settings.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the Settings class.

#include "Pythia8/Settings.h"

// Allow string and character manipulation.
#include <cctype>
#include <regex>
#include <exception>

namespace Pythia8 {


// ----- name->index maps -----

const map<string, Flag> FlagMap = 
{
  { "Alpgen:setHeavyMasses", Flag::Alpgen_setHeavyMasses },
  { "Alpgen:setLightMasses", Flag::Alpgen_setLightMasses },
  { "Alpgen:setMLM", Flag::Alpgen_setMLM },
  { "Alpgen:setNjet", Flag::Alpgen_setNjet },
  { "BeamRemnants:allowBeamJunction", Flag::BeamRemnants_allowBeamJunction },
  { "BeamRemnants:allowJunction", Flag::BeamRemnants_allowJunction },
  { "BeamRemnants:beamJunction", Flag::BeamRemnants_beamJunction },
  { "BeamRemnants:primordialKT", Flag::BeamRemnants_primordialKT },
  { "BeamRemnants:rescatterRestoreY", Flag::BeamRemnants_rescatterRestoreY },
  { "Beams:allowMomentumSpread", Flag::Beams_allowMomentumSpread },
  { "Beams:allowVertexSpread", Flag::Beams_allowVertexSpread },
  { "Beams:newLHEFsameInit", Flag::Beams_newLHEFsameInit },
  { "Beams:readLHEFheaders", Flag::Beams_readLHEFheaders },
  { "Beams:setProductionScalesFromLHEF", Flag::Beams_setProductionScalesFromLHEF },
  { "Beams:strictLHEFscale", Flag::Beams_strictLHEFscale },
  { "BoseEinstein:Eta", Flag::BoseEinstein_Eta },
  { "BoseEinstein:Kaon", Flag::BoseEinstein_Kaon },
  { "BoseEinstein:Pion", Flag::BoseEinstein_Pion },
  { "Bottomonium:all", Flag::Bottomonium_all },
  { "Charmonium:all", Flag::Charmonium_all },
  { "Check:abortIfVeto", Flag::Check_abortIfVeto },
  { "Check:event", Flag::Check_event },
  { "Check:history", Flag::Check_history },
  { "Check:particleData", Flag::Check_particleData },
  { "ColourReconnection:allowDoubleJunRem", Flag::ColourReconnection_allowDoubleJunRem },
  { "ColourReconnection:allowJunctions", Flag::ColourReconnection_allowJunctions },
  { "ColourReconnection:forceHadronLevelCR", Flag::ColourReconnection_forceHadronLevelCR },
  { "ColourReconnection:forceResonance", Flag::ColourReconnection_forceResonance },
  { "ColourReconnection:lowerLambdaOnly", Flag::ColourReconnection_lowerLambdaOnly },
  { "ColourReconnection:reconnect", Flag::ColourReconnection_reconnect },
  { "ColourReconnection:sameNeighbourColours", Flag::ColourReconnection_sameNeighbourColours },
  { "ColourReconnection:singleReconnection", Flag::ColourReconnection_singleReconnection },
  { "ContactInteractions:QCffbar2eebar", Flag::ContactInteractions_QCffbar2eebar },
  { "ContactInteractions:QCffbar2mumubar", Flag::ContactInteractions_QCffbar2mumubar },
  { "ContactInteractions:QCffbar2tautaubar", Flag::ContactInteractions_QCffbar2tautaubar },
  { "ContactInteractions:QCqq2qq", Flag::ContactInteractions_QCqq2qq },
  { "ContactInteractions:QCqqbar2qqbar", Flag::ContactInteractions_QCqqbar2qqbar },
  { "Diffraction:doHard", Flag::Diffraction_doHard },
  { "ExcitedFermion:all", Flag::ExcitedFermion_all },
  { "ExcitedFermion:bg2bStar", Flag::ExcitedFermion_bg2bStar },
  { "ExcitedFermion:cg2cStar", Flag::ExcitedFermion_cg2cStar },
  { "ExcitedFermion:dg2dStar", Flag::ExcitedFermion_dg2dStar },
  { "ExcitedFermion:egm2eStar", Flag::ExcitedFermion_egm2eStar },
  { "ExcitedFermion:mugm2muStar", Flag::ExcitedFermion_mugm2muStar },
  { "ExcitedFermion:qq2bStarq", Flag::ExcitedFermion_qq2bStarq },
  { "ExcitedFermion:qq2cStarq", Flag::ExcitedFermion_qq2cStarq },
  { "ExcitedFermion:qq2dStarq", Flag::ExcitedFermion_qq2dStarq },
  { "ExcitedFermion:qq2sStarq", Flag::ExcitedFermion_qq2sStarq },
  { "ExcitedFermion:qq2uStarq", Flag::ExcitedFermion_qq2uStarq },
  { "ExcitedFermion:qqbar2eStare", Flag::ExcitedFermion_qqbar2eStare },
  { "ExcitedFermion:qqbar2eStareStar", Flag::ExcitedFermion_qqbar2eStareStar },
  { "ExcitedFermion:qqbar2muStarmu", Flag::ExcitedFermion_qqbar2muStarmu },
  { "ExcitedFermion:qqbar2muStarmuStar", Flag::ExcitedFermion_qqbar2muStarmuStar },
  { "ExcitedFermion:qqbar2nueStarnue", Flag::ExcitedFermion_qqbar2nueStarnue },
  { "ExcitedFermion:qqbar2nueStarnueStar", Flag::ExcitedFermion_qqbar2nueStarnueStar },
  { "ExcitedFermion:qqbar2numuStarnumu", Flag::ExcitedFermion_qqbar2numuStarnumu },
  { "ExcitedFermion:qqbar2numuStarnumuStar", Flag::ExcitedFermion_qqbar2numuStarnumuStar },
  { "ExcitedFermion:qqbar2nutauStarnutau", Flag::ExcitedFermion_qqbar2nutauStarnutau },
  { "ExcitedFermion:qqbar2nutauStarnutauStar", Flag::ExcitedFermion_qqbar2nutauStarnutauStar },
  { "ExcitedFermion:qqbar2tauStartau", Flag::ExcitedFermion_qqbar2tauStartau },
  { "ExcitedFermion:qqbar2tauStartauStar", Flag::ExcitedFermion_qqbar2tauStartauStar },
  { "ExcitedFermion:sg2sStar", Flag::ExcitedFermion_sg2sStar },
  { "ExcitedFermion:taugm2tauStar", Flag::ExcitedFermion_taugm2tauStar },
  { "ExcitedFermion:ug2uStar", Flag::ExcitedFermion_ug2uStar },
  { "ExtraDimensionsG*:all", Flag::ExtraDimensionsGstar_all },
  { "ExtraDimensionsG*:ffbar2G*", Flag::ExtraDimensionsGstar_ffbar2Gstar },
  { "ExtraDimensionsG*:gg2G*", Flag::ExtraDimensionsGstar_gg2Gstar },
  { "ExtraDimensionsG*:gg2G*g", Flag::ExtraDimensionsGstar_gg2Gstarg },
  { "ExtraDimensionsG*:qg2G*q", Flag::ExtraDimensionsGstar_qg2Gstarq },
  { "ExtraDimensionsG*:qqbar2G*g", Flag::ExtraDimensionsGstar_qqbar2Gstarg },
  { "ExtraDimensionsG*:qqbar2KKgluon*", Flag::ExtraDimensionsGstar_qqbar2KKgluonstar },
  { "ExtraDimensionsG*:SMinBulk", Flag::ExtraDimensionsGstar_SMinBulk },
  { "ExtraDimensionsG*:VLVL", Flag::ExtraDimensionsGstar_VLVL },
  { "ExtraDimensionsLED:dijets", Flag::ExtraDimensionsLED_dijets },
  { "ExtraDimensionsLED:ffbar2gammagamma", Flag::ExtraDimensionsLED_ffbar2gammagamma },
  { "ExtraDimensionsLED:ffbar2Ggamma", Flag::ExtraDimensionsLED_ffbar2Ggamma },
  { "ExtraDimensionsLED:ffbar2GZ", Flag::ExtraDimensionsLED_ffbar2GZ },
  { "ExtraDimensionsLED:ffbar2llbar", Flag::ExtraDimensionsLED_ffbar2llbar },
  { "ExtraDimensionsLED:gg2DJgg", Flag::ExtraDimensionsLED_gg2DJgg },
  { "ExtraDimensionsLED:gg2DJqqbar", Flag::ExtraDimensionsLED_gg2DJqqbar },
  { "ExtraDimensionsLED:gg2gammagamma", Flag::ExtraDimensionsLED_gg2gammagamma },
  { "ExtraDimensionsLED:gg2Gg", Flag::ExtraDimensionsLED_gg2Gg },
  { "ExtraDimensionsLED:gg2llbar", Flag::ExtraDimensionsLED_gg2llbar },
  { "ExtraDimensionsLED:GravScalar", Flag::ExtraDimensionsLED_GravScalar },
  { "ExtraDimensionsLED:monojet", Flag::ExtraDimensionsLED_monojet },
  { "ExtraDimensionsLED:qg2DJqg", Flag::ExtraDimensionsLED_qg2DJqg },
  { "ExtraDimensionsLED:qg2Gq", Flag::ExtraDimensionsLED_qg2Gq },
  { "ExtraDimensionsLED:qq2DJqq", Flag::ExtraDimensionsLED_qq2DJqq },
  { "ExtraDimensionsLED:qqbar2DJgg", Flag::ExtraDimensionsLED_qqbar2DJgg },
  { "ExtraDimensionsLED:qqbar2DJqqbarNew", Flag::ExtraDimensionsLED_qqbar2DJqqbarNew },
  { "ExtraDimensionsLED:qqbar2Gg", Flag::ExtraDimensionsLED_qqbar2Gg },
  { "ExtraDimensionsTEV:ffbar2bbbar", Flag::ExtraDimensionsTEV_ffbar2bbbar },
  { "ExtraDimensionsTEV:ffbar2ccbar", Flag::ExtraDimensionsTEV_ffbar2ccbar },
  { "ExtraDimensionsTEV:ffbar2ddbar", Flag::ExtraDimensionsTEV_ffbar2ddbar },
  { "ExtraDimensionsTEV:ffbar2e+e-", Flag::ExtraDimensionsTEV_ffbar2epluseminus },
  { "ExtraDimensionsTEV:ffbar2mu+mu-", Flag::ExtraDimensionsTEV_ffbar2muplusmuminus },
  { "ExtraDimensionsTEV:ffbar2nuenuebar", Flag::ExtraDimensionsTEV_ffbar2nuenuebar },
  { "ExtraDimensionsTEV:ffbar2numunumubar", Flag::ExtraDimensionsTEV_ffbar2numunumubar },
  { "ExtraDimensionsTEV:ffbar2nutaunutaubar", Flag::ExtraDimensionsTEV_ffbar2nutaunutaubar },
  { "ExtraDimensionsTEV:ffbar2ssbar", Flag::ExtraDimensionsTEV_ffbar2ssbar },
  { "ExtraDimensionsTEV:ffbar2tau+tau-", Flag::ExtraDimensionsTEV_ffbar2tauplustauminus },
  { "ExtraDimensionsTEV:ffbar2ttbar", Flag::ExtraDimensionsTEV_ffbar2ttbar },
  { "ExtraDimensionsTEV:ffbar2uubar", Flag::ExtraDimensionsTEV_ffbar2uubar },
  { "ExtraDimensionsUnpart:ffbar2gammagamma", Flag::ExtraDimensionsUnpart_ffbar2gammagamma },
  { "ExtraDimensionsUnpart:ffbar2llbar", Flag::ExtraDimensionsUnpart_ffbar2llbar },
  { "ExtraDimensionsUnpart:ffbar2Ugamma", Flag::ExtraDimensionsUnpart_ffbar2Ugamma },
  { "ExtraDimensionsUnpart:ffbar2UZ", Flag::ExtraDimensionsUnpart_ffbar2UZ },
  { "ExtraDimensionsUnpart:gg2gammagamma", Flag::ExtraDimensionsUnpart_gg2gammagamma },
  { "ExtraDimensionsUnpart:gg2llbar", Flag::ExtraDimensionsUnpart_gg2llbar },
  { "ExtraDimensionsUnpart:gg2Ug", Flag::ExtraDimensionsUnpart_gg2Ug },
  { "ExtraDimensionsUnpart:monojet", Flag::ExtraDimensionsUnpart_monojet },
  { "ExtraDimensionsUnpart:qg2Uq", Flag::ExtraDimensionsUnpart_qg2Uq },
  { "ExtraDimensionsUnpart:qqbar2Ug", Flag::ExtraDimensionsUnpart_qqbar2Ug },
  { "FourthBottom:all", Flag::FourthBottom_all },
  { "FourthBottom:ffbar2bPrimebPrimebar(s:gmZ)", Flag::FourthBottom_ffbar2bPrimebPrimebar_s_gmZ_ },
  { "FourthBottom:ffbar2bPrimeqbar(s:W)", Flag::FourthBottom_ffbar2bPrimeqbar_s_W_ },
  { "FourthBottom:ffbar2bPrimetbar(s:W)", Flag::FourthBottom_ffbar2bPrimetbar_s_W_ },
  { "FourthBottom:gg2bPrimebPrimebar", Flag::FourthBottom_gg2bPrimebPrimebar },
  { "FourthBottom:qq2bPrimeq(t:W)", Flag::FourthBottom_qq2bPrimeq_t_W_ },
  { "FourthBottom:qqbar2bPrimebPrimebar", Flag::FourthBottom_qqbar2bPrimebPrimebar },
  { "FourthPair:ffbar2tauPrimenuPrimebar(s:W)", Flag::FourthPair_ffbar2tauPrimenuPrimebar_s_W_ },
  { "FourthPair:ffbar2tPrimebPrimebar(s:W)", Flag::FourthPair_ffbar2tPrimebPrimebar_s_W_ },
  { "FourthTop:all", Flag::FourthTop_all },
  { "FourthTop:ffbar2tPrimeqbar(s:W)", Flag::FourthTop_ffbar2tPrimeqbar_s_W_ },
  { "FourthTop:ffbar2tPrimetPrimebar(s:gmZ)", Flag::FourthTop_ffbar2tPrimetPrimebar_s_gmZ_ },
  { "FourthTop:gg2tPrimetPrimebar", Flag::FourthTop_gg2tPrimetPrimebar },
  { "FourthTop:qq2tPrimeq(t:W)", Flag::FourthTop_qq2tPrimeq_t_W_ },
  { "FourthTop:qqbar2tPrimetPrimebar", Flag::FourthTop_qqbar2tPrimetPrimebar },
  { "HadronLevel:all", Flag::HadronLevel_all },
  { "HadronLevel:BoseEinstein", Flag::HadronLevel_BoseEinstein },
  { "HadronLevel:Decay", Flag::HadronLevel_Decay },
  { "HadronLevel:Hadronize", Flag::HadronLevel_Hadronize },
  { "HadronScatter:afterDecay", Flag::HadronScatter_afterDecay },
  { "HadronScatter:allowDecayProd", Flag::HadronScatter_allowDecayProd },
  { "HadronScatter:scatter", Flag::HadronScatter_scatter },
  { "HadronScatter:scatterRepeat", Flag::HadronScatter_scatterRepeat },
  { "HadronScatter:tile", Flag::HadronScatter_tile },
  { "HardQCD:3parton", Flag::HardQCD_3parton },
  { "HardQCD:all", Flag::HardQCD_all },
  { "HardQCD:gg2bbbar", Flag::HardQCD_gg2bbbar },
  { "HardQCD:gg2ccbar", Flag::HardQCD_gg2ccbar },
  { "HardQCD:gg2gg", Flag::HardQCD_gg2gg },
  { "HardQCD:gg2ggg", Flag::HardQCD_gg2ggg },
  { "HardQCD:gg2qqbar", Flag::HardQCD_gg2qqbar },
  { "HardQCD:gg2qqbarg", Flag::HardQCD_gg2qqbarg },
  { "HardQCD:hardbbbar", Flag::HardQCD_hardbbbar },
  { "HardQCD:hardccbar", Flag::HardQCD_hardccbar },
  { "HardQCD:qg2qg", Flag::HardQCD_qg2qg },
  { "HardQCD:qg2qgg", Flag::HardQCD_qg2qgg },
  { "HardQCD:qg2qqqbarDiff", Flag::HardQCD_qg2qqqbarDiff },
  { "HardQCD:qg2qqqbarSame", Flag::HardQCD_qg2qqqbarSame },
  { "HardQCD:qq2qq", Flag::HardQCD_qq2qq },
  { "HardQCD:qq2qqgDiff", Flag::HardQCD_qq2qqgDiff },
  { "HardQCD:qq2qqgSame", Flag::HardQCD_qq2qqgSame },
  { "HardQCD:qqbar2bbbar", Flag::HardQCD_qqbar2bbbar },
  { "HardQCD:qqbar2ccbar", Flag::HardQCD_qqbar2ccbar },
  { "HardQCD:qqbar2gg", Flag::HardQCD_qqbar2gg },
  { "HardQCD:qqbar2ggg", Flag::HardQCD_qqbar2ggg },
  { "HardQCD:qqbar2qqbargDiff", Flag::HardQCD_qqbar2qqbargDiff },
  { "HardQCD:qqbar2qqbargSame", Flag::HardQCD_qqbar2qqbargSame },
  { "HardQCD:qqbar2qqbarNew", Flag::HardQCD_qqbar2qqbarNew },
  { "HiddenValley:all", Flag::HiddenValley_all },
  { "HiddenValley:doKinMix", Flag::HiddenValley_doKinMix },
  { "HiddenValley:ffbar2BvBvbar", Flag::HiddenValley_ffbar2BvBvbar },
  { "HiddenValley:ffbar2CvCvbar", Flag::HiddenValley_ffbar2CvCvbar },
  { "HiddenValley:ffbar2DvDvbar", Flag::HiddenValley_ffbar2DvDvbar },
  { "HiddenValley:ffbar2EvEvbar", Flag::HiddenValley_ffbar2EvEvbar },
  { "HiddenValley:ffbar2MUvMUvbar", Flag::HiddenValley_ffbar2MUvMUvbar },
  { "HiddenValley:ffbar2nuEvnuEvbar", Flag::HiddenValley_ffbar2nuEvnuEvbar },
  { "HiddenValley:ffbar2nuMUvnuMUvbar", Flag::HiddenValley_ffbar2nuMUvnuMUvbar },
  { "HiddenValley:ffbar2nuTAUvnuTAUvbar", Flag::HiddenValley_ffbar2nuTAUvnuTAUvbar },
  { "HiddenValley:ffbar2SvSvbar", Flag::HiddenValley_ffbar2SvSvbar },
  { "HiddenValley:ffbar2TAUvTAUvbar", Flag::HiddenValley_ffbar2TAUvTAUvbar },
  { "HiddenValley:ffbar2TvTvbar", Flag::HiddenValley_ffbar2TvTvbar },
  { "HiddenValley:ffbar2UvUvbar", Flag::HiddenValley_ffbar2UvUvbar },
  { "HiddenValley:ffbar2Zv", Flag::HiddenValley_ffbar2Zv },
  { "HiddenValley:fragment", Flag::HiddenValley_fragment },
  { "HiddenValley:FSR", Flag::HiddenValley_FSR },
  { "HiddenValley:gg2BvBvbar", Flag::HiddenValley_gg2BvBvbar },
  { "HiddenValley:gg2CvCvbar", Flag::HiddenValley_gg2CvCvbar },
  { "HiddenValley:gg2DvDvbar", Flag::HiddenValley_gg2DvDvbar },
  { "HiddenValley:gg2SvSvbar", Flag::HiddenValley_gg2SvSvbar },
  { "HiddenValley:gg2TvTvbar", Flag::HiddenValley_gg2TvTvbar },
  { "HiddenValley:gg2UvUvbar", Flag::HiddenValley_gg2UvUvbar },
  { "HiddenValley:qqbar2BvBvbar", Flag::HiddenValley_qqbar2BvBvbar },
  { "HiddenValley:qqbar2CvCvbar", Flag::HiddenValley_qqbar2CvCvbar },
  { "HiddenValley:qqbar2DvDvbar", Flag::HiddenValley_qqbar2DvDvbar },
  { "HiddenValley:qqbar2SvSvbar", Flag::HiddenValley_qqbar2SvSvbar },
  { "HiddenValley:qqbar2TvTvbar", Flag::HiddenValley_qqbar2TvTvbar },
  { "HiddenValley:qqbar2UvUvbar", Flag::HiddenValley_qqbar2UvUvbar },
  { "Higgs:clipWings", Flag::Higgs_clipWings },
  { "Higgs:cubicWidth", Flag::Higgs_cubicWidth },
  { "Higgs:runningLoopMass", Flag::Higgs_runningLoopMass },
  { "Higgs:useBSM", Flag::Higgs_useBSM },
  { "HiggsBSM:all", Flag::HiggsBSM_all },
  { "HiggsBSM:allA3", Flag::HiggsBSM_allA3 },
  { "HiggsBSM:allH+-", Flag::HiggsBSM_allHplusminus },
  { "HiggsBSM:allH1", Flag::HiggsBSM_allH1 },
  { "HiggsBSM:allH2", Flag::HiggsBSM_allH2 },
  { "HiggsBSM:allHpair", Flag::HiggsBSM_allHpair },
  { "HiggsBSM:bg2H+-t", Flag::HiggsBSM_bg2Hplusminust },
  { "HiggsBSM:ff2A3ff(t:WW)", Flag::HiggsBSM_ff2A3ff_t_WW_ },
  { "HiggsBSM:ff2A3ff(t:ZZ)", Flag::HiggsBSM_ff2A3ff_t_ZZ_ },
  { "HiggsBSM:ff2H1ff(t:WW)", Flag::HiggsBSM_ff2H1ff_t_WW_ },
  { "HiggsBSM:ff2H1ff(t:ZZ)", Flag::HiggsBSM_ff2H1ff_t_ZZ_ },
  { "HiggsBSM:ff2H2ff(t:WW)", Flag::HiggsBSM_ff2H2ff_t_WW_ },
  { "HiggsBSM:ff2H2ff(t:ZZ)", Flag::HiggsBSM_ff2H2ff_t_ZZ_ },
  { "HiggsBSM:ffbar2A3", Flag::HiggsBSM_ffbar2A3 },
  { "HiggsBSM:ffbar2A3H1", Flag::HiggsBSM_ffbar2A3H1 },
  { "HiggsBSM:ffbar2A3H2", Flag::HiggsBSM_ffbar2A3H2 },
  { "HiggsBSM:ffbar2A3W", Flag::HiggsBSM_ffbar2A3W },
  { "HiggsBSM:ffbar2A3Z", Flag::HiggsBSM_ffbar2A3Z },
  { "HiggsBSM:ffbar2H+-", Flag::HiggsBSM_ffbar2Hplusminus },
  { "HiggsBSM:ffbar2H+-H1", Flag::HiggsBSM_ffbar2HplusminusH1 },
  { "HiggsBSM:ffbar2H+-H2", Flag::HiggsBSM_ffbar2HplusminusH2 },
  { "HiggsBSM:ffbar2H+H-", Flag::HiggsBSM_ffbar2HplusHminus },
  { "HiggsBSM:ffbar2H1", Flag::HiggsBSM_ffbar2H1 },
  { "HiggsBSM:ffbar2H1W", Flag::HiggsBSM_ffbar2H1W },
  { "HiggsBSM:ffbar2H1Z", Flag::HiggsBSM_ffbar2H1Z },
  { "HiggsBSM:ffbar2H2", Flag::HiggsBSM_ffbar2H2 },
  { "HiggsBSM:ffbar2H2W", Flag::HiggsBSM_ffbar2H2W },
  { "HiggsBSM:ffbar2H2Z", Flag::HiggsBSM_ffbar2H2Z },
  { "HiggsBSM:gg2A3", Flag::HiggsBSM_gg2A3 },
  { "HiggsBSM:gg2A3bbbar", Flag::HiggsBSM_gg2A3bbbar },
  { "HiggsBSM:gg2A3g(l:t)", Flag::HiggsBSM_gg2A3g_l_t_ },
  { "HiggsBSM:gg2A3ttbar", Flag::HiggsBSM_gg2A3ttbar },
  { "HiggsBSM:gg2H1", Flag::HiggsBSM_gg2H1 },
  { "HiggsBSM:gg2H1bbbar", Flag::HiggsBSM_gg2H1bbbar },
  { "HiggsBSM:gg2H1g(l:t)", Flag::HiggsBSM_gg2H1g_l_t_ },
  { "HiggsBSM:gg2H1ttbar", Flag::HiggsBSM_gg2H1ttbar },
  { "HiggsBSM:gg2H2", Flag::HiggsBSM_gg2H2 },
  { "HiggsBSM:gg2H2bbbar", Flag::HiggsBSM_gg2H2bbbar },
  { "HiggsBSM:gg2H2g(l:t)", Flag::HiggsBSM_gg2H2g_l_t_ },
  { "HiggsBSM:gg2H2ttbar", Flag::HiggsBSM_gg2H2ttbar },
  { "HiggsBSM:gmgm2A3", Flag::HiggsBSM_gmgm2A3 },
  { "HiggsBSM:gmgm2H1", Flag::HiggsBSM_gmgm2H1 },
  { "HiggsBSM:gmgm2H2", Flag::HiggsBSM_gmgm2H2 },
  { "HiggsBSM:qg2A3q", Flag::HiggsBSM_qg2A3q },
  { "HiggsBSM:qg2A3q(l:t)", Flag::HiggsBSM_qg2A3q_l_t_ },
  { "HiggsBSM:qg2H1q", Flag::HiggsBSM_qg2H1q },
  { "HiggsBSM:qg2H1q(l:t)", Flag::HiggsBSM_qg2H1q_l_t_ },
  { "HiggsBSM:qg2H2q", Flag::HiggsBSM_qg2H2q },
  { "HiggsBSM:qg2H2q(l:t)", Flag::HiggsBSM_qg2H2q_l_t_ },
  { "HiggsBSM:qqbar2A3bbbar", Flag::HiggsBSM_qqbar2A3bbbar },
  { "HiggsBSM:qqbar2A3g(l:t)", Flag::HiggsBSM_qqbar2A3g_l_t_ },
  { "HiggsBSM:qqbar2A3ttbar", Flag::HiggsBSM_qqbar2A3ttbar },
  { "HiggsBSM:qqbar2H1bbbar", Flag::HiggsBSM_qqbar2H1bbbar },
  { "HiggsBSM:qqbar2H1g(l:t)", Flag::HiggsBSM_qqbar2H1g_l_t_ },
  { "HiggsBSM:qqbar2H1ttbar", Flag::HiggsBSM_qqbar2H1ttbar },
  { "HiggsBSM:qqbar2H2bbbar", Flag::HiggsBSM_qqbar2H2bbbar },
  { "HiggsBSM:qqbar2H2g(l:t)", Flag::HiggsBSM_qqbar2H2g_l_t_ },
  { "HiggsBSM:qqbar2H2ttbar", Flag::HiggsBSM_qqbar2H2ttbar },
  { "HiggsSM:all", Flag::HiggsSM_all },
  { "HiggsSM:ff2Hff(t:WW)", Flag::HiggsSM_ff2Hff_t_WW_ },
  { "HiggsSM:ff2Hff(t:ZZ)", Flag::HiggsSM_ff2Hff_t_ZZ_ },
  { "HiggsSM:ffbar2H", Flag::HiggsSM_ffbar2H },
  { "HiggsSM:ffbar2HW", Flag::HiggsSM_ffbar2HW },
  { "HiggsSM:ffbar2HZ", Flag::HiggsSM_ffbar2HZ },
  { "HiggsSM:gg2H", Flag::HiggsSM_gg2H },
  { "HiggsSM:gg2Hbbbar", Flag::HiggsSM_gg2Hbbbar },
  { "HiggsSM:gg2Hg(l:t)", Flag::HiggsSM_gg2Hg_l_t_ },
  { "HiggsSM:gg2Httbar", Flag::HiggsSM_gg2Httbar },
  { "HiggsSM:gmgm2H", Flag::HiggsSM_gmgm2H },
  { "HiggsSM:NLOWidths", Flag::HiggsSM_NLOWidths },
  { "HiggsSM:qg2Hq", Flag::HiggsSM_qg2Hq },
  { "HiggsSM:qg2Hq(l:t)", Flag::HiggsSM_qg2Hq_l_t_ },
  { "HiggsSM:qqbar2Hbbbar", Flag::HiggsSM_qqbar2Hbbbar },
  { "HiggsSM:qqbar2Hg(l:t)", Flag::HiggsSM_qqbar2Hg_l_t_ },
  { "HiggsSM:qqbar2Httbar", Flag::HiggsSM_qqbar2Httbar },
  { "Init:showAllParticleData", Flag::Init_showAllParticleData },
  { "Init:showAllSettings", Flag::Init_showAllSettings },
  { "Init:showChangedParticleData", Flag::Init_showChangedParticleData },
  { "Init:showChangedResonanceData", Flag::Init_showChangedResonanceData },
  { "Init:showChangedSettings", Flag::Init_showChangedSettings },
  { "Init:showMultipartonInteractions", Flag::Init_showMultipartonInteractions },
  { "Init:showProcesses", Flag::Init_showProcesses },
  { "JetMatching:doFxFx", Flag::JetMatching_doFxFx },
  { "JetMatching:doShowerKt", Flag::JetMatching_doShowerKt },
  { "JetMatching:merge", Flag::JetMatching_merge },
  { "JetMatching:setMad", Flag::JetMatching_setMad },
  { "LeftRightSymmmetry:all", Flag::LeftRightSymmmetry_all },
  { "LeftRightSymmmetry:ff2HLff", Flag::LeftRightSymmmetry_ff2HLff },
  { "LeftRightSymmmetry:ff2HRff", Flag::LeftRightSymmmetry_ff2HRff },
  { "LeftRightSymmmetry:ffbar2HLHL", Flag::LeftRightSymmmetry_ffbar2HLHL },
  { "LeftRightSymmmetry:ffbar2HRHR", Flag::LeftRightSymmmetry_ffbar2HRHR },
  { "LeftRightSymmmetry:ffbar2WR", Flag::LeftRightSymmmetry_ffbar2WR },
  { "LeftRightSymmmetry:ffbar2ZR", Flag::LeftRightSymmmetry_ffbar2ZR },
  { "LeftRightSymmmetry:lgm2HLe", Flag::LeftRightSymmmetry_lgm2HLe },
  { "LeftRightSymmmetry:lgm2HLmu", Flag::LeftRightSymmmetry_lgm2HLmu },
  { "LeftRightSymmmetry:lgm2HLtau", Flag::LeftRightSymmmetry_lgm2HLtau },
  { "LeftRightSymmmetry:lgm2HRe", Flag::LeftRightSymmmetry_lgm2HRe },
  { "LeftRightSymmmetry:lgm2HRmu", Flag::LeftRightSymmmetry_lgm2HRmu },
  { "LeftRightSymmmetry:lgm2HRtau", Flag::LeftRightSymmmetry_lgm2HRtau },
  { "LeftRightSymmmetry:ll2HL", Flag::LeftRightSymmmetry_ll2HL },
  { "LeftRightSymmmetry:ll2HR", Flag::LeftRightSymmmetry_ll2HR },
  { "LeptoQuark:all", Flag::LeptoQuark_all },
  { "LeptoQuark:gg2LQLQbar", Flag::LeptoQuark_gg2LQLQbar },
  { "LeptoQuark:qg2LQl", Flag::LeptoQuark_qg2LQl },
  { "LeptoQuark:ql2LQ", Flag::LeptoQuark_ql2LQ },
  { "LeptoQuark:qqbar2LQLQbar", Flag::LeptoQuark_qqbar2LQLQbar },
  { "LesHouches:matchInOut", Flag::LesHouches_matchInOut },
  { "Main:LHEFskipInit", Flag::Main_LHEFskipInit },
  { "Main:spareFlag1", Flag::Main_spareFlag1 },
  { "Main:spareFlag2", Flag::Main_spareFlag2 },
  { "Main:spareFlag3", Flag::Main_spareFlag3 },
  { "Merging:allowColourShuffling", Flag::Merging_allowColourShuffling },
  { "Merging:allowIncompleteHistoriesInReal", Flag::Merging_allowIncompleteHistoriesInReal },
  { "Merging:allowSQCDClustering", Flag::Merging_allowSQCDClustering },
  { "Merging:allowWClustering", Flag::Merging_allowWClustering },
  { "Merging:doCutBasedMerging", Flag::Merging_doCutBasedMerging },
  { "Merging:doKTMerging", Flag::Merging_doKTMerging },
  { "Merging:doMGMerging", Flag::Merging_doMGMerging },
  { "Merging:doNL3Loop", Flag::Merging_doNL3Loop },
  { "Merging:doNL3Subt", Flag::Merging_doNL3Subt },
  { "Merging:doNL3Tree", Flag::Merging_doNL3Tree },
  { "Merging:doPTLundMerging", Flag::Merging_doPTLundMerging },
  { "Merging:doUMEPSSubt", Flag::Merging_doUMEPSSubt },
  { "Merging:doUMEPSTree", Flag::Merging_doUMEPSTree },
  { "Merging:doUNLOPSLoop", Flag::Merging_doUNLOPSLoop },
  { "Merging:doUNLOPSSubt", Flag::Merging_doUNLOPSSubt },
  { "Merging:doUNLOPSSubtNLO", Flag::Merging_doUNLOPSSubtNLO },
  { "Merging:doUNLOPSTilde", Flag::Merging_doUNLOPSTilde },
  { "Merging:doUNLOPSTree", Flag::Merging_doUNLOPSTree },
  { "Merging:doUserMerging", Flag::Merging_doUserMerging },
  { "Merging:doXSectionEstimate", Flag::Merging_doXSectionEstimate },
  { "Merging:enforceCutOnLHE", Flag::Merging_enforceCutOnLHE },
  { "Merging:enforceStrongOrdering", Flag::Merging_enforceStrongOrdering },
  { "Merging:includeMassive", Flag::Merging_includeMassive },
  { "Merging:includeRedundant", Flag::Merging_includeRedundant },
  { "Merging:mayRemoveDecayProducts", Flag::Merging_mayRemoveDecayProducts },
  { "Merging:orderInRapidity", Flag::Merging_orderInRapidity },
  { "Merging:pickByFullP", Flag::Merging_pickByFullP },
  { "Merging:pickByPoPT2", Flag::Merging_pickByPoPT2 },
  { "Merging:pickBySumPT", Flag::Merging_pickBySumPT },
  { "Merging:usePythiaQFacHard", Flag::Merging_usePythiaQFacHard },
  { "Merging:usePythiaQRenHard", Flag::Merging_usePythiaQRenHard },
  { "Merging:useShowerPlugin", Flag::Merging_useShowerPlugin },
  { "MultipartonInteractions:allowDoubleRescatter", Flag::MultipartonInteractions_allowDoubleRescatter },
  { "MultipartonInteractions:allowRescatter", Flag::MultipartonInteractions_allowRescatter },
  { "NewGaugeBoson:ffbar2gmZZprime", Flag::NewGaugeBoson_ffbar2gmZZprime },
  { "NewGaugeBoson:ffbar2R0", Flag::NewGaugeBoson_ffbar2R0 },
  { "NewGaugeBoson:ffbar2Wprime", Flag::NewGaugeBoson_ffbar2Wprime },
  { "Next:showMothersAndDaughters", Flag::Next_showMothersAndDaughters },
  { "Next:showScaleAndVertex", Flag::Next_showScaleAndVertex },
  { "Onia:all", Flag::Onia_all },
  { "Onia:all(3DJ)", Flag::Onia_all_3DJ_ },
  { "Onia:all(3PJ)", Flag::Onia_all_3PJ_ },
  { "Onia:all(3S1)", Flag::Onia_all_3S1_ },
  { "Onia:forceMassSplit", Flag::Onia_forceMassSplit },
  { "ParticleDecays:allowPhotonRadiation", Flag::ParticleDecays_allowPhotonRadiation },
  { "ParticleDecays:FSRinDecays", Flag::ParticleDecays_FSRinDecays },
  { "ParticleDecays:limitCylinder", Flag::ParticleDecays_limitCylinder },
  { "ParticleDecays:limitRadius", Flag::ParticleDecays_limitRadius },
  { "ParticleDecays:limitTau", Flag::ParticleDecays_limitTau },
  { "ParticleDecays:limitTau0", Flag::ParticleDecays_limitTau0 },
  { "ParticleDecays:mixB", Flag::ParticleDecays_mixB },
  { "PartonLevel:all", Flag::PartonLevel_all },
  { "PartonLevel:earlyResDec", Flag::PartonLevel_earlyResDec },
  { "PartonLevel:FSR", Flag::PartonLevel_FSR },
  { "PartonLevel:FSRinProcess", Flag::PartonLevel_FSRinProcess },
  { "PartonLevel:FSRinResonances", Flag::PartonLevel_FSRinResonances },
  { "PartonLevel:ISR", Flag::PartonLevel_ISR },
  { "PartonLevel:MPI", Flag::PartonLevel_MPI },
  { "PartonLevel:Remnants", Flag::PartonLevel_Remnants },
  { "PDF:extrapolate", Flag::PDF_extrapolate },
  { "PDF:lepton", Flag::PDF_lepton },
  { "PDF:useHard", Flag::PDF_useHard },
  { "PhaseSpace:bias2Selection", Flag::PhaseSpace_bias2Selection },
  { "PhaseSpace:increaseMaximum", Flag::PhaseSpace_increaseMaximum },
  { "PhaseSpace:sameForSecond", Flag::PhaseSpace_sameForSecond },
  { "PhaseSpace:showSearch", Flag::PhaseSpace_showSearch },
  { "PhaseSpace:showViolation", Flag::PhaseSpace_showViolation },
  { "PhaseSpace:useBreitWigners", Flag::PhaseSpace_useBreitWigners },
  { "PhotonCollision:all", Flag::PhotonCollision_all },
  { "PhotonCollision:gmgm2bbbar", Flag::PhotonCollision_gmgm2bbbar },
  { "PhotonCollision:gmgm2ccbar", Flag::PhotonCollision_gmgm2ccbar },
  { "PhotonCollision:gmgm2ee", Flag::PhotonCollision_gmgm2ee },
  { "PhotonCollision:gmgm2mumu", Flag::PhotonCollision_gmgm2mumu },
  { "PhotonCollision:gmgm2qqbar", Flag::PhotonCollision_gmgm2qqbar },
  { "PhotonCollision:gmgm2tautau", Flag::PhotonCollision_gmgm2tautau },
  { "POWHEG:pythiaRandom", Flag::POWHEG_pythiaRandom },
  { "Print:quiet", Flag::Print_quiet },
  { "ProcessLevel:all", Flag::ProcessLevel_all },
  { "ProcessLevel:resonanceDecays", Flag::ProcessLevel_resonanceDecays },
  { "PromptPhoton:all", Flag::PromptPhoton_all },
  { "PromptPhoton:ffbar2gammagamma", Flag::PromptPhoton_ffbar2gammagamma },
  { "PromptPhoton:gg2gammagamma", Flag::PromptPhoton_gg2gammagamma },
  { "PromptPhoton:gg2ggamma", Flag::PromptPhoton_gg2ggamma },
  { "PromptPhoton:qg2qgamma", Flag::PromptPhoton_qg2qgamma },
  { "PromptPhoton:qqbar2ggamma", Flag::PromptPhoton_qqbar2ggamma },
  { "Random:setSeed", Flag::Random_setSeed },
  { "RHadrons:allow", Flag::RHadrons_allow },
  { "RHadrons:allowDecay", Flag::RHadrons_allowDecay },
  { "RHadrons:setMasses", Flag::RHadrons_setMasses },
  { "SecondHard:Bottomonium", Flag::SecondHard_Bottomonium },
  { "SecondHard:Charmonium", Flag::SecondHard_Charmonium },
  { "SecondHard:generate", Flag::SecondHard_generate },
  { "SecondHard:GmZAndJet", Flag::SecondHard_GmZAndJet },
  { "SecondHard:PhotonAndJet", Flag::SecondHard_PhotonAndJet },
  { "SecondHard:SingleGmZ", Flag::SecondHard_SingleGmZ },
  { "SecondHard:SingleTop", Flag::SecondHard_SingleTop },
  { "SecondHard:SingleW", Flag::SecondHard_SingleW },
  { "SecondHard:TopPair", Flag::SecondHard_TopPair },
  { "SecondHard:TwoBJets", Flag::SecondHard_TwoBJets },
  { "SecondHard:TwoJets", Flag::SecondHard_TwoJets },
  { "SecondHard:TwoPhotons", Flag::SecondHard_TwoPhotons },
  { "SecondHard:WAndJet", Flag::SecondHard_WAndJet },
  { "SigmaDiffractive:dampen", Flag::SigmaDiffractive_dampen },
  { "SigmaElastic:setOwn", Flag::SigmaElastic_setOwn },
  { "SigmaProcess:bMassiveME", Flag::SigmaProcess_bMassiveME },
  { "SigmaProcess:cMassiveME", Flag::SigmaProcess_cMassiveME },
  { "SigmaProcess:muMassiveME", Flag::SigmaProcess_muMassiveME },
  { "SigmaProcess:tauMassiveME", Flag::SigmaProcess_tauMassiveME },
  { "SigmaTotal:setOwn", Flag::SigmaTotal_setOwn },
  { "SigmaTotal:zeroAXB", Flag::SigmaTotal_zeroAXB },
  { "SLHA:allowUserOverride", Flag::SLHA_allowUserOverride },
  { "SLHA:keepSM", Flag::SLHA_keepSM },
  { "SLHA:NMSSM", Flag::SLHA_NMSSM },
  { "SLHA:useDecayTable", Flag::SLHA_useDecayTable },
  { "SoftQCD:all", Flag::SoftQCD_all },
  { "SoftQCD:centralDiffractive", Flag::SoftQCD_centralDiffractive },
  { "SoftQCD:doubleDiffractive", Flag::SoftQCD_doubleDiffractive },
  { "SoftQCD:elastic", Flag::SoftQCD_elastic },
  { "SoftQCD:inelastic", Flag::SoftQCD_inelastic },
  { "SoftQCD:nonDiffractive", Flag::SoftQCD_nonDiffractive },
  { "SoftQCD:singleDiffractive", Flag::SoftQCD_singleDiffractive },
  { "SpaceShower:alphaSuseCMW", Flag::SpaceShower_alphaSuseCMW },
  { "SpaceShower:MEafterFirst", Flag::SpaceShower_MEafterFirst },
  { "SpaceShower:MEcorrections", Flag::SpaceShower_MEcorrections },
  { "SpaceShower:phiIntAsym", Flag::SpaceShower_phiIntAsym },
  { "SpaceShower:phiPolAsym", Flag::SpaceShower_phiPolAsym },
  { "SpaceShower:QCDshower", Flag::SpaceShower_QCDshower },
  { "SpaceShower:QEDshowerByL", Flag::SpaceShower_QEDshowerByL },
  { "SpaceShower:QEDshowerByQ", Flag::SpaceShower_QEDshowerByQ },
  { "SpaceShower:rapidityOrder", Flag::SpaceShower_rapidityOrder },
  { "SpaceShower:samePTasMPI", Flag::SpaceShower_samePTasMPI },
  { "SpaceShower:useFixedFacScale", Flag::SpaceShower_useFixedFacScale },
  { "SpaceShower:weakShower", Flag::SpaceShower_weakShower },
  { "Stat:reset", Flag::Stat_reset },
  { "Stat:showErrors", Flag::Stat_showErrors },
  { "Stat:showPartonLevel", Flag::Stat_showPartonLevel },
  { "Stat:showProcessLevel", Flag::Stat_showProcessLevel },
  { "StringFlav:suppressLeadingB", Flag::StringFlav_suppressLeadingB },
  { "StringZ:useNonstandardB", Flag::StringZ_useNonstandardB },
  { "StringZ:useNonstandardC", Flag::StringZ_useNonstandardC },
  { "StringZ:useNonstandardH", Flag::StringZ_useNonstandardH },
  { "StringZ:usePetersonB", Flag::StringZ_usePetersonB },
  { "StringZ:usePetersonC", Flag::StringZ_usePetersonC },
  { "StringZ:usePetersonH", Flag::StringZ_usePetersonH },
  { "SUSY:all", Flag::SUSY_all },
  { "SUSY:gg2gluinogluino", Flag::SUSY_gg2gluinogluino },
  { "SUSY:gg2squarkantisquark", Flag::SUSY_gg2squarkantisquark },
  { "SUSY:qg2chi+-squark", Flag::SUSY_qg2chiplusminussquark },
  { "SUSY:qg2chi0squark", Flag::SUSY_qg2chi0squark },
  { "SUSY:qg2squarkgluino", Flag::SUSY_qg2squarkgluino },
  { "SUSY:qq2antisquark", Flag::SUSY_qq2antisquark },
  { "SUSY:qq2squarksquark", Flag::SUSY_qq2squarksquark },
  { "SUSY:qq2squarksquark:onlyQCD", Flag::SUSY_qq2squarksquark_onlyQCD },
  { "SUSY:qqbar2chi+-chi0", Flag::SUSY_qqbar2chiplusminuschi0 },
  { "SUSY:qqbar2chi+-gluino", Flag::SUSY_qqbar2chiplusminusgluino },
  { "SUSY:qqbar2chi+chi-", Flag::SUSY_qqbar2chipluschiminus },
  { "SUSY:qqbar2chi0chi0", Flag::SUSY_qqbar2chi0chi0 },
  { "SUSY:qqbar2chi0gluino", Flag::SUSY_qqbar2chi0gluino },
  { "SUSY:qqbar2gluinogluino", Flag::SUSY_qqbar2gluinogluino },
  { "SUSY:qqbar2sleptonantislepton", Flag::SUSY_qqbar2sleptonantislepton },
  { "SUSY:qqbar2squarkantisquark", Flag::SUSY_qqbar2squarkantisquark },
  { "SUSY:qqbar2squarkantisquark:onlyQCD", Flag::SUSY_qqbar2squarkantisquark_onlyQCD },
  { "SUSYResonance:3BodyMatrixElement", Flag::SUSYResonance_3BodyMatrixElement },
  { "TimeShower:allowBeamRecoil", Flag::TimeShower_allowBeamRecoil },
  { "TimeShower:allowMPIdipole", Flag::TimeShower_allowMPIdipole },
  { "TimeShower:alphaSuseCMW", Flag::TimeShower_alphaSuseCMW },
  { "TimeShower:dampenBeamRecoil", Flag::TimeShower_dampenBeamRecoil },
  { "TimeShower:globalRecoil", Flag::TimeShower_globalRecoil },
  { "TimeShower:interleave", Flag::TimeShower_interleave },
  { "TimeShower:limitPTmaxGlobal", Flag::TimeShower_limitPTmaxGlobal },
  { "TimeShower:MEafterFirst", Flag::TimeShower_MEafterFirst },
  { "TimeShower:MEcorrections", Flag::TimeShower_MEcorrections },
  { "TimeShower:phiPolAsym", Flag::TimeShower_phiPolAsym },
  { "TimeShower:QCDshower", Flag::TimeShower_QCDshower },
  { "TimeShower:QEDshowerByGamma", Flag::TimeShower_QEDshowerByGamma },
  { "TimeShower:QEDshowerByL", Flag::TimeShower_QEDshowerByL },
  { "TimeShower:QEDshowerByQ", Flag::TimeShower_QEDshowerByQ },
  { "TimeShower:recoilToColoured", Flag::TimeShower_recoilToColoured },
  { "TimeShower:useFixedFacScale", Flag::TimeShower_useFixedFacScale },
  { "TimeShower:weakShower", Flag::TimeShower_weakShower },
  { "Top:all", Flag::Top_all },
  { "Top:ffbar2tqbar(s:W)", Flag::Top_ffbar2tqbar_s_W_ },
  { "Top:ffbar2ttbar(s:gmZ)", Flag::Top_ffbar2ttbar_s_gmZ_ },
  { "Top:gg2ttbar", Flag::Top_gg2ttbar },
  { "Top:gmgm2ttbar", Flag::Top_gmgm2ttbar },
  { "Top:qq2tq(t:W)", Flag::Top_qq2tq_t_W_ },
  { "Top:qqbar2ttbar", Flag::Top_qqbar2ttbar },
  { "WeakBosonAndParton:all", Flag::WeakBosonAndParton_all },
  { "WeakBosonAndParton:ffbar2gmZgm", Flag::WeakBosonAndParton_ffbar2gmZgm },
  { "WeakBosonAndParton:ffbar2Wgm", Flag::WeakBosonAndParton_ffbar2Wgm },
  { "WeakBosonAndParton:fgm2gmZf", Flag::WeakBosonAndParton_fgm2gmZf },
  { "WeakBosonAndParton:fgm2Wf", Flag::WeakBosonAndParton_fgm2Wf },
  { "WeakBosonAndParton:qg2gmZq", Flag::WeakBosonAndParton_qg2gmZq },
  { "WeakBosonAndParton:qg2Wq", Flag::WeakBosonAndParton_qg2Wq },
  { "WeakBosonAndParton:qqbar2gmZg", Flag::WeakBosonAndParton_qqbar2gmZg },
  { "WeakBosonAndParton:qqbar2Wg", Flag::WeakBosonAndParton_qqbar2Wg },
  { "WeakBosonExchange:all", Flag::WeakBosonExchange_all },
  { "WeakBosonExchange:ff2ff(t:gmZ)", Flag::WeakBosonExchange_ff2ff_t_gmZ_ },
  { "WeakBosonExchange:ff2ff(t:W)", Flag::WeakBosonExchange_ff2ff_t_W_ },
  { "WeakDoubleBoson:all", Flag::WeakDoubleBoson_all },
  { "WeakDoubleBoson:ffbar2gmZgmZ", Flag::WeakDoubleBoson_ffbar2gmZgmZ },
  { "WeakDoubleBoson:ffbar2WW", Flag::WeakDoubleBoson_ffbar2WW },
  { "WeakDoubleBoson:ffbar2ZW", Flag::WeakDoubleBoson_ffbar2ZW },
  { "WeakShower:singleEmission", Flag::WeakShower_singleEmission },
  { "WeakShower:vetoQCDjets", Flag::WeakShower_vetoQCDjets },
  { "WeakShower:vetoWeakJets", Flag::WeakShower_vetoWeakJets },
  { "WeakSingleBoson:all", Flag::WeakSingleBoson_all },
  { "WeakSingleBoson:ffbar2ffbar(s:gm)", Flag::WeakSingleBoson_ffbar2ffbar_s_gm_ },
  { "WeakSingleBoson:ffbar2ffbar(s:gmZ)", Flag::WeakSingleBoson_ffbar2ffbar_s_gmZ_ },
  { "WeakSingleBoson:ffbar2ffbar(s:W)", Flag::WeakSingleBoson_ffbar2ffbar_s_W_ },
  { "WeakSingleBoson:ffbar2gmZ", Flag::WeakSingleBoson_ffbar2gmZ },
  { "WeakSingleBoson:ffbar2W", Flag::WeakSingleBoson_ffbar2W },
  { "Zprime:coup2gen4", Flag::Zprime_coup2gen4 },
  { "Zprime:universality", Flag::Zprime_universality },
};

const map<string, Mode> ModeMap = 
{
  { "BeamRemnants:companionPower", Mode::BeamRemnants_companionPower },
  { "BeamRemnants:maxValQuark", Mode::BeamRemnants_maxValQuark },
  { "BeamRemnants:remnantMode", Mode::BeamRemnants_remnantMode },
  { "Beams:frameType", Mode::Beams_frameType },
  { "Beams:idA", Mode::Beams_idA },
  { "Beams:idB", Mode::Beams_idB },
  { "Beams:nSkipLHEFatInit", Mode::Beams_nSkipLHEFatInit },
  { "Check:levelParticleData", Mode::Check_levelParticleData },
  { "Check:nErrList", Mode::Check_nErrList },
  { "ColourReconnection:flipMode", Mode::ColourReconnection_flipMode },
  { "ColourReconnection:lambdaForm", Mode::ColourReconnection_lambdaForm },
  { "ColourReconnection:mode", Mode::ColourReconnection_mode },
  { "ColourReconnection:nColours", Mode::ColourReconnection_nColours },
  { "ColourReconnection:timeDilationMode", Mode::ColourReconnection_timeDilationMode },
  { "ContactInteractions:etaLL", Mode::ContactInteractions_etaLL },
  { "ContactInteractions:etaLR", Mode::ContactInteractions_etaLR },
  { "ContactInteractions:etaRR", Mode::ContactInteractions_etaRR },
  { "ContactInteractions:nQuarkNew", Mode::ContactInteractions_nQuarkNew },
  { "Diffraction:bProfile", Mode::Diffraction_bProfile },
  { "Diffraction:PomFlux", Mode::Diffraction_PomFlux },
  { "Diffraction:sampleType", Mode::Diffraction_sampleType },
  { "Event:startColTag", Mode::Event_startColTag },
  { "ExtraDimensionsG*:KKintMode", Mode::ExtraDimensionsGstar_KKintMode },
  { "ExtraDimensionsLED:CutOffMode", Mode::ExtraDimensionsLED_CutOffMode },
  { "ExtraDimensionsLED:n", Mode::ExtraDimensionsLED_n },
  { "ExtraDimensionsLED:NegInt", Mode::ExtraDimensionsLED_NegInt },
  { "ExtraDimensionsLED:nQuarkNew", Mode::ExtraDimensionsLED_nQuarkNew },
  { "ExtraDimensionsLED:opMode", Mode::ExtraDimensionsLED_opMode },
  { "ExtraDimensionsTEV:gmZmode", Mode::ExtraDimensionsTEV_gmZmode },
  { "ExtraDimensionsTEV:nMax", Mode::ExtraDimensionsTEV_nMax },
  { "ExtraDimensionsUnpart:CutOffMode", Mode::ExtraDimensionsUnpart_CutOffMode },
  { "ExtraDimensionsUnpart:gXX", Mode::ExtraDimensionsUnpart_gXX },
  { "ExtraDimensionsUnpart:gXY", Mode::ExtraDimensionsUnpart_gXY },
  { "ExtraDimensionsUnpart:spinU", Mode::ExtraDimensionsUnpart_spinU },
  { "HadronScatter:hadronSelect", Mode::HadronScatter_hadronSelect },
  { "HadronScatter:scatterProb", Mode::HadronScatter_scatterProb },
  { "HardQCD:nQuarkNew", Mode::HardQCD_nQuarkNew },
  { "HiddenValley:nFlav", Mode::HiddenValley_nFlav },
  { "HiddenValley:Ngauge", Mode::HiddenValley_Ngauge },
  { "HiddenValley:spinFv", Mode::HiddenValley_spinFv },
  { "HiddenValley:spinqv", Mode::HiddenValley_spinqv },
  { "HiggsA3:parity", Mode::HiggsA3_parity },
  { "HiggsH1:parity", Mode::HiggsH1_parity },
  { "HiggsH2:parity", Mode::HiggsH2_parity },
  { "Init:showOneParticleData", Mode::Init_showOneParticleData },
  { "JetMatching:exclusive", Mode::JetMatching_exclusive },
  { "JetMatching:jetAlgorithm", Mode::JetMatching_jetAlgorithm },
  { "JetMatching:jetAllow", Mode::JetMatching_jetAllow },
  { "JetMatching:jetMatch", Mode::JetMatching_jetMatch },
  { "JetMatching:nEta", Mode::JetMatching_nEta },
  { "JetMatching:nJet", Mode::JetMatching_nJet },
  { "JetMatching:nJetMax", Mode::JetMatching_nJetMax },
  { "JetMatching:nPartonsNow", Mode::JetMatching_nPartonsNow },
  { "JetMatching:nPhi", Mode::JetMatching_nPhi },
  { "JetMatching:nQmatch", Mode::JetMatching_nQmatch },
  { "JetMatching:scheme", Mode::JetMatching_scheme },
  { "JetMatching:slowJetPower", Mode::JetMatching_slowJetPower },
  { "LesHouches:idRenameBeams", Mode::LesHouches_idRenameBeams },
  { "LesHouches:setLeptonMass", Mode::LesHouches_setLeptonMass },
  { "LesHouches:setLifetime", Mode::LesHouches_setLifetime },
  { "Main:numberOfEvents", Mode::Main_numberOfEvents },
  { "Main:numberOfSubruns", Mode::Main_numberOfSubruns },
  { "Main:spareMode1", Mode::Main_spareMode1 },
  { "Main:spareMode2", Mode::Main_spareMode2 },
  { "Main:spareMode3", Mode::Main_spareMode3 },
  { "Main:subrun", Mode::Main_subrun },
  { "Main:timesAllowErrors", Mode::Main_timesAllowErrors },
  { "Merging:incompleteScalePrescrip", Mode::Merging_incompleteScalePrescrip },
  { "Merging:ktType", Mode::Merging_ktType },
  { "Merging:nJetMax", Mode::Merging_nJetMax },
  { "Merging:nJetMaxNLO", Mode::Merging_nJetMaxNLO },
  { "Merging:nQuarksMerge", Mode::Merging_nQuarksMerge },
  { "Merging:nRecluster", Mode::Merging_nRecluster },
  { "Merging:nRequested", Mode::Merging_nRequested },
  { "Merging:unorderedASscalePrescrip", Mode::Merging_unorderedASscalePrescrip },
  { "Merging:unorderedPDFscalePrescrip", Mode::Merging_unorderedPDFscalePrescrip },
  { "Merging:unorderedScalePrescrip", Mode::Merging_unorderedScalePrescrip },
  { "MiniStringFragmentation:nTry", Mode::MiniStringFragmentation_nTry },
  { "MultipartonInteractions:alphaEMorder", Mode::MultipartonInteractions_alphaEMorder },
  { "MultipartonInteractions:alphaSorder", Mode::MultipartonInteractions_alphaSorder },
  { "MultipartonInteractions:bProfile", Mode::MultipartonInteractions_bProfile },
  { "MultipartonInteractions:bSelScale", Mode::MultipartonInteractions_bSelScale },
  { "MultipartonInteractions:enhanceScreening", Mode::MultipartonInteractions_enhanceScreening },
  { "MultipartonInteractions:nQuarkIn", Mode::MultipartonInteractions_nQuarkIn },
  { "MultipartonInteractions:nSample", Mode::MultipartonInteractions_nSample },
  { "MultipartonInteractions:processLevel", Mode::MultipartonInteractions_processLevel },
  { "MultipartonInteractions:pTmaxMatch", Mode::MultipartonInteractions_pTmaxMatch },
  { "MultipartonInteractions:rescatterMode", Mode::MultipartonInteractions_rescatterMode },
  { "Next:numberCount", Mode::Next_numberCount },
  { "Next:numberShowEvent", Mode::Next_numberShowEvent },
  { "Next:numberShowInfo", Mode::Next_numberShowInfo },
  { "Next:numberShowLHA", Mode::Next_numberShowLHA },
  { "Next:numberShowProcess", Mode::Next_numberShowProcess },
  { "ParticleData:modeBreitWigner", Mode::ParticleData_modeBreitWigner },
  { "PDF:PomSet", Mode::PDF_PomSet },
  { "PDFinProcess:nQuarkIn", Mode::PDFinProcess_nQuarkIn },
  { "POWHEG:emitted", Mode::POWHEG_emitted },
  { "POWHEG:MPIveto", Mode::POWHEG_MPIveto },
  { "POWHEG:nFinal", Mode::POWHEG_nFinal },
  { "POWHEG:pTdef", Mode::POWHEG_pTdef },
  { "POWHEG:pTemt", Mode::POWHEG_pTemt },
  { "POWHEG:pThard", Mode::POWHEG_pThard },
  { "POWHEG:QEDveto", Mode::POWHEG_QEDveto },
  { "POWHEG:veto", Mode::POWHEG_veto },
  { "POWHEG:vetoCount", Mode::POWHEG_vetoCount },
  { "PromptPhoton:nQuarkLoop", Mode::PromptPhoton_nQuarkLoop },
  { "Pythia:versionDate", Mode::Pythia_versionDate },
  { "Random:seed", Mode::Random_seed },
  { "RHadrons:idGluino", Mode::RHadrons_idGluino },
  { "RHadrons:idSbottom", Mode::RHadrons_idSbottom },
  { "RHadrons:idStop", Mode::RHadrons_idStop },
  { "SigmaProcess:alphaEMorder", Mode::SigmaProcess_alphaEMorder },
  { "SigmaProcess:alphaSorder", Mode::SigmaProcess_alphaSorder },
  { "SigmaProcess:factorScale1", Mode::SigmaProcess_factorScale1 },
  { "SigmaProcess:factorScale2", Mode::SigmaProcess_factorScale2 },
  { "SigmaProcess:factorScale3", Mode::SigmaProcess_factorScale3 },
  { "SigmaProcess:factorScale3VV", Mode::SigmaProcess_factorScale3VV },
  { "SigmaProcess:renormScale1", Mode::SigmaProcess_renormScale1 },
  { "SigmaProcess:renormScale2", Mode::SigmaProcess_renormScale2 },
  { "SigmaProcess:renormScale3", Mode::SigmaProcess_renormScale3 },
  { "SigmaProcess:renormScale3VV", Mode::SigmaProcess_renormScale3VV },
  { "SLHA:meMode", Mode::SLHA_meMode },
  { "SLHA:readFrom", Mode::SLHA_readFrom },
  { "SLHA:verbose", Mode::SLHA_verbose },
  { "SpaceShower:alphaEMorder", Mode::SpaceShower_alphaEMorder },
  { "SpaceShower:alphaSorder", Mode::SpaceShower_alphaSorder },
  { "SpaceShower:nQuarkIn", Mode::SpaceShower_nQuarkIn },
  { "SpaceShower:pTdampMatch", Mode::SpaceShower_pTdampMatch },
  { "SpaceShower:pTmaxMatch", Mode::SpaceShower_pTmaxMatch },
  { "SpaceShower:weakShowerMode", Mode::SpaceShower_weakShowerMode },
  { "StandardModel:alphaSnfmax", Mode::StandardModel_alphaSnfmax },
  { "SUSY:idA", Mode::SUSY_idA },
  { "SUSY:idB", Mode::SUSY_idB },
  { "SUSY:sin2thetaWMode", Mode::SUSY_sin2thetaWMode },
  { "TauDecays:externalMode", Mode::TauDecays_externalMode },
  { "TauDecays:mode", Mode::TauDecays_mode },
  { "TauDecays:tauMother", Mode::TauDecays_tauMother },
  { "TimeShower:alphaEMorder", Mode::TimeShower_alphaEMorder },
  { "TimeShower:alphaSorder", Mode::TimeShower_alphaSorder },
  { "TimeShower:globalRecoilMode", Mode::TimeShower_globalRecoilMode },
  { "TimeShower:nGammaToLepton", Mode::TimeShower_nGammaToLepton },
  { "TimeShower:nGammaToQuark", Mode::TimeShower_nGammaToQuark },
  { "TimeShower:nGluonToQuark", Mode::TimeShower_nGluonToQuark },
  { "TimeShower:nMaxGlobalBranch", Mode::TimeShower_nMaxGlobalBranch },
  { "TimeShower:nMaxGlobalRecoil", Mode::TimeShower_nMaxGlobalRecoil },
  { "TimeShower:nPartonsInBorn", Mode::TimeShower_nPartonsInBorn },
  { "TimeShower:pTdampMatch", Mode::TimeShower_pTdampMatch },
  { "TimeShower:pTmaxMatch", Mode::TimeShower_pTmaxMatch },
  { "TimeShower:weakShowerMode", Mode::TimeShower_weakShowerMode },
  { "TimeShower:weightGluonToQuark", Mode::TimeShower_weightGluonToQuark },
  { "Tune:ee", Mode::Tune_ee },
  { "Tune:pp", Mode::Tune_pp },
  { "Tune:preferLHAPDF", Mode::Tune_preferLHAPDF },
  { "WeakZ0:gmZmode", Mode::WeakZ0_gmZmode },
  { "Zprime:gmZmode", Mode::Zprime_gmZmode },
};

const map<string, Param> ParamMap = 
{
  { "BeamRemnants:gluonPower", Param::BeamRemnants_gluonPower },
  { "BeamRemnants:halfMassForKT", Param::BeamRemnants_halfMassForKT },
  { "BeamRemnants:halfScaleForKT", Param::BeamRemnants_halfScaleForKT },
  { "BeamRemnants:primordialKThard", Param::BeamRemnants_primordialKThard },
  { "BeamRemnants:primordialKTremnant", Param::BeamRemnants_primordialKTremnant },
  { "BeamRemnants:primordialKTsoft", Param::BeamRemnants_primordialKTsoft },
  { "BeamRemnants:reducedKTatHighY", Param::BeamRemnants_reducedKTatHighY },
  { "BeamRemnants:saturation", Param::BeamRemnants_saturation },
  { "BeamRemnants:valenceDiqEnhance", Param::BeamRemnants_valenceDiqEnhance },
  { "BeamRemnants:valencePowerDinP", Param::BeamRemnants_valencePowerDinP },
  { "BeamRemnants:valencePowerMeson", Param::BeamRemnants_valencePowerMeson },
  { "BeamRemnants:valencePowerUinP", Param::BeamRemnants_valencePowerUinP },
  { "BeamRemnants:xGluonCutoff", Param::BeamRemnants_xGluonCutoff },
  { "Beams:eA", Param::Beams_eA },
  { "Beams:eB", Param::Beams_eB },
  { "Beams:eCM", Param::Beams_eCM },
  { "Beams:maxDevA", Param::Beams_maxDevA },
  { "Beams:maxDevB", Param::Beams_maxDevB },
  { "Beams:maxDevTime", Param::Beams_maxDevTime },
  { "Beams:maxDevVertex", Param::Beams_maxDevVertex },
  { "Beams:offsetTime", Param::Beams_offsetTime },
  { "Beams:offsetVertexX", Param::Beams_offsetVertexX },
  { "Beams:offsetVertexY", Param::Beams_offsetVertexY },
  { "Beams:offsetVertexZ", Param::Beams_offsetVertexZ },
  { "Beams:pxA", Param::Beams_pxA },
  { "Beams:pxB", Param::Beams_pxB },
  { "Beams:pyA", Param::Beams_pyA },
  { "Beams:pyB", Param::Beams_pyB },
  { "Beams:pzA", Param::Beams_pzA },
  { "Beams:pzB", Param::Beams_pzB },
  { "Beams:sigmaPxA", Param::Beams_sigmaPxA },
  { "Beams:sigmaPxB", Param::Beams_sigmaPxB },
  { "Beams:sigmaPyA", Param::Beams_sigmaPyA },
  { "Beams:sigmaPyB", Param::Beams_sigmaPyB },
  { "Beams:sigmaPzA", Param::Beams_sigmaPzA },
  { "Beams:sigmaPzB", Param::Beams_sigmaPzB },
  { "Beams:sigmaTime", Param::Beams_sigmaTime },
  { "Beams:sigmaVertexX", Param::Beams_sigmaVertexX },
  { "Beams:sigmaVertexY", Param::Beams_sigmaVertexY },
  { "Beams:sigmaVertexZ", Param::Beams_sigmaVertexZ },
  { "BoseEinstein:lambda", Param::BoseEinstein_lambda },
  { "BoseEinstein:QRef", Param::BoseEinstein_QRef },
  { "BoseEinstein:widthSep", Param::BoseEinstein_widthSep },
  { "Check:epTolErr", Param::Check_epTolErr },
  { "Check:epTolWarn", Param::Check_epTolWarn },
  { "Check:mTolErr", Param::Check_mTolErr },
  { "Check:mTolWarn", Param::Check_mTolWarn },
  { "ColourReconnection:blowR", Param::ColourReconnection_blowR },
  { "ColourReconnection:blowT", Param::ColourReconnection_blowT },
  { "ColourReconnection:dLambdaCut", Param::ColourReconnection_dLambdaCut },
  { "ColourReconnection:fracGluon", Param::ColourReconnection_fracGluon },
  { "ColourReconnection:fragmentationTime", Param::ColourReconnection_fragmentationTime },
  { "ColourReconnection:junctionCorrection", Param::ColourReconnection_junctionCorrection },
  { "ColourReconnection:kI", Param::ColourReconnection_kI },
  { "ColourReconnection:m0", Param::ColourReconnection_m0 },
  { "ColourReconnection:m2Lambda", Param::ColourReconnection_m2Lambda },
  { "ColourReconnection:range", Param::ColourReconnection_range },
  { "ColourReconnection:rHadron", Param::ColourReconnection_rHadron },
  { "ColourReconnection:timeDilationPar", Param::ColourReconnection_timeDilationPar },
  { "ContactInteractions:Lambda", Param::ContactInteractions_Lambda },
  { "Diffraction:coreFraction", Param::Diffraction_coreFraction },
  { "Diffraction:coreRadius", Param::Diffraction_coreRadius },
  { "Diffraction:expPow", Param::Diffraction_expPow },
  { "Diffraction:largeMassSuppress", Param::Diffraction_largeMassSuppress },
  { "Diffraction:MBRalpha", Param::Diffraction_MBRalpha },
  { "Diffraction:MBRbeta0", Param::Diffraction_MBRbeta0 },
  { "Diffraction:MBRdyminCD", Param::Diffraction_MBRdyminCD },
  { "Diffraction:MBRdyminCDflux", Param::Diffraction_MBRdyminCDflux },
  { "Diffraction:MBRdyminDD", Param::Diffraction_MBRdyminDD },
  { "Diffraction:MBRdyminDDflux", Param::Diffraction_MBRdyminDDflux },
  { "Diffraction:MBRdyminSD", Param::Diffraction_MBRdyminSD },
  { "Diffraction:MBRdyminSDflux", Param::Diffraction_MBRdyminSDflux },
  { "Diffraction:MBRdyminSigCD", Param::Diffraction_MBRdyminSigCD },
  { "Diffraction:MBRdyminSigDD", Param::Diffraction_MBRdyminSigDD },
  { "Diffraction:MBRdyminSigSD", Param::Diffraction_MBRdyminSigSD },
  { "Diffraction:MBRepsilon", Param::Diffraction_MBRepsilon },
  { "Diffraction:MBRm2Min", Param::Diffraction_MBRm2Min },
  { "Diffraction:MBRsigma0", Param::Diffraction_MBRsigma0 },
  { "Diffraction:mMinPert", Param::Diffraction_mMinPert },
  { "Diffraction:mPowPomP", Param::Diffraction_mPowPomP },
  { "Diffraction:mRefPomP", Param::Diffraction_mRefPomP },
  { "Diffraction:mWidthPert", Param::Diffraction_mWidthPert },
  { "Diffraction:pickQuarkNorm", Param::Diffraction_pickQuarkNorm },
  { "Diffraction:pickQuarkPower", Param::Diffraction_pickQuarkPower },
  { "Diffraction:PomFluxAlphaPrime", Param::Diffraction_PomFluxAlphaPrime },
  { "Diffraction:PomFluxEpsilon", Param::Diffraction_PomFluxEpsilon },
  { "Diffraction:primKTwidth", Param::Diffraction_primKTwidth },
  { "Diffraction:probMaxPert", Param::Diffraction_probMaxPert },
  { "Diffraction:sigmaRefPomP", Param::Diffraction_sigmaRefPomP },
  { "ExcitedFermion:contactDec", Param::ExcitedFermion_contactDec },
  { "ExcitedFermion:coupF", Param::ExcitedFermion_coupF },
  { "ExcitedFermion:coupFcol", Param::ExcitedFermion_coupFcol },
  { "ExcitedFermion:coupFprime", Param::ExcitedFermion_coupFprime },
  { "ExcitedFermion:Lambda", Param::ExcitedFermion_Lambda },
  { "ExtraDimensionsG*:Gbb", Param::ExtraDimensionsGstar_Gbb },
  { "ExtraDimensionsG*:Ggg", Param::ExtraDimensionsGstar_Ggg },
  { "ExtraDimensionsG*:Ggmgm", Param::ExtraDimensionsGstar_Ggmgm },
  { "ExtraDimensionsG*:Ghh", Param::ExtraDimensionsGstar_Ghh },
  { "ExtraDimensionsG*:Gll", Param::ExtraDimensionsGstar_Gll },
  { "ExtraDimensionsG*:Gqq", Param::ExtraDimensionsGstar_Gqq },
  { "ExtraDimensionsG*:Gtt", Param::ExtraDimensionsGstar_Gtt },
  { "ExtraDimensionsG*:GWW", Param::ExtraDimensionsGstar_GWW },
  { "ExtraDimensionsG*:GZZ", Param::ExtraDimensionsGstar_GZZ },
  { "ExtraDimensionsG*:kappaMG", Param::ExtraDimensionsGstar_kappaMG },
  { "ExtraDimensionsG*:KKgbL", Param::ExtraDimensionsGstar_KKgbL },
  { "ExtraDimensionsG*:KKgbR", Param::ExtraDimensionsGstar_KKgbR },
  { "ExtraDimensionsG*:KKgqL", Param::ExtraDimensionsGstar_KKgqL },
  { "ExtraDimensionsG*:KKgqR", Param::ExtraDimensionsGstar_KKgqR },
  { "ExtraDimensionsG*:KKgtL", Param::ExtraDimensionsGstar_KKgtL },
  { "ExtraDimensionsG*:KKgtR", Param::ExtraDimensionsGstar_KKgtR },
  { "ExtraDimensionsLED:c", Param::ExtraDimensionsLED_c },
  { "ExtraDimensionsLED:g", Param::ExtraDimensionsLED_g },
  { "ExtraDimensionsLED:LambdaT", Param::ExtraDimensionsLED_LambdaT },
  { "ExtraDimensionsLED:MD", Param::ExtraDimensionsLED_MD },
  { "ExtraDimensionsLED:t", Param::ExtraDimensionsLED_t },
  { "ExtraDimensionsTEV:mStar", Param::ExtraDimensionsTEV_mStar },
  { "ExtraDimensionsUnpart:dU", Param::ExtraDimensionsUnpart_dU },
  { "ExtraDimensionsUnpart:lambda", Param::ExtraDimensionsUnpart_lambda },
  { "ExtraDimensionsUnpart:LambdaU", Param::ExtraDimensionsUnpart_LambdaU },
  { "ExtraDimensionsUnpart:ratio", Param::ExtraDimensionsUnpart_ratio },
  { "FourthGeneration:VcbPrime", Param::FourthGeneration_VcbPrime },
  { "FourthGeneration:VtbPrime", Param::FourthGeneration_VtbPrime },
  { "FourthGeneration:VtPrimeb", Param::FourthGeneration_VtPrimeb },
  { "FourthGeneration:VtPrimebPrime", Param::FourthGeneration_VtPrimebPrime },
  { "FourthGeneration:VtPrimed", Param::FourthGeneration_VtPrimed },
  { "FourthGeneration:VtPrimes", Param::FourthGeneration_VtPrimes },
  { "FourthGeneration:VubPrime", Param::FourthGeneration_VubPrime },
  { "FragmentationSystems:mJoin", Param::FragmentationSystems_mJoin },
  { "FragmentationSystems:mJoinJunction", Param::FragmentationSystems_mJoinJunction },
  { "HadronLevel:mStringMin", Param::HadronLevel_mStringMin },
  { "HadronScatter:j", Param::HadronScatter_j },
  { "HadronScatter:k", Param::HadronScatter_k },
  { "HadronScatter:N", Param::HadronScatter_N },
  { "HadronScatter:p", Param::HadronScatter_p },
  { "HadronScatter:rMax", Param::HadronScatter_rMax },
  { "HiddenValley:alphaFSR", Param::HiddenValley_alphaFSR },
  { "HiddenValley:aLund", Param::HiddenValley_aLund },
  { "HiddenValley:bmqv2", Param::HiddenValley_bmqv2 },
  { "HiddenValley:kappa", Param::HiddenValley_kappa },
  { "HiddenValley:kinMix", Param::HiddenValley_kinMix },
  { "HiddenValley:probVector", Param::HiddenValley_probVector },
  { "HiddenValley:pTminFSR", Param::HiddenValley_pTminFSR },
  { "HiddenValley:rFactqv", Param::HiddenValley_rFactqv },
  { "HiddenValley:sigmamqv", Param::HiddenValley_sigmamqv },
  { "Higgs:wingsFac", Param::Higgs_wingsFac },
  { "HiggsA3:coup2d", Param::HiggsA3_coup2d },
  { "HiggsA3:coup2H1H1", Param::HiggsA3_coup2H1H1 },
  { "HiggsA3:coup2H1Z", Param::HiggsA3_coup2H1Z },
  { "HiggsA3:coup2H2Z", Param::HiggsA3_coup2H2Z },
  { "HiggsA3:coup2Hchg", Param::HiggsA3_coup2Hchg },
  { "HiggsA3:coup2HchgW", Param::HiggsA3_coup2HchgW },
  { "HiggsA3:coup2l", Param::HiggsA3_coup2l },
  { "HiggsA3:coup2u", Param::HiggsA3_coup2u },
  { "HiggsA3:coup2W", Param::HiggsA3_coup2W },
  { "HiggsA3:coup2Z", Param::HiggsA3_coup2Z },
  { "HiggsA3:etaParity", Param::HiggsA3_etaParity },
  { "HiggsA3:phiParity", Param::HiggsA3_phiParity },
  { "HiggsH1:coup2d", Param::HiggsH1_coup2d },
  { "HiggsH1:coup2Hchg", Param::HiggsH1_coup2Hchg },
  { "HiggsH1:coup2l", Param::HiggsH1_coup2l },
  { "HiggsH1:coup2u", Param::HiggsH1_coup2u },
  { "HiggsH1:coup2W", Param::HiggsH1_coup2W },
  { "HiggsH1:coup2Z", Param::HiggsH1_coup2Z },
  { "HiggsH1:etaParity", Param::HiggsH1_etaParity },
  { "HiggsH1:phiParity", Param::HiggsH1_phiParity },
  { "HiggsH2:coup2A3A3", Param::HiggsH2_coup2A3A3 },
  { "HiggsH2:coup2A3H1", Param::HiggsH2_coup2A3H1 },
  { "HiggsH2:coup2d", Param::HiggsH2_coup2d },
  { "HiggsH2:coup2H1H1", Param::HiggsH2_coup2H1H1 },
  { "HiggsH2:coup2H1Z", Param::HiggsH2_coup2H1Z },
  { "HiggsH2:coup2Hchg", Param::HiggsH2_coup2Hchg },
  { "HiggsH2:coup2HchgW", Param::HiggsH2_coup2HchgW },
  { "HiggsH2:coup2l", Param::HiggsH2_coup2l },
  { "HiggsH2:coup2u", Param::HiggsH2_coup2u },
  { "HiggsH2:coup2W", Param::HiggsH2_coup2W },
  { "HiggsH2:coup2Z", Param::HiggsH2_coup2Z },
  { "HiggsH2:etaParity", Param::HiggsH2_etaParity },
  { "HiggsH2:phiParity", Param::HiggsH2_phiParity },
  { "HiggsHchg:coup2H1W", Param::HiggsHchg_coup2H1W },
  { "HiggsHchg:coup2H2W", Param::HiggsHchg_coup2H2W },
  { "HiggsHchg:tanBeta", Param::HiggsHchg_tanBeta },
  { "JetMatching:clFact", Param::JetMatching_clFact },
  { "JetMatching:coneMatchHeavy", Param::JetMatching_coneMatchHeavy },
  { "JetMatching:coneMatchLight", Param::JetMatching_coneMatchLight },
  { "JetMatching:coneRadius", Param::JetMatching_coneRadius },
  { "JetMatching:coneRadiusHeavy", Param::JetMatching_coneRadiusHeavy },
  { "JetMatching:etaJetMax", Param::JetMatching_etaJetMax },
  { "JetMatching:eTjetMin", Param::JetMatching_eTjetMin },
  { "JetMatching:eTseed", Param::JetMatching_eTseed },
  { "JetMatching:eTthreshold", Param::JetMatching_eTthreshold },
  { "JetMatching:qCut", Param::JetMatching_qCut },
  { "JetMatching:qCutME", Param::JetMatching_qCutME },
  { "LeftRightSymmmetry:coupHee", Param::LeftRightSymmmetry_coupHee },
  { "LeftRightSymmmetry:coupHmue", Param::LeftRightSymmmetry_coupHmue },
  { "LeftRightSymmmetry:coupHmumu", Param::LeftRightSymmmetry_coupHmumu },
  { "LeftRightSymmmetry:coupHtaue", Param::LeftRightSymmmetry_coupHtaue },
  { "LeftRightSymmmetry:coupHtaumu", Param::LeftRightSymmmetry_coupHtaumu },
  { "LeftRightSymmmetry:coupHtautau", Param::LeftRightSymmmetry_coupHtautau },
  { "LeftRightSymmmetry:gL", Param::LeftRightSymmmetry_gL },
  { "LeftRightSymmmetry:gR", Param::LeftRightSymmmetry_gR },
  { "LeftRightSymmmetry:vL", Param::LeftRightSymmmetry_vL },
  { "LeptoQuark:kCoup", Param::LeptoQuark_kCoup },
  { "LesHouches:mRecalculate", Param::LesHouches_mRecalculate },
  { "Main:spareParm1", Param::Main_spareParm1 },
  { "Main:spareParm2", Param::Main_spareParm2 },
  { "Main:spareParm3", Param::Main_spareParm3 },
  { "Merging:aCollFSR", Param::Merging_aCollFSR },
  { "Merging:aCollISR", Param::Merging_aCollISR },
  { "Merging:Dparameter", Param::Merging_Dparameter },
  { "Merging:dRijMS", Param::Merging_dRijMS },
  { "Merging:fsrInRecNorm", Param::Merging_fsrInRecNorm },
  { "Merging:kFactor0j", Param::Merging_kFactor0j },
  { "Merging:kFactor1j", Param::Merging_kFactor1j },
  { "Merging:kFactor2j", Param::Merging_kFactor2j },
  { "Merging:muFac", Param::Merging_muFac },
  { "Merging:muFacInME", Param::Merging_muFacInME },
  { "Merging:muRen", Param::Merging_muRen },
  { "Merging:muRenInME", Param::Merging_muRenInME },
  { "Merging:nonJoinedNorm", Param::Merging_nonJoinedNorm },
  { "Merging:pTiMS", Param::Merging_pTiMS },
  { "Merging:QijMS", Param::Merging_QijMS },
  { "Merging:scaleSeparationFactor", Param::Merging_scaleSeparationFactor },
  { "Merging:TMS", Param::Merging_TMS },
  { "MultipartonInteractions:a1", Param::MultipartonInteractions_a1 },
  { "MultipartonInteractions:alphaSvalue", Param::MultipartonInteractions_alphaSvalue },
  { "MultipartonInteractions:coreFraction", Param::MultipartonInteractions_coreFraction },
  { "MultipartonInteractions:coreRadius", Param::MultipartonInteractions_coreRadius },
  { "MultipartonInteractions:deltaYRescatter", Param::MultipartonInteractions_deltaYRescatter },
  { "MultipartonInteractions:ecmPow", Param::MultipartonInteractions_ecmPow },
  { "MultipartonInteractions:ecmRef", Param::MultipartonInteractions_ecmRef },
  { "MultipartonInteractions:expPow", Param::MultipartonInteractions_expPow },
  { "MultipartonInteractions:Kfactor", Param::MultipartonInteractions_Kfactor },
  { "MultipartonInteractions:pT0Ref", Param::MultipartonInteractions_pT0Ref },
  { "MultipartonInteractions:pTmin", Param::MultipartonInteractions_pTmin },
  { "MultipartonInteractions:ySepRescatter", Param::MultipartonInteractions_ySepRescatter },
  { "Onia:massSplit", Param::Onia_massSplit },
  { "ParticleData:alphaSvalueMRun", Param::ParticleData_alphaSvalueMRun },
  { "ParticleData:maxEnhanceBW", Param::ParticleData_maxEnhanceBW },
  { "ParticleData:mbRun", Param::ParticleData_mbRun },
  { "ParticleData:mcRun", Param::ParticleData_mcRun },
  { "ParticleData:mdRun", Param::ParticleData_mdRun },
  { "ParticleData:msRun", Param::ParticleData_msRun },
  { "ParticleData:mtRun", Param::ParticleData_mtRun },
  { "ParticleData:muRun", Param::ParticleData_muRun },
  { "ParticleDecays:colRearrange", Param::ParticleDecays_colRearrange },
  { "ParticleDecays:mSafety", Param::ParticleDecays_mSafety },
  { "ParticleDecays:multGoffset", Param::ParticleDecays_multGoffset },
  { "ParticleDecays:multIncrease", Param::ParticleDecays_multIncrease },
  { "ParticleDecays:multIncreaseWeak", Param::ParticleDecays_multIncreaseWeak },
  { "ParticleDecays:multRefMass", Param::ParticleDecays_multRefMass },
  { "ParticleDecays:rMax", Param::ParticleDecays_rMax },
  { "ParticleDecays:sigmaSoft", Param::ParticleDecays_sigmaSoft },
  { "ParticleDecays:tau0Max", Param::ParticleDecays_tau0Max },
  { "ParticleDecays:tauMax", Param::ParticleDecays_tauMax },
  { "ParticleDecays:xBdMix", Param::ParticleDecays_xBdMix },
  { "ParticleDecays:xBsMix", Param::ParticleDecays_xBsMix },
  { "ParticleDecays:xyMax", Param::ParticleDecays_xyMax },
  { "ParticleDecays:zMax", Param::ParticleDecays_zMax },
  { "PDF:PomGluonA", Param::PDF_PomGluonA },
  { "PDF:PomGluonB", Param::PDF_PomGluonB },
  { "PDF:PomQuarkA", Param::PDF_PomQuarkA },
  { "PDF:PomQuarkB", Param::PDF_PomQuarkB },
  { "PDF:PomQuarkFrac", Param::PDF_PomQuarkFrac },
  { "PDF:PomRescale", Param::PDF_PomRescale },
  { "PDF:PomStrangeSupp", Param::PDF_PomStrangeSupp },
  { "PhaseSpace:bias2SelectionPow", Param::PhaseSpace_bias2SelectionPow },
  { "PhaseSpace:bias2SelectionRef", Param::PhaseSpace_bias2SelectionRef },
  { "PhaseSpace:mHatMax", Param::PhaseSpace_mHatMax },
  { "PhaseSpace:mHatMaxSecond", Param::PhaseSpace_mHatMaxSecond },
  { "PhaseSpace:mHatMin", Param::PhaseSpace_mHatMin },
  { "PhaseSpace:mHatMinSecond", Param::PhaseSpace_mHatMinSecond },
  { "PhaseSpace:minWidthBreitWigners", Param::PhaseSpace_minWidthBreitWigners },
  { "PhaseSpace:pTHat3Max", Param::PhaseSpace_pTHat3Max },
  { "PhaseSpace:pTHat3Min", Param::PhaseSpace_pTHat3Min },
  { "PhaseSpace:pTHat5Max", Param::PhaseSpace_pTHat5Max },
  { "PhaseSpace:pTHat5Min", Param::PhaseSpace_pTHat5Min },
  { "PhaseSpace:pTHatMax", Param::PhaseSpace_pTHatMax },
  { "PhaseSpace:pTHatMaxSecond", Param::PhaseSpace_pTHatMaxSecond },
  { "PhaseSpace:pTHatMin", Param::PhaseSpace_pTHatMin },
  { "PhaseSpace:pTHatMinDiverge", Param::PhaseSpace_pTHatMinDiverge },
  { "PhaseSpace:pTHatMinSecond", Param::PhaseSpace_pTHatMinSecond },
  { "PhaseSpace:RsepMin", Param::PhaseSpace_RsepMin },
  { "Pythia:versionNumber", Param::Pythia_versionNumber },
  { "ResonanceWidths:minThreshold", Param::ResonanceWidths_minThreshold },
  { "ResonanceWidths:minWidth", Param::ResonanceWidths_minWidth },
  { "RHadrons:diquarkSpin1", Param::RHadrons_diquarkSpin1 },
  { "RHadrons:maxWidth", Param::RHadrons_maxWidth },
  { "RHadrons:mCollapse", Param::RHadrons_mCollapse },
  { "RHadrons:mOffsetCloud", Param::RHadrons_mOffsetCloud },
  { "RHadrons:probGluinoball", Param::RHadrons_probGluinoball },
  { "SigmaDiffractive:maxAX", Param::SigmaDiffractive_maxAX },
  { "SigmaDiffractive:maxAXB", Param::SigmaDiffractive_maxAXB },
  { "SigmaDiffractive:maxXB", Param::SigmaDiffractive_maxXB },
  { "SigmaDiffractive:maxXX", Param::SigmaDiffractive_maxXX },
  { "SigmaElastic:bSlope", Param::SigmaElastic_bSlope },
  { "SigmaElastic:lambda", Param::SigmaElastic_lambda },
  { "SigmaElastic:phaseConst", Param::SigmaElastic_phaseConst },
  { "SigmaElastic:rho", Param::SigmaElastic_rho },
  { "SigmaElastic:tAbsMin", Param::SigmaElastic_tAbsMin },
  { "SigmaProcess:alphaSvalue", Param::SigmaProcess_alphaSvalue },
  { "SigmaProcess:factorFixScale", Param::SigmaProcess_factorFixScale },
  { "SigmaProcess:factorMultFac", Param::SigmaProcess_factorMultFac },
  { "SigmaProcess:Kfactor", Param::SigmaProcess_Kfactor },
  { "SigmaProcess:renormFixScale", Param::SigmaProcess_renormFixScale },
  { "SigmaProcess:renormMultFac", Param::SigmaProcess_renormMultFac },
  { "SigmaTotal:sigmaAX", Param::SigmaTotal_sigmaAX },
  { "SigmaTotal:sigmaAXB", Param::SigmaTotal_sigmaAXB },
  { "SigmaTotal:sigmaAXB2TeV", Param::SigmaTotal_sigmaAXB2TeV },
  { "SigmaTotal:sigmaEl", Param::SigmaTotal_sigmaEl },
  { "SigmaTotal:sigmaTot", Param::SigmaTotal_sigmaTot },
  { "SigmaTotal:sigmaXB", Param::SigmaTotal_sigmaXB },
  { "SigmaTotal:sigmaXX", Param::SigmaTotal_sigmaXX },
  { "SLHA:minMassSM", Param::SLHA_minMassSM },
  { "SpaceShower:alphaSuseCMW", Param::SpaceShower_alphaSuseCMW },
  { "SpaceShower:alphaSvalue", Param::SpaceShower_alphaSvalue },
  { "SpaceShower:ecmPow", Param::SpaceShower_ecmPow },
  { "SpaceShower:ecmRef", Param::SpaceShower_ecmRef },
  { "SpaceShower:factorMultFac", Param::SpaceShower_factorMultFac },
  { "SpaceShower:fixedFacScale", Param::SpaceShower_fixedFacScale },
  { "SpaceShower:pT0Ref", Param::SpaceShower_pT0Ref },
  { "SpaceShower:pTdampFudge", Param::SpaceShower_pTdampFudge },
  { "SpaceShower:pTmaxFudge", Param::SpaceShower_pTmaxFudge },
  { "SpaceShower:pTmaxFudgeMPI", Param::SpaceShower_pTmaxFudgeMPI },
  { "SpaceShower:pTmin", Param::SpaceShower_pTmin },
  { "SpaceShower:pTminChgL", Param::SpaceShower_pTminChgL },
  { "SpaceShower:pTminChgQ", Param::SpaceShower_pTminChgQ },
  { "SpaceShower:pTminWeak", Param::SpaceShower_pTminWeak },
  { "SpaceShower:renormMultFac", Param::SpaceShower_renormMultFac },
  { "SpaceShower:strengthIntAsym", Param::SpaceShower_strengthIntAsym },
  { "StandardModel:alphaEM0", Param::StandardModel_alphaEM0 },
  { "StandardModel:alphaEMmZ", Param::StandardModel_alphaEMmZ },
  { "StandardModel:GF", Param::StandardModel_GF },
  { "StandardModel:sin2thetaW", Param::StandardModel_sin2thetaW },
  { "StandardModel:sin2thetaWbar", Param::StandardModel_sin2thetaWbar },
  { "StandardModel:Vcb", Param::StandardModel_Vcb },
  { "StandardModel:Vcd", Param::StandardModel_Vcd },
  { "StandardModel:Vcs", Param::StandardModel_Vcs },
  { "StandardModel:Vtb", Param::StandardModel_Vtb },
  { "StandardModel:Vtd", Param::StandardModel_Vtd },
  { "StandardModel:Vts", Param::StandardModel_Vts },
  { "StandardModel:Vub", Param::StandardModel_Vub },
  { "StandardModel:Vud", Param::StandardModel_Vud },
  { "StandardModel:Vus", Param::StandardModel_Vus },
  { "StringFlav:decupletSup", Param::StringFlav_decupletSup },
  { "StringFlav:etaPrimeSup", Param::StringFlav_etaPrimeSup },
  { "StringFlav:etaSup", Param::StringFlav_etaSup },
  { "StringFlav:heavyLeadingBSup", Param::StringFlav_heavyLeadingBSup },
  { "StringFlav:lightLeadingBSup", Param::StringFlav_lightLeadingBSup },
  { "StringFlav:mesonBL1S0J1", Param::StringFlav_mesonBL1S0J1 },
  { "StringFlav:mesonBL1S1J0", Param::StringFlav_mesonBL1S1J0 },
  { "StringFlav:mesonBL1S1J1", Param::StringFlav_mesonBL1S1J1 },
  { "StringFlav:mesonBL1S1J2", Param::StringFlav_mesonBL1S1J2 },
  { "StringFlav:mesonBvector", Param::StringFlav_mesonBvector },
  { "StringFlav:mesonCL1S0J1", Param::StringFlav_mesonCL1S0J1 },
  { "StringFlav:mesonCL1S1J0", Param::StringFlav_mesonCL1S1J0 },
  { "StringFlav:mesonCL1S1J1", Param::StringFlav_mesonCL1S1J1 },
  { "StringFlav:mesonCL1S1J2", Param::StringFlav_mesonCL1S1J2 },
  { "StringFlav:mesonCvector", Param::StringFlav_mesonCvector },
  { "StringFlav:mesonSL1S0J1", Param::StringFlav_mesonSL1S0J1 },
  { "StringFlav:mesonSL1S1J0", Param::StringFlav_mesonSL1S1J0 },
  { "StringFlav:mesonSL1S1J1", Param::StringFlav_mesonSL1S1J1 },
  { "StringFlav:mesonSL1S1J2", Param::StringFlav_mesonSL1S1J2 },
  { "StringFlav:mesonSvector", Param::StringFlav_mesonSvector },
  { "StringFlav:mesonUDL1S0J1", Param::StringFlav_mesonUDL1S0J1 },
  { "StringFlav:mesonUDL1S1J0", Param::StringFlav_mesonUDL1S1J0 },
  { "StringFlav:mesonUDL1S1J1", Param::StringFlav_mesonUDL1S1J1 },
  { "StringFlav:mesonUDL1S1J2", Param::StringFlav_mesonUDL1S1J2 },
  { "StringFlav:mesonUDvector", Param::StringFlav_mesonUDvector },
  { "StringFlav:popcornRate", Param::StringFlav_popcornRate },
  { "StringFlav:popcornSmeson", Param::StringFlav_popcornSmeson },
  { "StringFlav:popcornSpair", Param::StringFlav_popcornSpair },
  { "StringFlav:probQQ1toQQ0", Param::StringFlav_probQQ1toQQ0 },
  { "StringFlav:probQQtoQ", Param::StringFlav_probQQtoQ },
  { "StringFlav:probSQtoQQ", Param::StringFlav_probSQtoQQ },
  { "StringFlav:probStoUD", Param::StringFlav_probStoUD },
  { "StringFlav:thetaL1S0J1", Param::StringFlav_thetaL1S0J1 },
  { "StringFlav:thetaL1S1J0", Param::StringFlav_thetaL1S1J0 },
  { "StringFlav:thetaL1S1J1", Param::StringFlav_thetaL1S1J1 },
  { "StringFlav:thetaL1S1J2", Param::StringFlav_thetaL1S1J2 },
  { "StringFlav:thetaPS", Param::StringFlav_thetaPS },
  { "StringFlav:thetaV", Param::StringFlav_thetaV },
  { "StringFragmentation:eBothLeftJunction", Param::StringFragmentation_eBothLeftJunction },
  { "StringFragmentation:eMaxLeftJunction", Param::StringFragmentation_eMaxLeftJunction },
  { "StringFragmentation:eMinLeftJunction", Param::StringFragmentation_eMinLeftJunction },
  { "StringFragmentation:eNormJunction", Param::StringFragmentation_eNormJunction },
  { "StringFragmentation:stopMass", Param::StringFragmentation_stopMass },
  { "StringFragmentation:stopNewFlav", Param::StringFragmentation_stopNewFlav },
  { "StringFragmentation:stopSmear", Param::StringFragmentation_stopSmear },
  { "StringPT:enhancedFraction", Param::StringPT_enhancedFraction },
  { "StringPT:enhancedWidth", Param::StringPT_enhancedWidth },
  { "StringPT:sigma", Param::StringPT_sigma },
  { "StringZ:aExtraDiquark", Param::StringZ_aExtraDiquark },
  { "StringZ:aExtraSQuark", Param::StringZ_aExtraSQuark },
  { "StringZ:aLund", Param::StringZ_aLund },
  { "StringZ:aNonstandardB", Param::StringZ_aNonstandardB },
  { "StringZ:aNonstandardC", Param::StringZ_aNonstandardC },
  { "StringZ:aNonstandardH", Param::StringZ_aNonstandardH },
  { "StringZ:bLund", Param::StringZ_bLund },
  { "StringZ:bNonstandardB", Param::StringZ_bNonstandardB },
  { "StringZ:bNonstandardC", Param::StringZ_bNonstandardC },
  { "StringZ:bNonstandardH", Param::StringZ_bNonstandardH },
  { "StringZ:epsilonB", Param::StringZ_epsilonB },
  { "StringZ:epsilonC", Param::StringZ_epsilonC },
  { "StringZ:epsilonH", Param::StringZ_epsilonH },
  { "StringZ:rFactB", Param::StringZ_rFactB },
  { "StringZ:rFactC", Param::StringZ_rFactC },
  { "StringZ:rFactH", Param::StringZ_rFactH },
  { "TauDecays:tauPolarization", Param::TauDecays_tauPolarization },
  { "TimeShower:alphaSvalue", Param::TimeShower_alphaSvalue },
  { "TimeShower:factorMultFac", Param::TimeShower_factorMultFac },
  { "TimeShower:fixedFacScale", Param::TimeShower_fixedFacScale },
  { "TimeShower:mMaxGamma", Param::TimeShower_mMaxGamma },
  { "TimeShower:octetOniumColFac", Param::TimeShower_octetOniumColFac },
  { "TimeShower:octetOniumFraction", Param::TimeShower_octetOniumFraction },
  { "TimeShower:pTdampFudge", Param::TimeShower_pTdampFudge },
  { "TimeShower:pTmaxFudge", Param::TimeShower_pTmaxFudge },
  { "TimeShower:pTmaxFudgeMPI", Param::TimeShower_pTmaxFudgeMPI },
  { "TimeShower:pTmin", Param::TimeShower_pTmin },
  { "TimeShower:pTminChgL", Param::TimeShower_pTminChgL },
  { "TimeShower:pTminChgQ", Param::TimeShower_pTminChgQ },
  { "TimeShower:pTminWeak", Param::TimeShower_pTminWeak },
  { "TimeShower:renormMultFac", Param::TimeShower_renormMultFac },
  { "TimeShower:scaleGluonToQuark", Param::TimeShower_scaleGluonToQuark },
  { "WeakShower:enhancement", Param::WeakShower_enhancement },
  { "WeakShower:vetoWeakDeltaR", Param::WeakShower_vetoWeakDeltaR },
  { "Wprime:al", Param::Wprime_al },
  { "Wprime:anglesWZ", Param::Wprime_anglesWZ },
  { "Wprime:aq", Param::Wprime_aq },
  { "Wprime:coup2WZ", Param::Wprime_coup2WZ },
  { "Wprime:vl", Param::Wprime_vl },
  { "Wprime:vq", Param::Wprime_vq },
  { "Zprime:ab", Param::Zprime_ab },
  { "Zprime:abPrime", Param::Zprime_abPrime },
  { "Zprime:ac", Param::Zprime_ac },
  { "Zprime:ad", Param::Zprime_ad },
  { "Zprime:ae", Param::Zprime_ae },
  { "Zprime:amu", Param::Zprime_amu },
  { "Zprime:anglesWW", Param::Zprime_anglesWW },
  { "Zprime:anue", Param::Zprime_anue },
  { "Zprime:anumu", Param::Zprime_anumu },
  { "Zprime:anutau", Param::Zprime_anutau },
  { "Zprime:anutauPrime", Param::Zprime_anutauPrime },
  { "Zprime:as", Param::Zprime_as },
  { "Zprime:at", Param::Zprime_at },
  { "Zprime:atau", Param::Zprime_atau },
  { "Zprime:atauPrime", Param::Zprime_atauPrime },
  { "Zprime:atPrime", Param::Zprime_atPrime },
  { "Zprime:au", Param::Zprime_au },
  { "Zprime:coup2WW", Param::Zprime_coup2WW },
  { "Zprime:vb", Param::Zprime_vb },
  { "Zprime:vbPrime", Param::Zprime_vbPrime },
  { "Zprime:vc", Param::Zprime_vc },
  { "Zprime:vd", Param::Zprime_vd },
  { "Zprime:ve", Param::Zprime_ve },
  { "Zprime:vmu", Param::Zprime_vmu },
  { "Zprime:vnue", Param::Zprime_vnue },
  { "Zprime:vnumu", Param::Zprime_vnumu },
  { "Zprime:vnutau", Param::Zprime_vnutau },
  { "Zprime:vnutauPrime", Param::Zprime_vnutauPrime },
  { "Zprime:vs", Param::Zprime_vs },
  { "Zprime:vt", Param::Zprime_vt },
  { "Zprime:vtau", Param::Zprime_vtau },
  { "Zprime:vtauPrime", Param::Zprime_vtauPrime },
  { "Zprime:vtPrime", Param::Zprime_vtPrime },
  { "Zprime:vu", Param::Zprime_vu },
};

const map<string, Word> WordMap = 
{
  { "Alpgen:file", Word::Alpgen_file },
  { "Beams:LHEF", Word::Beams_LHEF },
  { "Beams:LHEFheader", Word::Beams_LHEFheader },
  { "Main:spareWord1", Word::Main_spareWord1 },
  { "Main:spareWord2", Word::Main_spareWord2 },
  { "Main:spareWord3", Word::Main_spareWord3 },
  { "Merging:Process", Word::Merging_Process },
  { "PDF:p", Word::PDF_p },
  { "PDF:pHardSet", Word::PDF_pHardSet },
  { "PDF:pHardSetB", Word::PDF_pHardSetB },
  { "PDF:piSet", Word::PDF_piSet },
  { "PDF:piSetB", Word::PDF_piSetB },
  { "PDF:pSet", Word::PDF_pSet },
  { "PDF:pSetB", Word::PDF_pSetB },
  { "POWHEG:dir", Word::POWHEG_dir },
  { "SLHA:file", Word::SLHA_file },
};

const map<string, FlagList> FlagListMap = 
{
  { "Bottomonium:gg2bbbar(3DJ)[3DJ(1)]g", FlagList::Bottomonium_gg2bbbar_3DJ__3DJ_1__g },
  { "Bottomonium:gg2bbbar(3DJ)[3PJ(8)]g", FlagList::Bottomonium_gg2bbbar_3DJ__3PJ_8__g },
  { "Bottomonium:gg2bbbar(3PJ)[3PJ(1)]g", FlagList::Bottomonium_gg2bbbar_3PJ__3PJ_1__g },
  { "Bottomonium:gg2bbbar(3PJ)[3S1(8)]g", FlagList::Bottomonium_gg2bbbar_3PJ__3S1_8__g },
  { "Bottomonium:gg2bbbar(3S1)[1S0(8)]g", FlagList::Bottomonium_gg2bbbar_3S1__1S0_8__g },
  { "Bottomonium:gg2bbbar(3S1)[3PJ(8)]g", FlagList::Bottomonium_gg2bbbar_3S1__3PJ_8__g },
  { "Bottomonium:gg2bbbar(3S1)[3S1(1)]g", FlagList::Bottomonium_gg2bbbar_3S1__3S1_1__g },
  { "Bottomonium:gg2bbbar(3S1)[3S1(8)]g", FlagList::Bottomonium_gg2bbbar_3S1__3S1_8__g },
  { "Bottomonium:qg2bbbar(3DJ)[3PJ(8)]q", FlagList::Bottomonium_qg2bbbar_3DJ__3PJ_8__q },
  { "Bottomonium:qg2bbbar(3PJ)[3PJ(1)]q", FlagList::Bottomonium_qg2bbbar_3PJ__3PJ_1__q },
  { "Bottomonium:qg2bbbar(3PJ)[3S1(8)]q", FlagList::Bottomonium_qg2bbbar_3PJ__3S1_8__q },
  { "Bottomonium:qg2bbbar(3S1)[1S0(8)]q", FlagList::Bottomonium_qg2bbbar_3S1__1S0_8__q },
  { "Bottomonium:qg2bbbar(3S1)[3PJ(8)]q", FlagList::Bottomonium_qg2bbbar_3S1__3PJ_8__q },
  { "Bottomonium:qg2bbbar(3S1)[3S1(8)]q", FlagList::Bottomonium_qg2bbbar_3S1__3S1_8__q },
  { "Bottomonium:qqbar2bbbar(3DJ)[3PJ(8)]g", FlagList::Bottomonium_qqbar2bbbar_3DJ__3PJ_8__g },
  { "Bottomonium:qqbar2bbbar(3PJ)[3PJ(1)]g", FlagList::Bottomonium_qqbar2bbbar_3PJ__3PJ_1__g },
  { "Bottomonium:qqbar2bbbar(3PJ)[3S1(8)]g", FlagList::Bottomonium_qqbar2bbbar_3PJ__3S1_8__g },
  { "Bottomonium:qqbar2bbbar(3S1)[1S0(8)]g", FlagList::Bottomonium_qqbar2bbbar_3S1__1S0_8__g },
  { "Bottomonium:qqbar2bbbar(3S1)[3PJ(8)]g", FlagList::Bottomonium_qqbar2bbbar_3S1__3PJ_8__g },
  { "Bottomonium:qqbar2bbbar(3S1)[3S1(8)]g", FlagList::Bottomonium_qqbar2bbbar_3S1__3S1_8__g },
  { "Charmonium:gg2ccbar(3DJ)[3DJ(1)]g", FlagList::Charmonium_gg2ccbar_3DJ__3DJ_1__g },
  { "Charmonium:gg2ccbar(3DJ)[3PJ(8)]g", FlagList::Charmonium_gg2ccbar_3DJ__3PJ_8__g },
  { "Charmonium:gg2ccbar(3PJ)[3PJ(1)]g", FlagList::Charmonium_gg2ccbar_3PJ__3PJ_1__g },
  { "Charmonium:gg2ccbar(3PJ)[3S1(8)]g", FlagList::Charmonium_gg2ccbar_3PJ__3S1_8__g },
  { "Charmonium:gg2ccbar(3S1)[1S0(8)]g", FlagList::Charmonium_gg2ccbar_3S1__1S0_8__g },
  { "Charmonium:gg2ccbar(3S1)[3PJ(8)]g", FlagList::Charmonium_gg2ccbar_3S1__3PJ_8__g },
  { "Charmonium:gg2ccbar(3S1)[3S1(1)]g", FlagList::Charmonium_gg2ccbar_3S1__3S1_1__g },
  { "Charmonium:gg2ccbar(3S1)[3S1(8)]g", FlagList::Charmonium_gg2ccbar_3S1__3S1_8__g },
  { "Charmonium:qg2ccbar(3DJ)[3PJ(8)]q", FlagList::Charmonium_qg2ccbar_3DJ__3PJ_8__q },
  { "Charmonium:qg2ccbar(3PJ)[3PJ(1)]q", FlagList::Charmonium_qg2ccbar_3PJ__3PJ_1__q },
  { "Charmonium:qg2ccbar(3PJ)[3S1(8)]q", FlagList::Charmonium_qg2ccbar_3PJ__3S1_8__q },
  { "Charmonium:qg2ccbar(3S1)[1S0(8)]q", FlagList::Charmonium_qg2ccbar_3S1__1S0_8__q },
  { "Charmonium:qg2ccbar(3S1)[3PJ(8)]q", FlagList::Charmonium_qg2ccbar_3S1__3PJ_8__q },
  { "Charmonium:qg2ccbar(3S1)[3S1(8)]q", FlagList::Charmonium_qg2ccbar_3S1__3S1_8__q },
  { "Charmonium:qqbar2ccbar(3DJ)[3PJ(8)]g", FlagList::Charmonium_qqbar2ccbar_3DJ__3PJ_8__g },
  { "Charmonium:qqbar2ccbar(3PJ)[3PJ(1)]g", FlagList::Charmonium_qqbar2ccbar_3PJ__3PJ_1__g },
  { "Charmonium:qqbar2ccbar(3PJ)[3S1(8)]g", FlagList::Charmonium_qqbar2ccbar_3PJ__3S1_8__g },
  { "Charmonium:qqbar2ccbar(3S1)[1S0(8)]g", FlagList::Charmonium_qqbar2ccbar_3S1__1S0_8__g },
  { "Charmonium:qqbar2ccbar(3S1)[3PJ(8)]g", FlagList::Charmonium_qqbar2ccbar_3S1__3PJ_8__g },
  { "Charmonium:qqbar2ccbar(3S1)[3S1(8)]g", FlagList::Charmonium_qqbar2ccbar_3S1__3S1_8__g },
};

const map<string, ModeList> ModeListMap = 
{
  { "Bottomonium:states(3DJ)", ModeList::Bottomonium_states_3DJ_ },
  { "Bottomonium:states(3PJ)", ModeList::Bottomonium_states_3PJ_ },
  { "Bottomonium:states(3S1)", ModeList::Bottomonium_states_3S1_ },
  { "Charmonium:states(3DJ)", ModeList::Charmonium_states_3DJ_ },
  { "Charmonium:states(3PJ)", ModeList::Charmonium_states_3PJ_ },
  { "Charmonium:states(3S1)", ModeList::Charmonium_states_3S1_ },
  { "SUSY:idVecA", ModeList::SUSY_idVecA },
  { "SUSY:idVecB", ModeList::SUSY_idVecB },
};

const map<string, ParamList> ParamListMap = 
{
  { "Bottomonium:O(3DJ)[3D1(1)]", ParamList::Bottomonium_O_3DJ__3D1_1__ },
  { "Bottomonium:O(3DJ)[3P0(8)]", ParamList::Bottomonium_O_3DJ__3P0_8__ },
  { "Bottomonium:O(3PJ)[3P0(1)]", ParamList::Bottomonium_O_3PJ__3P0_1__ },
  { "Bottomonium:O(3PJ)[3S1(8)]", ParamList::Bottomonium_O_3PJ__3S1_8__ },
  { "Bottomonium:O(3S1)[1S0(8)]", ParamList::Bottomonium_O_3S1__1S0_8__ },
  { "Bottomonium:O(3S1)[3P0(8)]", ParamList::Bottomonium_O_3S1__3P0_8__ },
  { "Bottomonium:O(3S1)[3S1(1)]", ParamList::Bottomonium_O_3S1__3S1_1__ },
  { "Bottomonium:O(3S1)[3S1(8)]", ParamList::Bottomonium_O_3S1__3S1_8__ },
  { "Charmonium:O(3DJ)[3D1(1)]", ParamList::Charmonium_O_3DJ__3D1_1__ },
  { "Charmonium:O(3DJ)[3P0(8)]", ParamList::Charmonium_O_3DJ__3P0_8__ },
  { "Charmonium:O(3PJ)[3P0(1)]", ParamList::Charmonium_O_3PJ__3P0_1__ },
  { "Charmonium:O(3PJ)[3S1(8)]", ParamList::Charmonium_O_3PJ__3S1_8__ },
  { "Charmonium:O(3S1)[1S0(8)]", ParamList::Charmonium_O_3S1__1S0_8__ },
  { "Charmonium:O(3S1)[3P0(8)]", ParamList::Charmonium_O_3S1__3P0_8__ },
  { "Charmonium:O(3S1)[3S1(1)]", ParamList::Charmonium_O_3S1__3S1_1__ },
  { "Charmonium:O(3S1)[3S1(8)]", ParamList::Charmonium_O_3S1__3S1_8__ },
  { "StringFlav:probQQ1toQQ0join", ParamList::StringFlav_probQQ1toQQ0join },
};

enum class Setting
{
  FLAG,
  MODE,
  PARAM,
  WORD,
  FLAGLIST,
  MODELIST,
  PARAMLIST
};

static map<string,Setting> settings;

// ----- helper functions -----

// neatly format a floating-point number
string format(double x)
{
  stringstream tmp;
  if ( x == 0. )
    tmp << fixed << setprecision(1) << setw(12) << x;
  else if ( abs(x) < 0.001 )
    tmp << scientific << setprecision(4) << setw(12) << x;
  else if ( abs(x) < 0.1 )
    tmp << fixed << setprecision(7) << setw(12) << x;
  else if ( abs(x) < 1000. )
    tmp << fixed << setprecision(5) << setw(12) << x;
  else if ( abs(x) < 1000000. )
    tmp << fixed << setprecision(3) << setw(12) << x;
  else
    tmp << scientific << setprecision(4) << setw(12) << x;
  return tmp.str();
}

// @OVERHEAD
// convert string to bool
bool stringToBool(const string& tag)
{
  string tagLow = toLower(tag);
  return ( tagLow == "true" || tagLow == "1" || tagLow == "on" || tagLow == "yes" || tagLow == "ok" );
}

// @OVERHEAD
// Extract XML value string following XML attribute.
// Format: attribute..."..."
string getAttributeValue(const string& line, const string& attribute)
{
  if (line.find(attribute) == string::npos) return "";
  int iBegAttri = line.find(attribute);
  int iBegQuote = line.find("\"", iBegAttri + 1);
  int iEndQuote = line.find("\"", iBegQuote + 1);
  return line.substr(iBegQuote + 1, iEndQuote - iBegQuote - 1);
}

// @OVERHEAD
// Extract XML bool value following XML attribute.
bool getBoolAttributeValue(const string& line, const string& attribute) 
{
  string valString = getAttributeValue(line, attribute);
  if (valString == "") return false;
  return stringToBool(valString);
}

// @OVERHEAD
// Extract XML int value following XML attribute.
int getIntAttributeValue(const string& line, const string& attribute)
{
  string valString = getAttributeValue(line, attribute);
  if (valString == "") return 0;
  istringstream valStream(valString);
  int intVal;
  valStream >> intVal;
  return intVal;
}

// @OVERHEAD
// Extract XML double value following XML attribute.
double getDoubleAttributeValue(const string& line, const string& attribute) 
{
  string valString = getAttributeValue(line, attribute);
  if (valString == "") return 0.;
  istringstream valStream(valString);
  double doubleVal;
  valStream >> doubleVal;
  return doubleVal;
}

// @OVERHEAD
// Extract XML bool vector value following XML attribute.
// Format: 
vector<bool> getBoolVectorAttributeValue(const string& line, const string& attribute)
{
  string valString = getAttributeValue(line, attribute);
  if (valString == "") return vector<bool>(1, false);
  vector<bool> vectorVal;
  size_t       stringPos(0);
  while (stringPos != string::npos) 
  {
    stringPos = valString.find(",");
    istringstream  valStream(valString.substr(0, stringPos));
    valString = valString.substr(stringPos + 1);
    vectorVal.push_back(stringToBool(valStream.str()));
  }
  return vectorVal;
}

// @OVERHEAD
// Extract XML int vector value following XML attribute.
vector<int> getIntVectorAttributeValue(const string& line, const string& attribute)
{
  string valString = getAttributeValue(line, attribute);
  if (valString == "") return vector<int>(1, 0);
  int         intVal;
  vector<int> vectorVal;
  size_t      stringPos(0);
  while (stringPos != string::npos) 
  {
    stringPos = valString.find(",");
    istringstream  valStream(valString.substr(0, stringPos));
    valString = valString.substr(stringPos + 1);
    valStream >> intVal;
    vectorVal.push_back(intVal);
  }
  return vectorVal;

}

// @OVERHEAD
// Extract XML double vector value following XML attribute.
vector<double> getDoubleVectorAttributeValue(const string& line, const string& attribute)
{
  string valString = getAttributeValue(line, attribute);
  if (valString == "") return vector<double>(1, 0.);
  double         doubleVal;
  vector<double> vectorVal;
  size_t         stringPos(0);
  while (stringPos != string::npos) 
  {
    stringPos = valString.find(",");
    istringstream  valStream(valString.substr(0, stringPos));
    valString = valString.substr(stringPos + 1);
    valStream >> doubleVal;
    vectorVal.push_back(doubleVal);
  }
  return vectorVal;
}

// ----- public functions -----

std::vector<bool>::reference Settings::operator[](const Flag flag)
{
  return flagValue[(int)flag];
}

int& Settings::operator[](const Mode mode)
{
  return modeValue[(int)mode];
}

double& Settings::operator[](const Param param)
{
  return paramValue[(int)param];
}

string& Settings::operator[](const Word word)
{
  return wordValue[(int)word];
}

vector<bool>& Settings::operator[](const FlagList flagList)
{
  return flagListValue[(int)flagList];
}

vector<int>& Settings::operator[](const ModeList modeList)
{
  return modeListValue[(int)modeList];
}

vector<double>& Settings::operator[](const ParamList paramList)
{
  return paramListValue[(int)paramList];
}


bool Settings::get(const Flag flag)
{
  return flagValue[(int)flag];
}

int Settings::get(const Mode mode)
{
  return modeValue[(int)mode];
}

double Settings::get(const Param param)
{
  return paramValue[(int)param];
}

const string& Settings::get(const Word word)
{
  return wordValue[(int)word];
}

const vector<bool>& Settings::get(const FlagList flagList)
{
  return flagListValue[(int)flagList];
}

const vector<int>& Settings::get(const ModeList modeList)
{
  return modeListValue[(int)modeList];
}

const vector<double>& Settings::get(const ParamList paramList)
{
  return paramListValue[(int)paramList];
}


bool Settings::getDefault(const Flag flag)
{
  return flagDefault[(int)flag];
}

int Settings::getDefault(const Mode mode)
{
  return modeDefault[(int)mode];
}

double Settings::getDefault(const Param param)
{
  return paramDefault[(int)param];
}

const string& Settings::getDefault(const Word word)
{
  return wordDefault[(int)word];
}

const vector<bool>& Settings::getDefault(const FlagList flagList)
{
  return flagListDefault[(int)flagList];
}

const vector<int>& Settings::getDefault(const ModeList modeList)
{
  return modeListDefault[(int)modeList];
}

const vector<double>& Settings::getDefault(const ParamList paramList)
{
  return paramListDefault[(int)paramList];
}


void Settings::set(Flag flag, bool value)
{
  flagValue[(int)flag] = value;
}

void Settings::set(Mode mode, int value)
{
  int minn = modeMin[(int)mode];
  int maxx = modeMax[(int)mode];
  if (modeFailOutOfRange[(int)mode] && (value < minn || value > maxx))
  {
    // ...
  }
  modeValue[(int)mode] = value;
  modeValue[(int)mode] = min(modeValue[(int)mode], maxx);
  modeValue[(int)mode] = max(modeValue[(int)mode], minn);
  // Tune:ee and Tune:pp each trigger a whole set of changes.
  if (mode == Mode::Tune_ee) initTuneEE();
  if (mode == Mode::Tune_pp) initTunePP();
}

void Settings::set(Param param, double value)
{
  paramValue[(int)param] = value;
  paramValue[(int)param] = min(paramValue[(int)param], paramMax[(int)param]);
  paramValue[(int)param] = max(paramValue[(int)param], paramMin[(int)param]);
}

void Settings::set(Word word, const string& value)
{
  wordValue[(int)word] = value;
}

void Settings::set(FlagList flagList, const vector<bool>& value)
{
  flagListValue[(int)flagList] = value;
}

void Settings::set(ModeList modeList, const vector<int>& value)
{
  modeListValue[(int)modeList] = value;
  for (size_t i=0; i<value.size(); ++i)
  {
    modeListValue[(int)modeList][i] = min(modeListValue[(int)modeList][i], modeListMax[(int)modeList]);
    modeListValue[(int)modeList][i] = max(modeListValue[(int)modeList][i], modeListMin[(int)modeList]);
  }
}

void Settings::set(ParamList paramList, const vector<double>& value)
{
  paramListValue[(int)paramList] = value;
  for (size_t i=0; i<value.size(); ++i)
  {
    paramListValue[(int)paramList][i] = min(paramListValue[(int)paramList][i], paramListMax[(int)paramList]);
    paramListValue[(int)paramList][i] = max(paramListValue[(int)paramList][i], paramListMin[(int)paramList]);
  }
}

void Settings::restoreDefault(Flag flag)
{
  flagValue[(int)flag] = flagDefault[(int)flag];
}

void Settings::restoreDefault(Mode mode)
{
  modeValue[(int)mode] = modeDefault[(int)mode];
  // Tune:ee and Tune:pp each trigger a whole set of changes.
  if (mode == Mode::Tune_ee) initTuneEE();
  if (mode == Mode::Tune_pp) initTunePP();
}

void Settings::restoreDefault(Param param)
{
  paramValue[(int)param] = paramDefault[(int)param];
}

void Settings::restoreDefault(Word word)
{
  wordValue[(int)word] = wordDefault[(int)word];
}

void Settings::restoreDefault(FlagList flagList)
{
  flagListValue[(int)flagList] = flagListDefault[(int)flagList];
}

void Settings::restoreDefault(ModeList modeList)
{
  modeListValue[(int)modeList] = modeListDefault[(int)modeList];
}

void Settings::restoreDefault(ParamList paramList)
{
  paramListValue[(int)paramList] = paramListDefault[(int)paramList];
}

void Settings::restoreDefault()
{
  for (size_t i=0; i<(int)Flag::END; ++i) restoreDefault((Flag)i);
  for (size_t i=0; i<(int)Mode::END; ++i) restoreDefault((Mode)i);
  for (size_t i=0; i<(int)Param::END; ++i) restoreDefault((Param)i);
  for (size_t i=0; i<(int)Word::END; ++i) restoreDefault((Word)i);
  for (size_t i=0; i<(int)FlagList::END; ++i) restoreDefault((FlagList)i);
  for (size_t i=0; i<(int)ModeList::END; ++i) restoreDefault((ModeList)i);
  for (size_t i=0; i<(int)ParamList::END; ++i) restoreDefault((ParamList)i);
}

// ----- backwards compatibility -----

  const vector<bool>& Settings::lookupFlagList(const string& name)
  {
    return flagListValue[(int)FlagListMap.at(name)];
  }

  const vector<double>& Settings::lookupParamList(const string& name)
  {
    return paramListValue[(int)ParamListMap.at(name)];
}

// ----- private functions -----

void Settings::addFlag(const string& name, bool defaultVal)
{
  flagDefault[int(FlagMap.at(name))] = defaultVal;
}

void Settings::addMode(const string& name, int defaultVal, bool hasMin, bool hasMax, int min, int max, bool boundsError)
{
  modeDefault[int(ModeMap.at(name))] = defaultVal;
  modeMin[int(ModeMap.at(name))] = hasMin ? min : std::numeric_limits<int>::min();
  modeMax[int(ModeMap.at(name))] = hasMax ? max : std::numeric_limits<int>::max();
  modeFailOutOfRange[int(ModeMap.at(name))] = boundsError;
}

void Settings::addParam(const string& name, double defaultVal, bool hasMin, bool hasMax, double min, double max)
{
  paramDefault[int(ParamMap.at(name))] = defaultVal;
  paramMin[int(ParamMap.at(name))] = hasMin ? min : std::numeric_limits<double>::min();
  paramMax[int(ParamMap.at(name))] = hasMax ? max : std::numeric_limits<double>::max();
}

void Settings::addWord(const string& name, const string& defaultVal)
{
  wordDefault[int(WordMap.at(name))] = defaultVal;
}

void Settings::addFlagList(const string& name, const vector<bool>& defaultVal)
{
  flagListDefault[int(FlagListMap.at(name))] = defaultVal;
}

void Settings::addModeList(const string& name, const vector<int>& defaultVal, bool hasMin, bool hasMax, int min, int max)
{
  modeListDefault[int(ModeListMap.at(name))] = defaultVal;
  modeListMin[int(ModeListMap.at(name))] = hasMin ? min : std::numeric_limits<int>::min();
  modeListMax[int(ModeListMap.at(name))] = hasMax ? max : std::numeric_limits<int>::max();
}

void Settings::addParamList(const string& name, const vector<double>& defaultVal, bool hasMin, bool hasMax, double min, double max)
{
  paramListDefault[int(ParamListMap.at(name))] = defaultVal;
  paramListMin[int(ParamListMap.at(name))] = hasMin ? min : std::numeric_limits<double>::min();
  paramListMax[int(ParamListMap.at(name))] = hasMax ? max : std::numeric_limits<double>::max();
}

// ----- more functions -----

bool Settings::init(const string& startFile, bool append, bool reinit, ostream& os)
{
  if (reinit) initialized = false;
  if (initialized) return true;
  int nError = 0;

  // fill in full array of settings
  if (settings.empty())
  {
    for (auto& i : FlagMap) settings[i.first] = Setting::FLAG;
    for (auto& i : ModeMap) settings[i.first] = Setting::MODE;
    for (auto& i : ParamMap) settings[i.first] = Setting::PARAM;
    for (auto& i : WordMap) settings[i.first] = Setting::WORD;
    for (auto& i : FlagListMap) settings[i.first] = Setting::FLAGLIST;
    for (auto& i : ModeListMap) settings[i.first] = Setting::MODELIST;
    for (auto& i : ParamListMap) settings[i.first] = Setting::PARAMLIST;
  }

  // clear setting lists
  flagValue.clear();
  flagDefault.clear();
  modeFailOutOfRange.clear();
  modeValue.clear();
  modeDefault.clear();
  modeMax.clear();
  modeMin.clear();
  paramValue.clear();
  paramDefault.clear();
  paramMax.clear();
  paramMin.clear();
  wordValue.clear();
  wordDefault.clear();
  flagListValue.clear();
  flagListDefault.clear();
  modeListValue.clear();
  modeListDefault.clear();
  modeListMax.clear();
  modeListMin.clear();
  paramListValue.clear();
  paramListDefault.clear();
  paramListMax.clear();
  paramListMin.clear();

  // resize setting lists and set a default value
  flagValue.resize((int)Flag::END, false);
  flagDefault.resize((int)Flag::END, false);
  modeFailOutOfRange.resize((int)Mode::END, false);
  modeValue.resize((int)Mode::END, 0);
  modeDefault.resize((int)Mode::END, 0);
  modeMax.resize((int)Mode::END, std::numeric_limits<int>::max());
  modeMin.resize((int)Mode::END, std::numeric_limits<int>::min());
  paramValue.resize((int)Param::END, 0.0);
  paramDefault.resize((int)Param::END, 0.0);
  paramMax.resize((int)Param::END, std::numeric_limits<double>::max());
  paramMin.resize((int)Param::END, std::numeric_limits<double>::min());
  wordValue.resize((int)Word::END);
  wordDefault.resize((int)Word::END);
  flagListValue.resize((int)FlagList::END);
  flagListDefault.resize((int)FlagList::END);
  modeListValue.resize((int)ModeList::END);
  modeListDefault.resize((int)ModeList::END);
  modeListMax.resize((int)ModeList::END, std::numeric_limits<int>::max());
  modeListMin.resize((int)ModeList::END, std::numeric_limits<int>::min());
  paramListValue.resize((int)ParamList::END);
  paramListDefault.resize((int)ParamList::END);
  paramListMax.resize((int)ParamList::END, std::numeric_limits<double>::max());
  paramListMin.resize((int)ParamList::END, std::numeric_limits<double>::min());

  // List of files to be checked. Start with input file.
  vector<string> files;
  files.push_back(startFile);

  // If nontrivial startfile path, then use that for other files as well.
  string pathName = "";
  if (startFile.rfind("/") != string::npos)
    pathName = startFile.substr(0, startFile.rfind("/") + 1);

  // Loop over files. Open them for read.
  for (int i = 0; i < int(files.size()); ++i)
  {
    const char* cstring = files[i].c_str();
    ifstream is(cstring);

    if (!is.good())
    {
      cerr << "\n PYTHIA Error: settings file " << files[i] << " not found" << endl;
      return false;
    }

    // @OVERHEAD
    // Read in one line at a time.
    string line;
    while (getline(is, line))
    {
      // @OVERHEAD
      // Get first word of a line, to interpret it as tag.
      istringstream getfirst(line);
      string tag;
      getfirst >> tag;

      // Skip ahead if not interesting. Only look for new files in startfile.
      if (tag != "<flag" && tag != "<flagfix" && tag != "<mode"
        && tag != "<modeopen" && tag != "<modepick" && tag != "<modefix"
        && tag != "<parm" && tag != "<parmfix" && tag != "<word"
        && tag != "<wordfix" && tag != "<fvec" && tag != "<fvecfix"
        && tag != "<mvec" && tag != "<mvecfix"
        && tag != "<pvec" && tag != "<pvecfix" && tag != "<aidx") continue;

      // @OVERHEAD
      // Read and append continuation line(s) if line does not contain >.
      while (line.find(">") == string::npos) 
      {
        string addLine;
        getline(is, addLine);
        line += " " + addLine;
      }

      // @OVERHEAD
      // Remove extra blanks before an = sign.
      while (line.find(" =") != string::npos) line.erase( line.find(" ="), 1);

      // process file tag.
      if (tag == "<aidx") 
      {
        string name = getAttributeValue(line, "href");
        if (name == "") 
        {
          cerr <<" PYTHIA Error: failed to find name attribute in line " << line << endl;
          ++nError;
          continue;
        }
        files.push_back(pathName + name + ".xml");
        continue;
      }

      // Find name attribute.
      string name = getAttributeValue( line, "name=");
      if (name == "") 
      {
        cerr << " PYTHIA Error: failed to find name attribute in line " << line << endl;
        ++nError;
        continue;
      }

      // Check that default value attribute present, and whether max and min.
      if (line.find("default=") == string::npos) 
      {
        cerr << " PYTHIA Error: failed to find default value token in line " << line << endl;
        ++nError;
        continue;
      }

      bool hasMin = (line.find("min=") != string::npos);
      bool hasMax = (line.find("max=") != string::npos);

      // Check for occurence of a bool and add to flag map.
      if (tag == "<flag" || tag == "<flagfix") 
      {
        bool value = getBoolAttributeValue( line, "default=");
        addFlag( name, value);
      }
      // Check for occurence of an int and add to mode map.
      else if (tag == "<mode" || tag == "<modeopen" || tag == "<modepick" || tag == "<modefix") 
      {
        int value    = getIntAttributeValue( line, "default=");
        int minVal   = getIntAttributeValue( line, "min=");
        int maxVal   = getIntAttributeValue( line, "max=");
        // Enforce check that only allowed options are accepted.
        bool optOnly = false;
        if (tag == "<modepick" && hasMin && hasMax) optOnly = true;
        if (tag == "<modefix") {
          hasMin  = true;
          hasMax  = true;
          minVal  = value;
          maxVal  = value;
          optOnly = true;
        }
        addMode( name, value, hasMin, hasMax, minVal, maxVal, optOnly);
      } 
      // Check for occurence of a double and add to parm map.
      else if (tag == "<parm" || tag == "<parmfix") 
      {
        double value  = getDoubleAttributeValue( line, "default=");
        double minVal = getDoubleAttributeValue( line, "min=");
        double maxVal = getDoubleAttributeValue( line, "max=");
        addParam( name, value, hasMin, hasMax, minVal, maxVal);

      } 
      // Check for occurence of a string and add to word map.
      else if (tag == "<word" || tag == "<wordfix") 
      {
        string value = getAttributeValue( line, "default=");
        addWord( name, value);

      } 
      // Check for occurence of a bool vector and add to fvec map.
      else if (tag == "<fvec" || tag == "<fvecfix") 
      {
        vector<bool> value = getBoolVectorAttributeValue( line, "default=");
        addFlagList( name, value);

      } 
      // Check for occurence of an int vector and add to mvec map.
      else if (tag == "<mvec" || tag == "<mvecfix") 
      {
        vector<int> value = getIntVectorAttributeValue( line, "default=");
        int minVal = getIntAttributeValue( line, "min=");
        int maxVal = getIntAttributeValue( line, "max=");
        addModeList( name, value, hasMin, hasMax, minVal, maxVal);

      } 
      // Check for occurence of a double vector and add to pvec map.
      else if (tag == "<pvec" || tag == "<pvecfix") 
      {
        vector<double> value = getDoubleVectorAttributeValue( line, "default=");
        double minVal = getDoubleAttributeValue( line, "min=");
        double maxVal = getDoubleAttributeValue( line, "max=");
        addParamList( name, value, hasMin, hasMax, minVal, maxVal);
      }
    }
  }

  // set values to defaults
  flagValue = flagDefault;
  modeValue = modeDefault;
  paramValue = paramDefault;
  wordValue = wordDefault;
  flagListValue = flagListDefault;
  modeListValue = modeListDefault;
  paramListValue = paramListDefault;

  // Set up default e+e- and pp tunes, if positive.
  initTuneEE();
  initTunePP();

  if (nError > 0) return false;
  initialized = true;
  return true;
}

bool Settings::readString(const string& line, bool warn, ostream& os) 
{
  // grab name and value strings
  std::smatch match;
  if (!std::regex_match(line, match, std::regex("\\s*(.*?)\\s*=\\s*(.*)"))) return true;
  string name = match[1].str();
  string value = match[2].str();
  
  // chop off anything at end
  // TODO: might break words...
  int i=0; while(value[i] >= 33 && value[i] <= 'z') ++i;
  value.resize(i);

  // If first character is not a letter, then taken to be a comment line.
  if (!isalpha(name[0])) return true;

  // get the setting type
  auto setting = settings.find(name);

  // check if no such setting
  if (setting == settings.end())
  {
    os << "\n PYTHIA Error: variable not regocnized:\n   " << line << endl;
    return false;
  }

  // missing: mode with illegal value
  // missing: name recognised but value not meaningful
  // missing: value is ?

  try
  {
    // update the appropriate setting
    if (setting->second == Setting::FLAG) set(FlagMap.at(name), stringToBool(value));
    if (setting->second == Setting::MODE) set(ModeMap.at(name), std::stoi(value));
    if (setting->second == Setting::PARAM) set(ParamMap.at(name), std::stod(value));
    if (setting->second == Setting::WORD) set(WordMap.at(name), value);

    if (setting->second == Setting::FLAGLIST)
    {
      set(FlagListMap.at(name), getBoolVectorAttributeValue("value=\"" + value + "\"", "value="));
    }
    if (setting->second == Setting::MODELIST)
    {
      set(ModeListMap.at(name), getIntVectorAttributeValue("value=\"" + value + "\"", "value="));
    }
    if (setting->second == Setting::PARAMLIST)
    {
      set(ParamListMap.at(name), getDoubleVectorAttributeValue("value=\"" + value + "\"", "value="));
    }
  }
  catch(std::exception& e)
  {
    if (warn) os << "\n PYTHIA Error: variable recognized, but its value"
      << " not meaningful:\n   " << line << endl;
    return false;
  }
  
  return true;
}

void Settings::printQuiet(bool quiet)
{
  // Switch off as much output as possible.
  if (quiet) {
    set(Flag::Init_showProcesses,               false );
    set(Flag::Init_showMultipartonInteractions, false );
    set(Flag::Init_showChangedSettings,         false );
    set(Flag::Init_showAllSettings,             false );
    set(Flag::Init_showChangedParticleData,     false );
    set(Flag::Init_showChangedResonanceData,    false );
    set(Flag::Init_showAllParticleData,         false );
    set(Mode::Init_showOneParticleData,             0 );
    set(Mode::Next_numberCount,                     0 );
    set(Mode::Next_numberShowLHA,                   0 );
    set(Mode::Next_numberShowInfo,                  0 );
    set(Mode::Next_numberShowProcess,               0 );
    set(Mode::Next_numberShowEvent,                 0 );

  // Restore ouput settings to default.
  } else {
    restoreDefault(Flag::Init_showProcesses);
    restoreDefault(Flag::Init_showMultipartonInteractions);
    restoreDefault(Flag::Init_showChangedSettings);
    restoreDefault(Flag::Init_showAllSettings);
    restoreDefault(Flag::Init_showChangedParticleData);
    restoreDefault(Flag::Init_showChangedResonanceData);
    restoreDefault(Flag::Init_showAllParticleData);
    restoreDefault(Mode::Init_showOneParticleData);
    restoreDefault(Mode::Next_numberCount);
    restoreDefault(Mode::Next_numberShowLHA);
    restoreDefault(Mode::Next_numberShowInfo);
    restoreDefault(Mode::Next_numberShowProcess);
    restoreDefault(Mode::Next_numberShowEvent);
  }
}

bool Settings::writeFile(const string& toFile, bool writeAll) 
{
  // Open file for writing.
  const char* cstring = toFile.c_str();
  ofstream os(cstring);
  if (!os) 
  {
    // Error in Settings::writeFile: could not open file
    return false;
  }

  // Hand over real work to next method.
  return writeFile( os, writeAll);
}

bool Settings::writeFile(ostream& os, bool writeAll)
{
  // so it looks like the format is
  // <setting> = <valueNow>

  // if we don't use writeAll then just those that changed from the default will be printed

  // Write simple header as comment.
  if (writeAll) os << "! List of all current PYTHIA ";
  else          os << "! List of all modified PYTHIA ";
  os << fixed << setprecision(3) << (*this)[Param::Pythia_versionNumber] << " settings.\n";

  for (auto& i : settings)
  {
    os << i.first << " = ";
    static string state[2] = {"off", "on"};
    if (i.second == Setting::FLAG)
    {
      int index = int(FlagMap.at(i.first));
      if (writeAll || flagValue[index] != flagDefault[index]) os << state[flagValue[index]];
    }
    if (i.second == Setting::MODE)
    {
      int index = int(ModeMap.at(i.first));
      if (writeAll || modeValue[index] != modeDefault[index]) os << modeValue[index];
    }
    if (i.second == Setting::PARAM)
    {
      int index = int(ParamMap.at(i.first));
      if (writeAll || paramValue[index] != paramDefault[index]) os << format(paramValue[index]);
    }
    if (i.second == Setting::FLAGLIST)
    {
      int index = int(FlagListMap.at(i.first));
      auto& val = flagListValue[index];
      if (writeAll || val != flagListDefault[index])
      {
        for (auto i = val.begin(); i != --val.end(); ++i) os << state[*i] << ",";
        os << state[*(--val.end())];
      }
    }
    if (i.second == Setting::MODELIST)
    {
      int index = int(ModeListMap.at(i.first));
      auto& val = modeListValue[index];
      if (writeAll || val != modeListDefault[index])
      {
        for (auto i = val.begin(); i != --val.end(); ++i) os << *i << ",";
        os << *(--val.end());
      }
    }
    if (i.second == Setting::PARAMLIST)
    {
      int index = int(ParamListMap.at(i.first));
      auto& val = paramListValue[index];
      if (writeAll || val != paramListDefault[index])
      {
        for (auto i = val.begin(); i != --val.end(); ++i) os << format(*i) << ",";
        os << format(*(--val.end()));
      }
    }
    os << "\n";
  }

  return true;
}

void Settings::list(bool onlyChanged, string filter, ostream& os) {

  // Table header; output for bool as off/on.
  if (!onlyChanged)
    os << "\n *-------  PYTHIA Flag + Mode + Parm + Word + FVec + MVec + PVec "
       << "Settings (all)  ----------------------------------* \n";
  else
    os << "\n *-------  PYTHIA Flag + Mode + Parm + Word + FVec + MVec + PVec "
       << "Settings (changes only)  -------------------------* \n" ;
  if (filter != "")
    os << "Settings (with requested string) -----------------* \n" ;
  os << " |                                                           "
     << "                                                      | \n"
     << " | Name                                          |           "
     << "           Now |      Default         Min         Max | \n"
     << " |                                               |           "
     << "               |                                      | \n";

  // Convert input string to lowercase for match.
  filter = toLower(filter);

  for (auto& i : settings)
  {
    static string state[2] = {"off", "on"};
    string nameLower = toLower(i.first);
    bool matches = (filter == "" || nameLower.find(filter) != string::npos);
    if (!matches) continue;

    if (i.second == Setting::FLAG)
    {
      int index = int(FlagMap.at(i.first));
      auto val = flagValue[index];
      auto valDefault = flagDefault[index];

      if (matches && (!onlyChanged || val != valDefault))
      {
        os << " | " << setw(45) << left; os << i.first << " | " << setw(24) << right;
        os << state[val] << " | " << setw(12) << state[valDefault];
        os << "                         | \n";
      }
    }

    if (i.second == Setting::MODE)
    {
      int index = int(ModeMap.at(i.first));
      auto val = modeValue[index];
      auto valDefault = modeDefault[index];
      auto valMin = modeMin[index];
      auto valMax = modeMax[index];
      if (!onlyChanged || val != valDefault)
      {
        os << " | " << setw(45) << left; os << i.first << " | " << setw(24) << right;
        os << val << " | " << setw(12) << valDefault;
        if (valMin != std::numeric_limits<int>::min()) os << setw(12) << valMin;
        else os << "            ";
        if (valMax != std::numeric_limits<int>::max()) os << setw(12) << valMax;
        else os << "            ";
        os << " | \n";
      }
    }

    if (i.second == Setting::PARAM)
    {
      int index = int(ParamMap.at(i.first));
      auto val = paramValue[index];
      auto valDefault = paramDefault[index];
      auto valMin = paramMin[index];
      auto valMax = paramMax[index];
      if (!onlyChanged || val != valDefault)
      {
        os << " | " << setw(45) << left; os << i.first << " | " << setw(24) << right;
        os << format(val) << " | " << setw(12) << format(valDefault);
        if (valMin != std::numeric_limits<double>::min()) os << setw(12) << format(valMin);
        else os << "            ";
        if (valMax != std::numeric_limits<double>::max()) os << setw(12) << format(valMax);
        else os << "            ";
        os << " | \n";
      }
    }

    if (i.second == Setting::WORD)
    {
      int index = int(WordMap.at(i.first));
      auto& val = wordValue[index];
      auto& valDefault = wordDefault[index];
      if (!onlyChanged || val != valDefault)
      {
        os << " | " << setw(45) << left; os << i.first << " | " << setw(24) << right;
        int blankLeft = max(0, 60 - max(24, int(val.length()) ) - max(12, int(valDefault.length()) ) );
        string blankPad( blankLeft, ' ');
        os << val << " | " << setw(12) << valDefault << blankPad << " | \n";
      }
    }

    if (i.second == Setting::FLAGLIST)
    {
      int index = int(FlagListMap.at(i.first));
      auto& val = flagListValue[index];
      auto& valDefault = flagListDefault[index];

      if (!onlyChanged || val != valDefault)
      {
        os << " | " << setw(45) << left; os << i.first << " | " << setw(24) << right;
        for (int i = 0; i < (int) val.size() || i < valDefault.size(); ++i) 
        {
          if (i != 0) os << " | " << setw(45) << ' ' << right << " |             ";
          if (i < val.size() ) os << setw(12) << state[val[i]] << " | ";
          else os << "            ";
          if (i < valDefault.size() ) os << setw(12) << state[valDefault[i]];
          else os << "            ";
          os << "            " << "            " << " | \n";
        }
      }
    }

    if (i.second == Setting::MODELIST)
    {
      int index = int(ModeListMap.at(i.first));
      auto& val = modeListValue[index];
      auto& valDefault = modeListDefault[index];
      auto valMin = modeListMin[index];
      auto valMax = modeListMax[index];

      if (!onlyChanged || val != valDefault)
      {
        os << " | " << setw(45) << left; os << i.first << " | " << setw(24) << right;
        for (int i = 0; i < (int) val.size() || i < valDefault.size(); ++i) 
        {
          if (i != 0) os << " | " << setw(45) << ' ' << right << " |             ";
          if (i < val.size() ) os << setw(12) << val[i] << " | ";
          else os << "            ";
          if (i < valDefault.size() ) os << setw(12) << valDefault[i];
          else os << "            ";
          if (valMin != std::numeric_limits<int>::min()) os << setw(12) << valMin;
          else os << "            ";
          if (valMax != std::numeric_limits<int>::max()) os << setw(12) << valMax;
          else os << "            ";
          os << " | \n";
        }
      }
    }

    if (i.second == Setting::PARAMLIST)
    {
      int index = int(ParamListMap.at(i.first));
      auto& val = paramListValue[index];
      auto& valDefault = paramListDefault[index];
      auto valMin = paramListMin[index];
      auto valMax = paramListMax[index];

      if (!onlyChanged || val != valDefault)
      {
        os << " | " << setw(45) << left; os << i.first << " | " << setw(24) << right;
        for (int i = 0; i < (int) val.size() || i < valDefault.size(); ++i) 
        {
          if (i != 0) os << " | " << setw(45) << ' ' << right << " |             ";
          if (i < val.size() ) os << setw(12) << format(val[i]) << " | ";
          else os << "            ";
          if (i < valDefault.size() ) os << setw(12) << format(valDefault[i]);
          else os << "            ";
          if (valMin != std::numeric_limits<double>::min()) os << setw(12) << format(valMin);
          else os << "            ";
          if (valMax != std::numeric_limits<double>::max()) os << setw(12) << format(valMax);
          else os << "            ";
          os << " | \n";
        }
      }
    }
  }

  // End of loop over database contents.
  os << " |                                                           "
     << "                                                      | \n"
     << " *-------  End PYTHIA Flag + Mode + Parm + Word + FVec + MVec + PVec "
     << "Settings  ------------------------------------* " << endl;

}

void Settings::resetTuneEE() 
{
  // Flavour composition.
  restoreDefault(Param::StringFlav_probStoUD);
  restoreDefault(Param::StringFlav_probQQtoQ);
  restoreDefault(Param::StringFlav_probSQtoQQ);
  restoreDefault(Param::StringFlav_probQQ1toQQ0);
  restoreDefault(Param::StringFlav_mesonUDvector);
  restoreDefault(Param::StringFlav_mesonSvector);
  restoreDefault(Param::StringFlav_mesonCvector);
  restoreDefault(Param::StringFlav_mesonBvector);
  restoreDefault(Param::StringFlav_etaSup);
  restoreDefault(Param::StringFlav_etaPrimeSup);
  restoreDefault(Param::StringFlav_popcornSpair);
  restoreDefault(Param::StringFlav_popcornSmeson);
  restoreDefault(Flag::StringFlav_suppressLeadingB);

  // String breaks: z.
  restoreDefault(Param::StringZ_aLund);
  restoreDefault(Param::StringZ_bLund);
  restoreDefault(Param::StringZ_aExtraSQuark);
  restoreDefault(Param::StringZ_aExtraDiquark);
  restoreDefault(Param::StringZ_rFactC);
  restoreDefault(Param::StringZ_rFactB);

  // String breaks: pT.
  restoreDefault(Param::StringPT_sigma);
  restoreDefault(Param::StringPT_enhancedFraction);
  restoreDefault(Param::StringPT_enhancedWidth);

  // FSR: strong coupling, IR cutoff.
  restoreDefault(Param::TimeShower_alphaSvalue);
  restoreDefault(Mode::TimeShower_alphaSorder);
  restoreDefault(Flag::TimeShower_alphaSuseCMW);
  restoreDefault(Param::TimeShower_pTmin);
  restoreDefault(Param::TimeShower_pTminChgQ);
}

void Settings::resetTunePP() 
{
  // PDF set.
  restoreDefault(Word::PDF_pSet);

  // Hard matrix elements alpha_s value.
  restoreDefault(Param::SigmaProcess_alphaSvalue);

  // Diffraction: cross sections and mass distributions.
  restoreDefault(Flag::SigmaTotal_zeroAXB);
  restoreDefault(Flag::SigmaDiffractive_dampen);
  restoreDefault(Param::SigmaDiffractive_maxXB);
  restoreDefault(Param::SigmaDiffractive_maxAX);
  restoreDefault(Param::SigmaDiffractive_maxXX);
  restoreDefault(Param::Diffraction_largeMassSuppress);

  // FSR: dipoles to beam, spin correlations.
  restoreDefault(Flag::TimeShower_dampenBeamRecoil);
  restoreDefault(Flag::TimeShower_phiPolAsym);

  // ISR: strong coupling, IR cutoff, coherence and spin correlations.
  restoreDefault(Param::SpaceShower_alphaSvalue);
  restoreDefault(Mode::SpaceShower_alphaSorder);
  restoreDefault(Param::SpaceShower_alphaSuseCMW);
  restoreDefault(Flag::SpaceShower_samePTasMPI);
  restoreDefault(Param::SpaceShower_pT0Ref);
  restoreDefault(Param::SpaceShower_ecmRef);
  restoreDefault(Param::SpaceShower_ecmPow);
  restoreDefault(Param::SpaceShower_pTmaxFudge);
  restoreDefault(Param::SpaceShower_pTdampFudge);
  restoreDefault(Flag::SpaceShower_rapidityOrder);
  restoreDefault(Flag::SpaceShower_phiPolAsym);
  restoreDefault(Flag::SpaceShower_phiIntAsym);

  // MPI: strong coupling, IR regularization, energy scaling.
  restoreDefault(Param::MultipartonInteractions_alphaSvalue);
  restoreDefault(Param::MultipartonInteractions_pT0Ref);
  restoreDefault(Param::MultipartonInteractions_ecmRef);
  restoreDefault(Param::MultipartonInteractions_ecmPow);
  restoreDefault(Mode::MultipartonInteractions_bProfile);
  restoreDefault(Param::MultipartonInteractions_expPow);
  restoreDefault(Param::MultipartonInteractions_a1);

  // Beam remnant parameters.
  restoreDefault(Param::BeamRemnants_primordialKTsoft);
  restoreDefault(Param::BeamRemnants_primordialKThard);
  restoreDefault(Param::BeamRemnants_halfScaleForKT);
  restoreDefault(Param::BeamRemnants_halfMassForKT);

  // Colour reconnection parameters.
  restoreDefault(Mode::ColourReconnection_mode);
  restoreDefault(Param::ColourReconnection_range);
}

void Settings::initTuneEE() 
{
  int eeTune = (*this)[Mode::Tune_ee];
  if (eeTune <= 0) return;

  // Restore all e+e- settings to their original values.
  // Is first step for setting up a specific tune.
  if (eeTune != 0) resetTuneEE();

  // Old flavour and FSR defaults carried over from very old JETSET tune,
  // only with alphaS roughly tuned for "new" pT-ordered shower.
  if (eeTune == 1) {
    set(Param::StringFlav_probStoUD,        0.30  );
    set(Param::StringFlav_probQQtoQ,        0.10  );
    set(Param::StringFlav_probSQtoQQ,       0.40  );
    set(Param::StringFlav_probQQ1toQQ0,     0.05  );
    set(Param::StringFlav_mesonUDvector,    1.00  );
    set(Param::StringFlav_mesonSvector,     1.50  );
    set(Param::StringFlav_mesonCvector,     2.50  );
    set(Param::StringFlav_mesonBvector,     3.00  );
    set(Param::StringFlav_etaSup,           1.00  );
    set(Param::StringFlav_etaPrimeSup,      0.40  );
    set(Param::StringFlav_popcornSpair,     0.50  );
    set(Param::StringFlav_popcornSmeson,    0.50  );
    set(Flag::StringFlav_suppressLeadingB,  false );
    set(Param::StringZ_aLund,               0.30  );
    set(Param::StringZ_bLund,               0.58  );
    set(Param::StringZ_aExtraSQuark,        0.00  );
    set(Param::StringZ_aExtraDiquark,       0.50  );
    set(Param::StringZ_rFactC,              1.00  );
    set(Param::StringZ_rFactB,              1.00  );
    set(Param::StringPT_sigma,              0.36  );
    set(Param::StringPT_enhancedFraction,   0.01  );
    set(Param::StringPT_enhancedWidth,      2.0   );
    set(Param::TimeShower_alphaSvalue,      0.137 );
    set(Mode::TimeShower_alphaSorder,       1     );
    set(Flag::TimeShower_alphaSuseCMW,      false );
    set(Param::TimeShower_pTmin,            0.5   );
    set(Param::TimeShower_pTminChgQ,        0.5   );
  }

  // Marc Montull's tune to particle composition at LEP1 (August 2007).
  else if (eeTune == 2) {
    set(Param::StringFlav_probStoUD,        0.22  );
    set(Param::StringFlav_probQQtoQ,        0.08  );
    set(Param::StringFlav_probSQtoQQ,       0.75  );
    set(Param::StringFlav_probQQ1toQQ0,     0.025 );
    set(Param::StringFlav_mesonUDvector,    0.5   );
    set(Param::StringFlav_mesonSvector,     0.6   );
    set(Param::StringFlav_mesonCvector,     1.5   );
    set(Param::StringFlav_mesonBvector,     2.5   );
    set(Param::StringFlav_etaSup,           0.60  );
    set(Param::StringFlav_etaPrimeSup,      0.15  );
    set(Param::StringFlav_popcornSpair,     1.0   );
    set(Param::StringFlav_popcornSmeson,    1.0   );
    set(Flag::StringFlav_suppressLeadingB,  false );   // kept fixed
    set(Param::StringZ_aLund,               0.76  );
    set(Param::StringZ_bLund,               0.58  );   // kept fixed
    set(Param::StringZ_aExtraSQuark,        0.00  );   // kept fixed
    set(Param::StringZ_aExtraDiquark,       0.50  );   // kept fixed
    set(Param::StringZ_rFactC,              1.00  );   // kept fixed
    set(Param::StringZ_rFactB,              1.00  );   // kept fixed
    set(Param::StringPT_sigma,              0.36  );   // kept fixed
    set(Param::StringPT_enhancedFraction,   0.01  );   // kept fixed
    set(Param::StringPT_enhancedWidth,      2.0   );   // kept fixed
    set(Param::TimeShower_alphaSvalue,      0.137 );   // kept fixed
    set(Mode::TimeShower_alphaSorder,       1     );   // kept fixed
    set(Flag::TimeShower_alphaSuseCMW,      false );   // kept fixed
    set(Param::TimeShower_pTmin,            0.5   );   // kept fixed
    set(Param::TimeShower_pTminChgQ,        0.5   );   // kept fixed
  }

  // Full e+e- tune of flavours and FSR to LEP1 data within the
  // Rivet + Professor framework, by Hendrik Hoeth (June 2009).
  else if (eeTune == 3) {
    set(Param::StringFlav_probStoUD,        0.19  );
    set(Param::StringFlav_probQQtoQ,        0.09  );
    set(Param::StringFlav_probSQtoQQ,       1.00  );
    set(Param::StringFlav_probQQ1toQQ0,     0.027 );
    set(Param::StringFlav_mesonUDvector,    0.62  );
    set(Param::StringFlav_mesonSvector,     0.725 );
    set(Param::StringFlav_mesonCvector,     1.06  );
    set(Param::StringFlav_mesonBvector,     3.0   );
    set(Param::StringFlav_etaSup,           0.63  );
    set(Param::StringFlav_etaPrimeSup,      0.12  );
    set(Param::StringFlav_popcornSpair,     0.5   );   // kept fixed
    set(Param::StringFlav_popcornSmeson,    0.5   );   // kept fixed
    set(Flag::StringFlav_suppressLeadingB,  false );   // kept fixed
    set(Param::StringZ_aLund,               0.3   );   // kept fixed
    set(Param::StringZ_bLund,               0.8   );
    set(Param::StringZ_aExtraSQuark,        0.00  );   // kept fixed
    set(Param::StringZ_aExtraDiquark,       0.50  );   // kept fixed
    set(Param::StringZ_rFactC,              1.00  );   // kept fixed
    set(Param::StringZ_rFactB,              0.67  );
    set(Param::StringPT_sigma,              0.304 );
    set(Param::StringPT_enhancedFraction,   0.01  );   // kept fixed
    set(Param::StringPT_enhancedWidth,      2.0   );   // kept fixed
    set(Param::TimeShower_alphaSvalue,      0.1383);
    set(Mode::TimeShower_alphaSorder,       1     );   // kept fixed
    set(Flag::TimeShower_alphaSuseCMW,      false );   // kept fixed
    set(Param::TimeShower_pTmin,            0.4   );   // kept fixed (near limit)
    set(Param::TimeShower_pTminChgQ,        0.4   );   // kept same as pTmin
  }

  // Full e+e- tune of flavours and FSR to LEP1 data, by Peter Skands
  // (September 2013). Note use of CMW convention for shower.
  else if (eeTune == 4) {
    set(Param::StringFlav_probStoUD,        0.21  );
    set(Param::StringFlav_probQQtoQ,        0.086 );
    set(Param::StringFlav_probSQtoQQ,       1.00  );
    set(Param::StringFlav_probQQ1toQQ0,     0.031 );
    set(Param::StringFlav_mesonUDvector,    0.45  );
    set(Param::StringFlav_mesonSvector,     0.60  );
    set(Param::StringFlav_mesonCvector,     0.95  );
    set(Param::StringFlav_mesonBvector,     3.0   );   // kept fixed
    set(Param::StringFlav_etaSup,           0.65  );
    set(Param::StringFlav_etaPrimeSup,      0.08  );
    set(Param::StringFlav_popcornSpair,     0.5   );   // kept fixed
    set(Param::StringFlav_popcornSmeson,    0.5   );   // kept fixed
    set(Flag::StringFlav_suppressLeadingB,  false );   // kept fixed
    set(Param::StringZ_aLund,               0.55  );
    set(Param::StringZ_bLund,               1.08  );
    set(Param::StringZ_aExtraSQuark,        0.00  );   // kept fixed
    set(Param::StringZ_aExtraDiquark,       1.00  );
    set(Param::StringZ_rFactC,              1.00  );   // kept fixed
    set(Param::StringZ_rFactB,              0.85  );
    set(Param::StringPT_sigma,              0.305 );
    set(Param::StringPT_enhancedFraction,   0.01  );   // kept fixed
    set(Param::StringPT_enhancedWidth,      2.0   );   // kept fixed
    set(Param::TimeShower_alphaSvalue,      0.127 );
    set(Mode::TimeShower_alphaSorder,       1     );   // kept fixed
    set(Flag::TimeShower_alphaSuseCMW,      true  );
    set(Param::TimeShower_pTmin,            0.4   );
    set(Param::TimeShower_pTminChgQ,        0.4   );   // kept same as pTmin
  }

  // First e+e- tune by Nadine Fischer, using eeTune = 3 for flavour
  // composition (September 2013).
  else if (eeTune == 5) {
    set(Param::StringFlav_probStoUD,        0.19  );   // kept fixed
    set(Param::StringFlav_probQQtoQ,        0.09  );   // kept fixed
    set(Param::StringFlav_probSQtoQQ,       1.00  );   // kept fixed
    set(Param::StringFlav_probQQ1toQQ0,     0.027 );   // kept fixed
    set(Param::StringFlav_mesonUDvector,    0.62  );   // kept fixed
    set(Param::StringFlav_mesonSvector,     0.725 );   // kept fixed
    set(Param::StringFlav_mesonCvector,     1.06  );   // kept fixed
    set(Param::StringFlav_mesonBvector,     3.0   );   // kept fixed
    set(Param::StringFlav_etaSup,           0.63  );   // kept fixed
    set(Param::StringFlav_etaPrimeSup,      0.12  );   // kept fixed
    set(Param::StringFlav_popcornSpair,     0.5   );   // kept fixed
    set(Param::StringFlav_popcornSmeson,    0.5   );   // kept fixed
    set(Flag::StringFlav_suppressLeadingB,  false );   // kept fixed
    set(Param::StringZ_aLund,               0.386 );
    set(Param::StringZ_bLund,               0.977 );
    set(Param::StringZ_aExtraSQuark,        0.00  );   // kept fixed
    set(Param::StringZ_aExtraDiquark,       0.940 );
    set(Param::StringZ_rFactC,              1.00  );   // kept fixed
    set(Param::StringZ_rFactB,              0.67  );   // kept fixed
    set(Param::StringPT_sigma,              0.286 );
    set(Param::StringPT_enhancedFraction,   0.01  );   // kept fixed
    set(Param::StringPT_enhancedWidth,      2.0   );   // kept fixed
    set(Param::TimeShower_alphaSvalue,      0.139 );
    set(Mode::TimeShower_alphaSorder,       1     );   // kept fixed
    set(Flag::TimeShower_alphaSuseCMW,      false );   // kept fixed
    set(Param::TimeShower_pTmin,            0.409 );
    set(Param::TimeShower_pTminChgQ,        0.409 );   // kept same as pTmin
  }

  // Second e+e- tune by Nadine Fischer, using eeTune = 3 for flavour
  // composition (September 2013).
  else if (eeTune == 6) {
    set(Param::StringFlav_probStoUD,        0.19  );   // kept fixed
    set(Param::StringFlav_probQQtoQ,        0.09  );   // kept fixed
    set(Param::StringFlav_probSQtoQQ,       1.00  );   // kept fixed
    set(Param::StringFlav_probQQ1toQQ0,     0.027 );   // kept fixed
    set(Param::StringFlav_mesonUDvector,    0.62  );   // kept fixed
    set(Param::StringFlav_mesonSvector,     0.725 );   // kept fixed
    set(Param::StringFlav_mesonCvector,     1.06  );   // kept fixed
    set(Param::StringFlav_mesonBvector,     3.0   );   // kept fixed
    set(Param::StringFlav_etaSup,           0.63  );   // kept fixed
    set(Param::StringFlav_etaPrimeSup,      0.12  );   // kept fixed
    set(Param::StringFlav_popcornSpair,     0.5   );   // kept fixed
    set(Param::StringFlav_popcornSmeson,    0.5   );   // kept fixed
    set(Flag::StringFlav_suppressLeadingB,  false );   // kept fixed
    set(Param::StringZ_aLund,               0.351 );
    set(Param::StringZ_bLund,               0.942 );
    set(Param::StringZ_aExtraSQuark,        0.00  );   // kept fixed
    set(Param::StringZ_aExtraDiquark,       0.547 );
    set(Param::StringZ_rFactC,              1.00  );   // kept fixed
    set(Param::StringZ_rFactB,              0.67  );   // kept fixed
    set(Param::StringPT_sigma,              0.283 );
    set(Param::StringPT_enhancedFraction,   0.01  );   // kept fixed
    set(Param::StringPT_enhancedWidth,      2.0   );   // kept fixed
    set(Param::TimeShower_alphaSvalue,      0.139);
    set(Mode::TimeShower_alphaSorder,       1     );   // kept fixed
    set(Flag::TimeShower_alphaSuseCMW,      false );   // kept fixed
    set(Param::TimeShower_pTmin,            0.406 );
    set(Param::TimeShower_pTminChgQ,        0.406 );   // kept same as pTmin
  }

  // The Monash 2013 tune by Peter Skands, the e+e- part (January 2014).
  else if (eeTune == 7) {
    set(Param::StringFlav_probStoUD,        0.217 );
    set(Param::StringFlav_probQQtoQ,        0.081 );
    set(Param::StringFlav_probSQtoQQ,       0.915 );
    set(Param::StringFlav_probQQ1toQQ0,     0.0275);
    set(Param::StringFlav_mesonUDvector,    0.50  );
    set(Param::StringFlav_mesonSvector,     0.55  );
    set(Param::StringFlav_mesonCvector,     0.88  );
    set(Param::StringFlav_mesonBvector,     2.20  );
    set(Param::StringFlav_etaSup,           0.60  );
    set(Param::StringFlav_etaPrimeSup,      0.12  );
    set(Param::StringFlav_popcornSpair,     0.90  );
    set(Param::StringFlav_popcornSmeson,    0.50  );
    set(Flag::StringFlav_suppressLeadingB,  false );   // kept fixed
    set(Param::StringZ_aLund,               0.68  );
    set(Param::StringZ_bLund,               0.98  );
    set(Param::StringZ_aExtraSQuark,        0.00  );   // kept fixed
    set(Param::StringZ_aExtraDiquark,       0.97  );
    set(Param::StringZ_rFactC,              1.32  );
    set(Param::StringZ_rFactB,              0.855 );
    set(Param::StringPT_sigma,              0.335 );
    set(Param::StringPT_enhancedFraction,   0.01  );   // kept fixed
    set(Param::StringPT_enhancedWidth,      2.0   );   // kept fixed
    set(Param::TimeShower_alphaSvalue,      0.1365);
    set(Mode::TimeShower_alphaSorder,       1     );   // kept fixed
    set(Flag::TimeShower_alphaSuseCMW,      false );   // kept fixed
    set(Param::TimeShower_pTmin,            0.5   );   // kept fixed
    set(Param::TimeShower_pTminChgQ,        0.5   );   // kept fixed
  }
}

void Settings::initTunePP() 
{

  int ppTune = (*this)[Mode::Tune_pp];
  if (ppTune <= 0) return;

  // Restore all pp/ppbar settings to their original values.
  // Is first step for setting up a specific tune.
  if (ppTune != 0) resetTunePP();

  // Set up e+e- tune that goes with the corresponding pp tune.
  if (ppTune > 0) 
  {
    int eeTune = 3;
    if (ppTune == 14 || ppTune >= 18) eeTune = 7;
    // The mode setting is for documentation, the real action is by initTuneEE.
    set(Mode::Tune_ee, eeTune);
    initTuneEE();
  }

  // Decide whether to use LHAPFD where possible.
  int preferLHAPDF = (*this)[Mode::Tune_preferLHAPDF];

  // Old ISR and MPI defaults from early and primitive comparisons with data.
  if (ppTune == 1) {
    set(Word::PDF_pSet,                            "2"   );
    set(Param::SigmaProcess_alphaSvalue,            0.1265);
    set(Flag::SigmaTotal_zeroAXB,                  true  );
    set(Flag::SigmaDiffractive_dampen,             false );
    set(Param::Diffraction_largeMassSuppress,       2.0   );
    set(Flag::TimeShower_dampenBeamRecoil,         false );
    set(Flag::TimeShower_phiPolAsym,               false );
    set(Param::SpaceShower_alphaSvalue,             0.127 );
    set(Mode::SpaceShower_alphaSorder,             1     );
    set(Flag::SpaceShower_alphaSuseCMW,            false );
    set(Flag::SpaceShower_samePTasMPI,             true  );
    set(Param::SpaceShower_pT0Ref,                  2.2   );
    set(Param::SpaceShower_ecmRef,                  1800.0);
    set(Param::SpaceShower_ecmPow,                  0.16  );
    set(Param::SpaceShower_pTmaxFudge,              1.0   );
    set(Param::SpaceShower_pTdampFudge,             1.0   );
    set(Flag::SpaceShower_rapidityOrder,           false );
    set(Flag::SpaceShower_phiPolAsym,              false );
    set(Flag::SpaceShower_phiIntAsym,              false );
    set(Param::MultipartonInteractions_alphaSvalue, 0.127 );
    set(Param::MultipartonInteractions_pT0Ref,      2.15  );
    set(Param::MultipartonInteractions_ecmRef,      1800. );
    set(Param::MultipartonInteractions_ecmPow,      0.16  );
    set(Mode::MultipartonInteractions_bProfile,    2     );
    set(Param::MultipartonInteractions_expPow,      1.0  );
    set(Param::MultipartonInteractions_a1,          0.15  );
    set(Param::BeamRemnants_primordialKTsoft,       0.4   );
    set(Param::BeamRemnants_primordialKThard,       2.1   );
    set(Param::BeamRemnants_halfScaleForKT,         7.0   );
    set(Param::BeamRemnants_halfMassForKT,          2.0   );
    set(Mode::ColourReconnection_mode,               0    );    
    set(Param::ColourReconnection_range,            2.5   );
  }

  // "Tune 1" simple first tune by Peter Skands to ISR and MPI, July 2009.
  else if (ppTune == 2) {
    set(Word::PDF_pSet,                            "2"   );
    set(Param::SigmaProcess_alphaSvalue,            0.1265);
    set(Flag::SigmaTotal_zeroAXB,                  true  );
    set(Flag::SigmaDiffractive_dampen,             false );
    set(Param::Diffraction_largeMassSuppress,       2.0   );
    set(Flag::TimeShower_dampenBeamRecoil,         false );
    set(Flag::TimeShower_phiPolAsym,               false );
    set(Param::SpaceShower_alphaSvalue,             0.137 );
    set(Mode::SpaceShower_alphaSorder,             1     );
    set(Flag::SpaceShower_alphaSuseCMW,            false );
    set(Flag::SpaceShower_samePTasMPI,             false );
    set(Param::SpaceShower_pT0Ref,                  2.0   );
    set(Param::SpaceShower_ecmRef,                  1800.0);
    set(Param::SpaceShower_ecmPow,                  0.0   );
    set(Param::SpaceShower_pTmaxFudge,              1.0   );
    set(Param::SpaceShower_pTdampFudge,             1.0   );
    set(Flag::SpaceShower_rapidityOrder,           false );
    set(Flag::SpaceShower_phiPolAsym,              false );
    set(Flag::SpaceShower_phiIntAsym,              false );
    set(Param::MultipartonInteractions_alphaSvalue, 0.127 );
    set(Param::MultipartonInteractions_pT0Ref,      2.25  );
    set(Param::MultipartonInteractions_ecmRef,      1800. );
    set(Param::MultipartonInteractions_ecmPow,      0.24  );
    set(Mode::MultipartonInteractions_bProfile,    1     );
    set(Param::MultipartonInteractions_expPow,      1.0  );
    set(Param::MultipartonInteractions_a1,          0.15  );
    set(Param::BeamRemnants_primordialKTsoft,       0.5   );
    set(Param::BeamRemnants_primordialKThard,       2.0   );
    set(Param::BeamRemnants_halfScaleForKT,         1.0   );
    set(Param::BeamRemnants_halfMassForKT,          1.0   );
    set(Mode::ColourReconnection_mode,               0    );    
    set(Param::ColourReconnection_range,            10.0  );
  }

  // Tune 2C, July 2010.
  else if (ppTune == 3) {
    set(Word::PDF_pSet,                            "8"   );
    set(Param::SigmaProcess_alphaSvalue,            0.135 );
    set(Flag::SigmaTotal_zeroAXB,                  true  );
    set(Flag::SigmaDiffractive_dampen,             false );
    set(Param::Diffraction_largeMassSuppress,       2.0   );
    set(Flag::TimeShower_dampenBeamRecoil,         true  );
    set(Flag::TimeShower_phiPolAsym,               true  );
    set(Param::SpaceShower_alphaSvalue,             0.137 );
    set(Mode::SpaceShower_alphaSorder,             1     );
    set(Flag::SpaceShower_alphaSuseCMW,            false );
    set(Flag::SpaceShower_samePTasMPI,             false );
    set(Param::SpaceShower_pT0Ref,                  2.0   );
    set(Param::SpaceShower_ecmRef,                  1800.0);
    set(Param::SpaceShower_ecmPow,                  0.0   );
    set(Param::SpaceShower_pTmaxFudge,              1.0   );
    set(Param::SpaceShower_pTdampFudge,             1.0   );
    set(Flag::SpaceShower_rapidityOrder,           true  );
    set(Flag::SpaceShower_phiPolAsym,              true  );
    set(Flag::SpaceShower_phiIntAsym,              true  );
    set(Param::MultipartonInteractions_alphaSvalue, 0.135 );
    set(Param::MultipartonInteractions_pT0Ref,      2.32  );
    set(Param::MultipartonInteractions_ecmRef,      1800. );
    set(Param::MultipartonInteractions_ecmPow,      0.21  );
    set(Mode::MultipartonInteractions_bProfile,    3     );
    set(Param::MultipartonInteractions_expPow,      1.6   );
    set(Param::MultipartonInteractions_a1,          0.15  );
    set(Param::BeamRemnants_primordialKTsoft,       0.5   );
    set(Param::BeamRemnants_primordialKThard,       2.0   );
    set(Param::BeamRemnants_halfScaleForKT,         1.0   );
    set(Param::BeamRemnants_halfMassForKT,          1.0   );
    set(Mode::ColourReconnection_mode,               0    );    
    set(Param::ColourReconnection_range,            3.0   );
  }

  // Tune 2M, July 2010.
  else if (ppTune == 4) {
    set(Word::PDF_pSet,                            "4"   );
    set(Param::SigmaProcess_alphaSvalue,            0.1265);
    set(Flag::SigmaTotal_zeroAXB,                  true  );
    set(Flag::SigmaDiffractive_dampen,             false );
    set(Param::Diffraction_largeMassSuppress,       2.0   );
    set(Flag::TimeShower_dampenBeamRecoil,         true  );
    set(Flag::TimeShower_phiPolAsym,               true  );
    set(Param::SpaceShower_alphaSvalue,             0.130 );
    set(Mode::SpaceShower_alphaSorder,             1     );
    set(Flag::SpaceShower_alphaSuseCMW,            false );
    set(Flag::SpaceShower_samePTasMPI,             false );
    set(Param::SpaceShower_pT0Ref,                  2.0   );
    set(Param::SpaceShower_ecmRef,                  1800.0);
    set(Param::SpaceShower_ecmPow,                  0.0   );
    set(Param::SpaceShower_pTmaxFudge,              1.0   );
    set(Param::SpaceShower_pTdampFudge,             1.0   );
    set(Flag::SpaceShower_rapidityOrder,           true  );
    set(Flag::SpaceShower_phiPolAsym,              true  );
    set(Flag::SpaceShower_phiIntAsym,              true  );
    set(Param::MultipartonInteractions_alphaSvalue, 0.127 );
    set(Param::MultipartonInteractions_pT0Ref,      2.455 );
    set(Param::MultipartonInteractions_ecmRef,      1800. );
    set(Param::MultipartonInteractions_ecmPow,      0.26  );
    set(Mode::MultipartonInteractions_bProfile,    3     );
    set(Param::MultipartonInteractions_expPow,      1.15  );
    set(Param::MultipartonInteractions_a1,          0.15  );
    set(Param::BeamRemnants_primordialKTsoft,       0.5   );
    set(Param::BeamRemnants_primordialKThard,       2.0   );
    set(Param::BeamRemnants_halfScaleForKT,         1.0   );
    set(Param::BeamRemnants_halfMassForKT,          1.0   );
    set(Mode::ColourReconnection_mode,               0    );    
    set(Param::ColourReconnection_range,            3.0   );
  }

  // Tune 4C, October 2010.
  else if (ppTune == 5) {
    set(Word::PDF_pSet,                            "8"   );
    set(Param::SigmaProcess_alphaSvalue,            0.135 );
    set(Flag::SigmaTotal_zeroAXB,                  true  );
    set(Flag::SigmaDiffractive_dampen,             true  );
    set(Param::SigmaDiffractive_maxXB,              65.0  );
    set(Param::SigmaDiffractive_maxAX,              65.0  );
    set(Param::SigmaDiffractive_maxXX,              65.0  );
    set(Param::Diffraction_largeMassSuppress,       2.0   );
    set(Flag::TimeShower_dampenBeamRecoil,         true  );
    set(Flag::TimeShower_phiPolAsym,               true  );
    set(Param::SpaceShower_alphaSvalue,             0.137 );
    set(Mode::SpaceShower_alphaSorder,             1     );
    set(Flag::SpaceShower_alphaSuseCMW,            false );
    set(Flag::SpaceShower_samePTasMPI,             false );
    set(Param::SpaceShower_pT0Ref,                  2.0   );
    set(Param::SpaceShower_ecmRef,                  1800.0);
    set(Param::SpaceShower_ecmPow,                  0.0   );
    set(Param::SpaceShower_pTmaxFudge,              1.0   );
    set(Param::SpaceShower_pTdampFudge,             1.0   );
    set(Flag::SpaceShower_rapidityOrder,           true  );
    set(Flag::SpaceShower_phiPolAsym,              true  );
    set(Flag::SpaceShower_phiIntAsym,              true  );
    set(Param::MultipartonInteractions_alphaSvalue, 0.135 );
    set(Param::MultipartonInteractions_pT0Ref,      2.085 );
    set(Param::MultipartonInteractions_ecmRef,      1800. );
    set(Param::MultipartonInteractions_ecmPow,      0.19  );
    set(Mode::MultipartonInteractions_bProfile,    3     );
    set(Param::MultipartonInteractions_expPow,      2.0   );
    set(Param::MultipartonInteractions_a1,          0.15  );
    set(Param::BeamRemnants_primordialKTsoft,       0.5   );
    set(Param::BeamRemnants_primordialKThard,       2.0   );
    set(Param::BeamRemnants_halfScaleForKT,         1.0   );
    set(Param::BeamRemnants_halfMassForKT,          1.0   );
    set(Mode::ColourReconnection_mode,               0    );    
    set(Param::ColourReconnection_range,            1.5   );
  }

  // Tune 4Cx, January 2011.
  else if (ppTune == 6) {
    set(Word::PDF_pSet,                            "8"   );
    set(Param::SigmaProcess_alphaSvalue,            0.135 );
    set(Flag::SigmaTotal_zeroAXB,                  true  );
    set(Flag::SigmaDiffractive_dampen,             true  );
    set(Param::SigmaDiffractive_maxXB,              65.0  );
    set(Param::SigmaDiffractive_maxAX,              65.0  );
    set(Param::SigmaDiffractive_maxXX,              65.0  );
    set(Param::Diffraction_largeMassSuppress,       2.0   );
    set(Flag::TimeShower_dampenBeamRecoil,         true  );
    set(Flag::TimeShower_phiPolAsym,               true  );
    set(Param::SpaceShower_alphaSvalue,             0.137 );
    set(Mode::SpaceShower_alphaSorder,             1     );
    set(Flag::SpaceShower_alphaSuseCMW,            false );
    set(Flag::SpaceShower_samePTasMPI,             false );
    set(Param::SpaceShower_pT0Ref,                  2.0   );
    set(Param::SpaceShower_ecmRef,                  1800.0);
    set(Param::SpaceShower_ecmPow,                  0.0   );
    set(Param::SpaceShower_pTmaxFudge,              1.0   );
    set(Param::SpaceShower_pTdampFudge,             1.0   );
    set(Flag::SpaceShower_rapidityOrder,           true  );
    set(Flag::SpaceShower_phiPolAsym,              true  );
    set(Flag::SpaceShower_phiIntAsym,              true  );
    set(Param::MultipartonInteractions_alphaSvalue, 0.135 );
    set(Param::MultipartonInteractions_pT0Ref,      2.15  );
    set(Param::MultipartonInteractions_ecmRef,      1800. );
    set(Param::MultipartonInteractions_ecmPow,      0.19  );
    set(Mode::MultipartonInteractions_bProfile,    4     );
    set(Param::MultipartonInteractions_expPow,      1.0   );
    set(Param::MultipartonInteractions_a1,          0.15  );
    set(Param::BeamRemnants_primordialKTsoft,       0.5   );
    set(Param::BeamRemnants_primordialKThard,       2.0   );
    set(Param::BeamRemnants_halfScaleForKT,         1.0   );
    set(Param::BeamRemnants_halfMassForKT,          1.0   );
    set(Mode::ColourReconnection_mode,               0    );    
    set(Param::ColourReconnection_range,            1.5   );
  }

  // The Monash 2013 tune by Peter Skands, the pp part (January 2014).
  else if (ppTune == 14) {
    set(Word::PDF_pSet,                            "13"  );   // NNPDF
    set(Param::SigmaProcess_alphaSvalue,            0.130 );   // same as PDF
    set(Flag::SigmaTotal_zeroAXB,                  true  );
    set(Flag::SigmaDiffractive_dampen,             true  );
    set(Param::SigmaDiffractive_maxXB,              65.0  );
    set(Param::SigmaDiffractive_maxAX,              65.0  );
    set(Param::SigmaDiffractive_maxXX,              65.0  );
    set(Param::Diffraction_largeMassSuppress,       4.0   );
    set(Flag::TimeShower_dampenBeamRecoil,         true  );
    set(Flag::TimeShower_phiPolAsym,               true  );
    set(Param::SpaceShower_alphaSvalue,             0.1365);   // same as FSR
    set(Mode::SpaceShower_alphaSorder,             1     );
    set(Flag::SpaceShower_alphaSuseCMW,            false );
    set(Flag::SpaceShower_samePTasMPI,             false );
    set(Param::SpaceShower_pT0Ref,                  2.0   );
    set(Param::SpaceShower_ecmRef,                  7000.0);
    set(Param::SpaceShower_ecmPow,                  0.0   );
    set(Param::SpaceShower_pTmaxFudge,              1.0   );
    set(Param::SpaceShower_pTdampFudge,             1.0   );
    set(Flag::SpaceShower_rapidityOrder,           true  );
    set(Flag::SpaceShower_phiPolAsym,              true  );
    set(Flag::SpaceShower_phiIntAsym,              true  );
    set(Param::MultipartonInteractions_alphaSvalue, 0.130 );   // same as PDF
    set(Param::MultipartonInteractions_pT0Ref,      2.28  );
    set(Param::MultipartonInteractions_ecmRef,      7000. );
    set(Param::MultipartonInteractions_ecmPow,      0.215 );
    set(Mode::MultipartonInteractions_bProfile,    3     );
    set(Param::MultipartonInteractions_expPow,      1.85  );
    set(Param::MultipartonInteractions_a1,          0.15  );
    set(Param::BeamRemnants_primordialKTsoft,       0.9   );
    set(Param::BeamRemnants_primordialKThard,       1.8   );
    set(Param::BeamRemnants_halfScaleForKT,         1.5   );
    set(Param::BeamRemnants_halfMassForKT,          1.0   );
    set(Mode::ColourReconnection_mode,               0    );    
    set(Param::ColourReconnection_range,            1.80  );
  }

  // Several ATLAS and CMS tunes start out from Tune 4C.
  else if (ppTune < 18) {
    set(Param::SigmaProcess_alphaSvalue,            0.135 );
    set(Flag::SigmaTotal_zeroAXB,                  true  );
    set(Flag::SigmaDiffractive_dampen,             true  );
    set(Param::SigmaDiffractive_maxXB,              65.0  );
    set(Param::SigmaDiffractive_maxAX,              65.0  );
    set(Param::SigmaDiffractive_maxXX,              65.0  );
    set(Param::Diffraction_largeMassSuppress,       2.0   );
    set(Flag::TimeShower_dampenBeamRecoil,         true  );
    set(Flag::TimeShower_phiPolAsym,               true  );
    set(Param::SpaceShower_alphaSvalue,             0.137 );
    set(Mode::SpaceShower_alphaSorder,             1     );
    set(Flag::SpaceShower_alphaSuseCMW,            false );
    set(Flag::SpaceShower_samePTasMPI,             false );
    set(Param::SpaceShower_pT0Ref,                  2.0   );
    set(Param::SpaceShower_ecmRef,                  1800.0);
    set(Param::SpaceShower_ecmPow,                  0.0   );
    set(Param::SpaceShower_pTmaxFudge,              1.0   );
    set(Param::SpaceShower_pTdampFudge,             1.0   );
    set(Flag::SpaceShower_rapidityOrder,           true );
    set(Flag::SpaceShower_phiPolAsym,              true  );
    set(Flag::SpaceShower_phiIntAsym,              true  );
    set(Param::MultipartonInteractions_alphaSvalue, 0.135 );
    set(Param::MultipartonInteractions_pT0Ref,      2.085 );
    set(Param::MultipartonInteractions_ecmRef,      1800. );
    set(Param::MultipartonInteractions_ecmPow,      0.19  );
    set(Mode::MultipartonInteractions_bProfile,    3     );
    set(Param::MultipartonInteractions_expPow,      2.0   );
    set(Param::MultipartonInteractions_a1,          0.15  );
    set(Param::BeamRemnants_primordialKTsoft,       0.5   );
    set(Param::BeamRemnants_primordialKThard,       2.0   );
    set(Param::BeamRemnants_halfScaleForKT,         1.0   );
    set(Param::BeamRemnants_halfMassForKT,          1.0   );
    set(Mode::ColourReconnection_mode,               0    );    
    set(Param::ColourReconnection_range,            1.5   );

    // Several ATLAS tunes in the A2 and AU2 series, see
    // ATLAS note ATL-PHYS-PUB-2012-003 (August 2012).
    // ATLAS MB tune A2-CTEQ6L1.
    if (ppTune == 7) {
      if (preferLHAPDF == 1)
        set(Word::PDF_pSet,       "LHAPDF5:cteq6ll.LHpdf");
      else if (preferLHAPDF == 2)
        set(Word::PDF_pSet,             "LHAPDF6:cteq6l1");
      else set(Word::PDF_pSet,                     "8"   );
      set(Flag::SpaceShower_rapidityOrder,         false );
      set(Param::MultipartonInteractions_pT0Ref,    2.18  );
      set(Param::MultipartonInteractions_ecmPow,    0.22  );
      set(Mode::MultipartonInteractions_bProfile,  4     );
      set(Param::MultipartonInteractions_expPow,    1.0   );
      set(Param::MultipartonInteractions_a1,        0.06  );
      set(Param::ColourReconnection_range,          1.55  );
    }

    // ATLAS MB tune A2-MSTW2008LO.
    else if (ppTune == 8) {
      if (preferLHAPDF == 1)
        set(Word::PDF_pSet, "LHAPDF5:MSTW2008lo68cl.LHgrid");
      else if (preferLHAPDF == 2)
        set(Word::PDF_pSet,      "LHAPDF6:MSTW2008lo68cl");
      else set(Word::PDF_pSet,                     "5"   );
      set(Flag::SpaceShower_rapidityOrder,         false );
      set(Param::MultipartonInteractions_pT0Ref,    1.90  );
      set(Param::MultipartonInteractions_ecmPow,    0.30  );
      set(Mode::MultipartonInteractions_bProfile,  4     );
      set(Param::MultipartonInteractions_expPow,    1.0   );
      set(Param::MultipartonInteractions_a1,        0.03  );
      set(Param::ColourReconnection_range,          2.28  );
    }

    // ATLAS UE tune AU2-CTEQ6L1.
    if (ppTune == 9) {
      if (preferLHAPDF == 1)
        set(Word::PDF_pSet,       "LHAPDF5:cteq6ll.LHpdf");
      else if (preferLHAPDF == 2)
        set(Word::PDF_pSet,             "LHAPDF6:cteq6l1");
      else set(Word::PDF_pSet,                     "8"   );
      set(Flag::SpaceShower_rapidityOrder,         false );
      set(Param::MultipartonInteractions_pT0Ref,    2.13  );
      set(Param::MultipartonInteractions_ecmPow,    0.21  );
      set(Mode::MultipartonInteractions_bProfile,  4     );
      set(Param::MultipartonInteractions_expPow,    1.0   );
      set(Param::MultipartonInteractions_a1,        0.00  );
      set(Param::ColourReconnection_range,          2.21  );
    }

    // ATLAS UE tune AU2-MSTW2008LO.
    else if (ppTune == 10) {
      if (preferLHAPDF == 1)
        set(Word::PDF_pSet, "LHAPDF5:MSTW2008lo68cl.LHgrid");
      else if (preferLHAPDF == 2)
        set(Word::PDF_pSet,      "LHAPDF6:MSTW2008lo68cl");
      else set(Word::PDF_pSet,                     "5"   );
      set(Flag::SpaceShower_rapidityOrder,         false );
      set(Param::MultipartonInteractions_pT0Ref,    1.87  );
      set(Param::MultipartonInteractions_ecmPow,    0.28  );
      set(Mode::MultipartonInteractions_bProfile,  4     );
      set(Param::MultipartonInteractions_expPow,    1.0   );
      set(Param::MultipartonInteractions_a1,        0.01  );
      set(Param::ColourReconnection_range,          5.32  );
    }

    // ATLAS UE tune AU2-CT10.
    else if (ppTune == 11) {
      if (preferLHAPDF == 2)
        set(Word::PDF_pSet,                "LHAPDF6:CT10");
      else
        set(Word::PDF_pSet,         "LHAPDF5:CT10.LHgrid");
      set(Flag::SpaceShower_rapidityOrder,         false );
      set(Param::MultipartonInteractions_pT0Ref,    1.70  );
      set(Param::MultipartonInteractions_ecmPow,    0.16  );
      set(Mode::MultipartonInteractions_bProfile,  4     );
      set(Param::MultipartonInteractions_expPow,    1.0   );
      set(Param::MultipartonInteractions_a1,        0.10  );
      set(Param::ColourReconnection_range,          4.67  );
    }

    // ATLAS UE tune AU2-MRST2007LO*.
    else if (ppTune == 12) {
      if (preferLHAPDF == 1)
        set(Word::PDF_pSet, "LHAPDF5:MRST2007lomod.LHgrid");
      else if (preferLHAPDF == 2)
        set(Word::PDF_pSet,       "LHAPDF6:MRST2007lomod");
      else set(Word::PDF_pSet,                     "3"   );
      set(Flag::SpaceShower_rapidityOrder,         false );
      set(Param::MultipartonInteractions_pT0Ref,    2.39  );
      set(Param::MultipartonInteractions_ecmPow,    0.24  );
      set(Mode::MultipartonInteractions_bProfile,  4     );
      set(Param::MultipartonInteractions_expPow,    1.0   );
      set(Param::MultipartonInteractions_a1,        0.01  );
      set(Param::ColourReconnection_range,          1.76  );
    }

    // ATLAS UE tune AU2-MRST2007LO**.
    else if (ppTune == 13) {
      if (preferLHAPDF == 1)
        set(Word::PDF_pSet,     "LHAPDF5:MRSTMCal.LHgrid");
      else if (preferLHAPDF == 2)
        set(Word::PDF_pSet,            "LHAPDF6:MRSTMCal");
      else set(Word::PDF_pSet,                     "4"   );
      set(Flag::SpaceShower_rapidityOrder,         false );
      set(Param::MultipartonInteractions_pT0Ref,    2.57  );
      set(Param::MultipartonInteractions_ecmPow,    0.23  );
      set(Mode::MultipartonInteractions_bProfile,  4     );
      set(Param::MultipartonInteractions_expPow,    1.0   );
      set(Param::MultipartonInteractions_a1,        0.01  );
      set(Param::ColourReconnection_range,          1.47  );
    }

    // The CMS UE tunes CUETP8S1-CTEQ6L1 and CUETP8S1-HERAPDF1.5LO,
    // see the note CMS PAS GEN-14-001 (April 2014).
    // CMS UE tune CUETP8S1-CTEQ6L1.
    else if (ppTune == 15) {
      if (preferLHAPDF == 1)
        set(Word::PDF_pSet,       "LHAPDF5:cteq6ll.LHpdf");
      else if (preferLHAPDF == 2)
        set(Word::PDF_pSet,             "LHAPDF6:cteq6l1");
      else set(Word::PDF_pSet,                     "8"   );
      set(Param::MultipartonInteractions_pT0Ref,    2.1006);
      set(Param::MultipartonInteractions_ecmPow,    0.2106);
      set(Param::MultipartonInteractions_expPow,    1.6089);
      set(Param::MultipartonInteractions_a1,        0.00  );
      set(Param::ColourReconnection_range,          3.3126);
    }

    // CMS UE tune CUETP8S1-HERAPDF1.5LO.
    else if (ppTune == 16) {
      if (preferLHAPDF == 2)
        set(Word::PDF_pSet,     "LHAPDF6:HERAPDF15LO_EIG");
      else
        set(Word::PDF_pSet, "LHAPDF5:HERAPDF1.5LO_EIG.LHgrid");
      set(Param::MultipartonInteractions_pT0Ref,    2.0001);
      set(Param::MultipartonInteractions_ecmPow,    0.2499);
      set(Param::MultipartonInteractions_expPow,    1.6905);
      set(Param::MultipartonInteractions_a1,        0.00  );
      set(Param::ColourReconnection_range,          6.0964);
    }

    // ATLAS tune AZ to the Z0/gamma* pTspectrum, see the note
    // CERN-PH-EP-2014-075 [arXiv:1406.3660 [hep-ex]] (June 2014).
    else if (ppTune == 17) {
      set(Param::SpaceShower_alphaSvalue,           0.1237);
      set(Param::SpaceShower_pT0Ref,                0.59  );
      set(Param::MultipartonInteractions_pT0Ref,    2.18  );
      set(Param::BeamRemnants_primordialKThard,     2.0   );
    }
  }

  // Several ATLAS and CMS tunes start out from Monash 2013 tune.
  else if (ppTune >= 18) {
    set(Word::PDF_pSet,                            "13"  );   // NNPDF
    set(Param::SigmaProcess_alphaSvalue,            0.130 );   // same as PDF
    set(Flag::SigmaTotal_zeroAXB,                  true  );
    set(Flag::SigmaDiffractive_dampen,             true  );
    set(Param::SigmaDiffractive_maxXB,              65.0  );
    set(Param::SigmaDiffractive_maxAX,              65.0  );
    set(Param::SigmaDiffractive_maxXX,              65.0  );
    set(Param::Diffraction_largeMassSuppress,       4.0   );
    set(Flag::TimeShower_dampenBeamRecoil,         true  );
    set(Flag::TimeShower_phiPolAsym,               true  );
    set(Param::SpaceShower_alphaSvalue,             0.1365);   // same as FSR
    set(Mode::SpaceShower_alphaSorder,             1     );
    set(Flag::SpaceShower_alphaSuseCMW,            false );
    set(Flag::SpaceShower_samePTasMPI,             false );
    set(Param::SpaceShower_pT0Ref,                  2.0   );
    set(Param::SpaceShower_ecmRef,                  7000.0);
    set(Param::SpaceShower_ecmPow,                  0.0   );
    set(Param::SpaceShower_pTmaxFudge,              1.0   );
    set(Param::SpaceShower_pTdampFudge,             1.0   );
    set(Flag::SpaceShower_rapidityOrder,           true  );
    set(Flag::SpaceShower_phiPolAsym,              true  );
    set(Flag::SpaceShower_phiIntAsym,              true  );
    set(Param::MultipartonInteractions_alphaSvalue, 0.130 );   // same as PDF
    set(Param::MultipartonInteractions_pT0Ref,      2.28  );
    set(Param::MultipartonInteractions_ecmRef,      7000. );
    set(Param::MultipartonInteractions_ecmPow,      0.215 );
    set(Mode::MultipartonInteractions_bProfile,    3     );
    set(Param::MultipartonInteractions_expPow,      1.85  );
    set(Param::MultipartonInteractions_a1,          0.15  );
    set(Param::BeamRemnants_primordialKTsoft,       0.9   );
    set(Param::BeamRemnants_primordialKThard,       1.8   );
    set(Param::BeamRemnants_halfScaleForKT,         1.5   );
    set(Param::BeamRemnants_halfMassForKT,          1.0   );
    set(Mode::ColourReconnection_mode,               0    );    
    set(Param::ColourReconnection_range,            1.80  );

    // CMS tune MonashStar = CUETP8M1-NNPDF2.3LO.
    // See R.D. Field, presentation at MPI@LHC 2014, Krakow, Poland.
    if (ppTune == 18) {
      set(Param::MultipartonInteractions_pT0Ref,    2.4024);
      set(Param::MultipartonInteractions_ecmPow,    0.2521);
      set(Param::MultipartonInteractions_expPow,    1.60  );
    }

    // The ATLAS A14 tunes, central tune with CTEQL1.
    // See ATL-PHYS-PUB-2014-021 (November 2014).
    // Warning: note that TimeShower:alphaSvalue is set here, although
    // normally it would be in the domain of ee tunes. This makes the
    // order of Tune:ee and Tune:pp commands relevant.
    else if (ppTune == 19) {
      if (preferLHAPDF == 1)
        set(Word::PDF_pSet,       "LHAPDF5:cteq6ll.LHpdf");
      else if (preferLHAPDF == 2)
        set(Word::PDF_pSet,             "LHAPDF6:cteq6l1");
      else set(Word::PDF_pSet,                     "8"   );
      set(Param::SigmaProcess_alphaSvalue,          0.144 );
      set(Param::TimeShower_alphaSvalue,            0.126 );
      set(Param::SpaceShower_alphaSvalue,           0.125 );
      set(Param::SpaceShower_pT0Ref,                1.3   );
      set(Param::SpaceShower_pTmaxFudge,            0.95   );
      set(Param::SpaceShower_pTdampFudge,           1.21  );
      set(Param::MultipartonInteractions_alphaSvalue,0.118);
      set(Param::MultipartonInteractions_pT0Ref,    1.98  );
      set(Param::BeamRemnants_primordialKThard,     1.72  );
      set(Param::ColourReconnection_range,          2.08  );
    }

    // The ATLAS A14 tunes, central tune with MSTW2008LO.
    else if (ppTune == 20) {
      if (preferLHAPDF == 1)
        set(Word::PDF_pSet, "LHAPDF5:MSTW2008lo68cl.LHgrid");
      else if (preferLHAPDF == 2)
        set(Word::PDF_pSet,      "LHAPDF6:MSTW2008lo68cl");
      else set(Word::PDF_pSet,                     "5"   );
      set(Param::SigmaProcess_alphaSvalue,          0.140 );
      set(Param::TimeShower_alphaSvalue,            0.129 );
      set(Param::SpaceShower_alphaSvalue,           0.129 );
      set(Param::SpaceShower_pT0Ref,                1.62  );
      set(Param::SpaceShower_pTmaxFudge,            0.92  );
      set(Param::SpaceShower_pTdampFudge,           1.14  );
      set(Param::MultipartonInteractions_alphaSvalue,0.130);
      set(Param::MultipartonInteractions_pT0Ref,    2.28  );
      set(Param::BeamRemnants_primordialKThard,     1.82  );
      set(Param::ColourReconnection_range,          1.87  );
    }

    // The ATLAS A14 tunes, central tune with NNPDF2.3LO.
    else if (ppTune == 21) {
      set(Word::PDF_pSet,                          "13"  );
      set(Param::SigmaProcess_alphaSvalue,          0.140 );
      set(Param::TimeShower_alphaSvalue,            0.127 );
      set(Param::SpaceShower_alphaSvalue,           0.127 );
      set(Param::SpaceShower_pT0Ref,                1.56  );
      set(Param::SpaceShower_pTmaxFudge,            0.91  );
      set(Param::SpaceShower_pTdampFudge,           1.05  );
      set(Param::MultipartonInteractions_alphaSvalue,0.126);
      set(Param::MultipartonInteractions_pT0Ref,    2.09  );
      set(Param::BeamRemnants_primordialKThard,     1.88  );
      set(Param::ColourReconnection_range,          1.71  );
    }

    // The ATLAS A14 tunes, central tune with HERAPDF1.5LO.
    else if (ppTune == 22) {
      if (preferLHAPDF == 2)
        set(Word::PDF_pSet,     "LHAPDF6:HERAPDF15LO_EIG");
      else
        set(Word::PDF_pSet, "LHAPDF5:HERAPDF1.5LO_EIG.LHgrid");
      set(Param::SigmaProcess_alphaSvalue,          0.141 );
      set(Param::TimeShower_alphaSvalue,            0.130 );
      set(Param::SpaceShower_alphaSvalue,           0.128);
      set(Param::SpaceShower_pT0Ref,                1.61  );
      set(Param::SpaceShower_pTmaxFudge,            0.95  );
      set(Param::SpaceShower_pTdampFudge,           1.10  );
      set(Param::MultipartonInteractions_alphaSvalue,0.123);
      set(Param::MultipartonInteractions_pT0Ref,    2.14  );
      set(Param::BeamRemnants_primordialKThard,     1.83  );
      set(Param::ColourReconnection_range,          1.78  );
    }

    // The ATLAS A14 tunes, variation 1+.
    else if (ppTune == 23) {
      set(Word::PDF_pSet,                          "13"  );
      set(Param::SigmaProcess_alphaSvalue,          0.140 );
      set(Param::TimeShower_alphaSvalue,            0.127 );
      set(Param::SpaceShower_alphaSvalue,           0.127 );
      set(Param::SpaceShower_pT0Ref,                1.56  );
      set(Param::SpaceShower_pTmaxFudge,            0.91  );
      set(Param::SpaceShower_pTdampFudge,           1.05  );
      set(Param::MultipartonInteractions_alphaSvalue,0.131);
      set(Param::MultipartonInteractions_pT0Ref,    2.09  );
      set(Param::BeamRemnants_primordialKThard,     1.88  );
      set(Param::ColourReconnection_range,          1.73  );
    }

    // The ATLAS A14 tunes, variation 1-.
    else if (ppTune == 24) {
      set(Word::PDF_pSet,                          "13"  );
      set(Param::SigmaProcess_alphaSvalue,          0.140 );
      set(Param::TimeShower_alphaSvalue,            0.127 );
      set(Param::SpaceShower_alphaSvalue,           0.127 );
      set(Param::SpaceShower_pT0Ref,                1.56  );
      set(Param::SpaceShower_pTmaxFudge,            0.91  );
      set(Param::SpaceShower_pTdampFudge,           1.05  );
      set(Param::MultipartonInteractions_alphaSvalue,0.121);
      set(Param::MultipartonInteractions_pT0Ref,    2.09  );
      set(Param::BeamRemnants_primordialKThard,     1.88  );
      set(Param::ColourReconnection_range,          1.69  );
    }

    // The ATLAS A14 tunes, variation 2+.
    else if (ppTune == 25) {
      set(Word::PDF_pSet,                          "13"  );
      set(Param::SigmaProcess_alphaSvalue,          0.140 );
      set(Param::TimeShower_alphaSvalue,            0.139 );
      set(Param::SpaceShower_alphaSvalue,           0.127 );
      set(Param::SpaceShower_pT0Ref,                1.60  );
      set(Param::SpaceShower_pTmaxFudge,            0.91  );
      set(Param::SpaceShower_pTdampFudge,           1.04  );
      set(Param::MultipartonInteractions_alphaSvalue,0.126);
      set(Param::MultipartonInteractions_pT0Ref,    2.09  );
      set(Param::BeamRemnants_primordialKThard,     1.88  );
      set(Param::ColourReconnection_range,          1.71  );
    }

    // The ATLAS A14 tunes, variation 2-.
    else if (ppTune == 26) {
      set(Word::PDF_pSet,                          "13"  );
      set(Param::SigmaProcess_alphaSvalue,          0.140 );
      set(Param::TimeShower_alphaSvalue,            0.111 );
      set(Param::SpaceShower_alphaSvalue,           0.127 );
      set(Param::SpaceShower_pT0Ref,                1.50  );
      set(Param::SpaceShower_pTmaxFudge,            0.91  );
      set(Param::SpaceShower_pTdampFudge,           1.08  );
      set(Param::MultipartonInteractions_alphaSvalue,0.126);
      set(Param::MultipartonInteractions_pT0Ref,    2.09  );
      set(Param::BeamRemnants_primordialKThard,     1.88  );
      set(Param::ColourReconnection_range,          1.71  );
    }

    // The ATLAS A14 tunes, variation 3a+.
    else if (ppTune == 27) {
      set(Word::PDF_pSet,                          "13"  );
      set(Param::SigmaProcess_alphaSvalue,          0.140 );
      set(Param::TimeShower_alphaSvalue,            0.136 );
      set(Param::SpaceShower_alphaSvalue,           0.127 );
      set(Param::SpaceShower_pT0Ref,                1.67  );
      set(Param::SpaceShower_pTmaxFudge,            0.98  );
      set(Param::SpaceShower_pTdampFudge,           1.36  );
      set(Param::MultipartonInteractions_alphaSvalue,0.125);
      set(Param::MultipartonInteractions_pT0Ref,    2.09  );
      set(Param::BeamRemnants_primordialKThard,     1.88  );
      set(Param::ColourReconnection_range,          1.71  );
    }

    // The ATLAS A14 tunes, variation 3a-.
    else if (ppTune == 28) {
      set(Word::PDF_pSet,                          "13"  );
      set(Param::SigmaProcess_alphaSvalue,          0.140 );
      set(Param::TimeShower_alphaSvalue,            0.124 );
      set(Param::SpaceShower_alphaSvalue,           0.127 );
      set(Param::SpaceShower_pT0Ref,                1.51  );
      set(Param::SpaceShower_pTmaxFudge,            0.88  );
      set(Param::SpaceShower_pTdampFudge,           0.93  );
      set(Param::MultipartonInteractions_alphaSvalue,0.127);
      set(Param::MultipartonInteractions_pT0Ref,    2.09  );
      set(Param::BeamRemnants_primordialKThard,     1.88  );
      set(Param::ColourReconnection_range,          1.71  );
    }

    // The ATLAS A14 tunes, variation 3b+.
    else if (ppTune == 29) {
      set(Word::PDF_pSet,                          "13"  );
      set(Param::SigmaProcess_alphaSvalue,          0.140 );
      set(Param::TimeShower_alphaSvalue,            0.114 );
      set(Param::SpaceShower_alphaSvalue,           0.129 );
      set(Param::SpaceShower_pT0Ref,                1.56  );
      set(Param::SpaceShower_pTmaxFudge,            1.00  );
      set(Param::SpaceShower_pTdampFudge,           1.04  );
      set(Param::MultipartonInteractions_alphaSvalue,0.126);
      set(Param::MultipartonInteractions_pT0Ref,    2.09  );
      set(Param::BeamRemnants_primordialKThard,     1.88  );
      set(Param::ColourReconnection_range,          1.71  );
    }

    // The ATLAS A14 tunes, variation 3b-.
    else if (ppTune == 30) {
      set(Word::PDF_pSet,                          "13"  );
      set(Param::SigmaProcess_alphaSvalue,          0.140 );
      set(Param::TimeShower_alphaSvalue,            0.138 );
      set(Param::SpaceShower_alphaSvalue,           0.126 );
      set(Param::SpaceShower_pT0Ref,                1.56  );
      set(Param::SpaceShower_pTmaxFudge,            0.83  );
      set(Param::SpaceShower_pTdampFudge,           1.07  );
      set(Param::MultipartonInteractions_alphaSvalue,0.126);
      set(Param::MultipartonInteractions_pT0Ref,    2.09  );
      set(Param::BeamRemnants_primordialKThard,     1.88  );
      set(Param::ColourReconnection_range,          1.71  );
    }

    // The ATLAS A14 tunes, variation 3c+.
    else if (ppTune == 31) {
      set(Word::PDF_pSet,                          "13"  );
      set(Param::SigmaProcess_alphaSvalue,          0.140 );
      set(Param::TimeShower_alphaSvalue,            0.127 );
      set(Param::SpaceShower_alphaSvalue,           0.140 );
      set(Param::SpaceShower_pT0Ref,                1.56  );
      set(Param::SpaceShower_pTmaxFudge,            0.91  );
      set(Param::SpaceShower_pTdampFudge,           1.05  );
      set(Param::MultipartonInteractions_alphaSvalue,0.126);
      set(Param::MultipartonInteractions_pT0Ref,    2.09  );
      set(Param::BeamRemnants_primordialKThard,     1.88  );
      set(Param::ColourReconnection_range,          1.71  );
    }

    // The ATLAS A14 tunes, variation 3c-.
    else if (ppTune == 32) {
      set(Word::PDF_pSet,                          "13"  );
      set(Param::SigmaProcess_alphaSvalue,          0.140 );
      set(Param::TimeShower_alphaSvalue,            0.127 );
      set(Param::SpaceShower_alphaSvalue,           0.115 );
      set(Param::SpaceShower_pT0Ref,                1.56  );
      set(Param::SpaceShower_pTmaxFudge,            0.91  );
      set(Param::SpaceShower_pTdampFudge,           1.05  );
      set(Param::MultipartonInteractions_alphaSvalue,0.126);
      set(Param::MultipartonInteractions_pT0Ref,    2.09  );
      set(Param::BeamRemnants_primordialKThard,     1.88  );
      set(Param::ColourReconnection_range,          1.71  );
    }

  }

}

} // end namespace Pythia8