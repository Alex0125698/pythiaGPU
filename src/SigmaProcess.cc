// SigmaProcess.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// SigmaProcess class, and classes derived from it.

#include "Pythia8/SigmaProcess.h"


namespace Pythia8 {

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Conversion of GeV^{-2} to mb for cross section.
const double SigmaProcess::CONVERT2MB    = 0.389380;

// The sum of outgoing masses must not be too close to the cm energy.
const double SigmaProcess::MASSMARGIN    = 0.1;

// Parameters of momentum rescaling procedure: maximally allowed
// relative energy error and number of iterations.
const double SigmaProcess::COMPRELERR = 1e-10;
const int    SigmaProcess::NCOMPSTEP  = 10;

// sets nFinal and a few other things
SigmaProcess::SigmaProcess(ProcessType type_)
  : type(type_)
{
  if (type == ProcessType::P2to0) nFinal = 2;
  if (type == ProcessType::P2to1) nFinal = 1;
  if (type == ProcessType::P2to2) nFinal = 2;
  if (type == ProcessType::P2to3) nFinal = 3;
  if (type == ProcessType::P2to0) convert2mb = false;
  if (type == ProcessType::LHA) convert2mb = false;
  if (type == ProcessType::LHA) isLHA = true;
  if (type == ProcessType::LHA) allowNegativeSigma = (pState->lhaUp->strategy() < 0);
  if (type == ProcessType::LHA) name = "Les Houches User Process(es)";
  if (type == ProcessType::LHA) code = 9999;
  if (type == ProcessType::LHA)
  {
    // At initialization size unknown, so return 0.
    if (pState->lhaUp->sizePart() <= 0) nFinal = 0;

    // TODO: does this actually work yet?

    // Sum up all particles that has first mother = 1.
    nFinal = 0;
    for (int i = 3; i < pState->lhaUp->sizePart(); ++i)
      if (pState->lhaUp->mother1(i) == 1) ++nFinal;
  }
}

// not used
SigmaProcess::~SigmaProcess()
{
  // TODO: should delete if nothing uses a destructor
}

// @OVERHEAD all these aliases should be removed
// creates aliases of beam info and some settings
void SigmaProcess::init(PythiaState* pState_)
{
  // Store pointers.
  pState = pState_;

  // --- alias of some settings / particleData
  // TODO: these need to be deleted

  // Medium heavy fermion masses set massless or not in ME expressions.
  mcME            = (pState->settings.get(Flag::SigmaProcess_cMassiveME))
                  ? pState->particleData.m0(4)  : 0.;
  mbME            = (pState->settings.get(Flag::SigmaProcess_bMassiveME))
                  ? pState->particleData.m0(5)  : 0.;
  mmuME           = (pState->settings.get(Flag::SigmaProcess_muMassiveME))
                  ? pState->particleData.m0(13) : 0.;
  mtauME          = (pState->settings.get(Flag::SigmaProcess_tauMassiveME))
                  ? pState->particleData.m0(15) : 0.;

  // Renormalization scale choice.
  renormScale1    = pState->settings.get(Mode::SigmaProcess_renormScale1);
  renormScale2    = pState->settings.get(Mode::SigmaProcess_renormScale2);
  renormScale3    = pState->settings.get(Mode::SigmaProcess_renormScale3);
  renormScale3VV  = pState->settings.get(Mode::SigmaProcess_renormScale3VV);
  renormMultFac   = pState->settings.get(Param::SigmaProcess_renormMultFac);
  renormFixScale  = pState->settings.get(Param::SigmaProcess_renormFixScale);

  // Factorization scale choice.
  factorScale1    = pState->settings.get(Mode::SigmaProcess_factorScale1);
  factorScale2    = pState->settings.get(Mode::SigmaProcess_factorScale2);
  factorScale3    = pState->settings.get(Mode::SigmaProcess_factorScale3);
  factorScale3VV  = pState->settings.get(Mode::SigmaProcess_factorScale3VV);
  factorMultFac   = pState->settings.get(Param::SigmaProcess_factorMultFac);
  factorFixScale  = pState->settings.get(Param::SigmaProcess_factorFixScale);

}

// adds all the possible partons (based on fluxtype) to inBeamA, inBeamB, inPair arrays
bool SigmaProcess::initFlux()
{
    // Maximum incoming quark flavour.
    int nQuarkIn        = pState->settings.get(Mode::PDFinProcess_nQuarkIn);

    // --- alias of beam info
    // TODO: these need to be deleted

    // Read out some properties of beams to allow shorthand.
    // TODO: how does this differ from the fluxType below?
    // it seems like these are mostly ignored except for lepton beams which are required to match the fluxType
    int idA             = (pState->beamA != 0) ? pState->beamA->id() : 0;
    int idB             = (pState->beamB != 0) ? pState->beamB->id() : 0;
    double mA              = (pState->beamA != 0) ? pState->beamA->m() : 0.;
    double mB              = (pState->beamB != 0) ? pState->beamB->m() : 0.;
    bool isLeptonA       = (pState->beamA != 0) ? pState->beamA->isLepton() : false;
    bool isLeptonB       = (pState->beamB != 0) ? pState->beamB->isLepton() : false;
    bool hasLeptonBeams  = isLeptonA || isLeptonB;

    // set the incoming particles list based on the flux type

    // addBeam: set up PDF's that need to be evaluated for the two beams.
    // addPair: set up pairs of incoming partons from the two beams.
    auto addBeamA = [&](int idInA) { inBeamA.push_back({idInA,0.}); };
    auto addBeamB = [&](int idInB) { inBeamB.push_back({idInB,0.}); };
    auto addPair = [&](int idInA, int idInB) { inPair.push_back({{idInA,0.},{idInB,0.},0.}); };

    // Reset arrays (in case of several init's in same run).
    inBeamA.clear();
    inBeamB.clear();
    inPair.clear();

    switch (fluxType)
    {
      case FluxType::NONE:
        pState->info.errorMsg("Error in SigmaProcess::initFlux; fluxType was not specified");
        return false;  

      case FluxType::GG:
        addBeamA(21);
        addBeamB(21);
        addPair(21, 21);
        break;

      case FluxType::QG:
        for (int i = -nQuarkIn; i <= nQuarkIn; ++i)
        {
          int idNow = (i == 0) ? 21 : i;
          addBeamA(idNow);
          addBeamB(idNow);
        }
        for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
        {
          if (idNow != 0)
          {
            addPair(idNow, 21);
            addPair(21, idNow);
          }
        }
        break;

      case FluxType::QQ:
        for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
        {
          if (idNow == 0) continue;
          addBeamA(idNow);
          addBeamB(idNow);
        }
        for (int id1Now = -nQuarkIn; id1Now <= nQuarkIn; ++id1Now)
        {
          if (id1Now == 0) continue;
          for (int id2Now = -nQuarkIn; id2Now <= nQuarkIn; ++id2Now)
          {
            if (id2Now != 0) addPair(id1Now, id2Now);
          }
        }
        break;

      case FluxType::QQBAR:
        for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
        {
          if (idNow == 0) continue;
          addBeamA(idNow);
          addBeamB(idNow);
        }
        for (int id1Now = -nQuarkIn; id1Now <= nQuarkIn; ++id1Now)
        {
          if (id1Now == 0) continue;
          for (int id2Now = -nQuarkIn; id2Now <= nQuarkIn; ++id2Now)
          {
            if (id2Now != 0 && id1Now * id2Now < 0) addPair(id1Now, id2Now);
          }
        }
        break;

      case FluxType::QQBARSAME:
        for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
        {
          if (idNow == 0) continue;
          addBeamA(idNow);
          addBeamB(idNow);
        }
        for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
        {
          if (idNow != 0) addPair(idNow, -idNow);
        }
        break;

      case FluxType::FF:
        // If beams are leptons then they are also the colliding partons.
        if ( isLeptonA && isLeptonB ) {
          addBeamA(idA);
          addBeamB(idB);
          addPair(idA, idB);
        // First beam is lepton and second is hadron.
        } else if ( isLeptonA ) {
          addBeamA(idA);
          for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
          if (idNow != 0) {
            addBeamB(idNow);
            addPair(idA, idNow);
          }
        // First beam is hadron and second is lepton.
        } else if ( isLeptonB ) {
          addBeamB(idB);
          for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
          if (idNow != 0) {
            addBeamA(idNow);
            addPair(idNow, idB);
          }
        // Hadron beams gives quarks.
        } else {
          for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
          if (idNow != 0) {
            addBeamA(idNow);
            addBeamB(idNow);
          }
          for (int id1Now = -nQuarkIn; id1Now <= nQuarkIn; ++id1Now)
          if (id1Now != 0)
          for (int id2Now = -nQuarkIn; id2Now <= nQuarkIn; ++id2Now)
          if (id2Now != 0)
            addPair(id1Now, id2Now);
        }
        break;

      case FluxType::FFBAR:
        // If beams are leptons then also colliding partons.
        if (isLeptonA && isLeptonB && idA * idB < 0) {
          addBeamA(idA);
          addBeamB(idB);
          addPair(idA, idB);
        // Hadron beams gives quarks.
        } else {
          for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
          if (idNow != 0) {
            addBeamA(idNow);
            addBeamB(idNow);
          }
          for (int id1Now = -nQuarkIn; id1Now <= nQuarkIn; ++id1Now)
          if (id1Now != 0)
          for (int id2Now = -nQuarkIn; id2Now <= nQuarkIn; ++id2Now)
          if (id2Now != 0 && id1Now * id2Now < 0)
            addPair(id1Now, id2Now);
        }
        break;

      case FluxType::FFBARSAME:
        // If beams are antiparticle pair and leptons then also colliding partons.
        if ( idA + idB == 0 && isLeptonA ) {
          addBeamA(idA);
          addBeamB(idB);
          addPair(idA, idB);
        // Else assume both to be hadrons, for better or worse.
        } else {
          for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
          if (idNow != 0) {
            addBeamA(idNow);
            addBeamB(idNow);
          }
          for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
          if (idNow != 0)
            addPair(idNow, -idNow);
        }
        break;

      case FluxType::FFBARCHG:
        // If beams are leptons then also colliding partons.
        if ( isLeptonA && isLeptonB && abs( pState->particleData.chargeType(idA)
                + pState->particleData.chargeType(idB) ) == 3 ) {
          addBeamA(idA);
          addBeamB(idB);
          addPair(idA, idB);
        // Hadron beams gives quarks.
        } else {
          for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
          if (idNow != 0) {
            addBeamA(idNow);
            addBeamB(idNow);
          }
          for (int id1Now = -nQuarkIn; id1Now <= nQuarkIn; ++id1Now)
          if (id1Now != 0)
          for (int id2Now = -nQuarkIn; id2Now <= nQuarkIn; ++id2Now)
          if (id2Now != 0 && id1Now * id2Now < 0
            && (abs(id1Now) + abs(id2Now))%2 == 1) addPair(id1Now, id2Now);
        }
        break;

      case FluxType::FGAMMA:
        // Fermion from incoming side A.
        if ( isLeptonA ) {
          addBeamA(idA);
          addPair(idA, 22);
        } else {
          for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
          if (idNow != 0) {
            addBeamA(idNow);
            addPair(idNow, 22);
          }
        }
        // Fermion from incoming side B.
        if ( isLeptonB ) {
          addBeamB( idB);
          addPair(22, idB);
        } else {
          for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
          if (idNow != 0) {
            addBeamB(idNow);
            addPair(22, idNow);
          }
        }
        // Photons in the beams.
        addBeamA(22);
        addBeamB(22);
        break;

      case FluxType::GGAMMA:
        addBeamA(21);
        addBeamA(22);
        addBeamB(21);
        addBeamB(22);
        addPair(21, 22);
        addPair(22, 21);
        break;

      case FluxType::GAMMAGAMMA:
        addBeamA(22);
        addBeamB(22);
        addPair(22, 22);
        break;
    };

  return true;
}

// (helper)
// stores some kinematic quantities; momentum fractions, Mandelstam variables, Q2 renormalization/factorization scales
// still not sure what these actually are...
void SigmaProcess::store1Kin(double x1in, double x2in, double sHin)
{
  // TODO: so are x1, x the parton momentum fractions??
  // TODO: what is sH? the proton energy? is it supposed to be sHat?

  // Default value only sensible for these processes.
  swapTU = false;

  // Incoming parton momentum fractions and sHat.
  x1Save = x1in;
  x2Save = x2in;
  sH     = sHin;
  mH     = sqrt(sH);
  sH2    = sH * sH;

  // Different options for renormalization scale, but normally sHat.
  Q2RenSave                        = renormMultFac * sH;
  if (renormScale1 == 2) Q2RenSave = renormFixScale;

  // Different options for factorization scale, but normally sHat.
  Q2FacSave                        = factorMultFac * sH;
  if (factorScale1 == 2) Q2FacSave = factorFixScale;

  // Evaluate alpha_strong and alpha_EM.
  alpS   = pState->couplings->alphaS(Q2RenSave);
  alpEM  = pState->couplings->alphaEM(Q2RenSave);
}

// (helper)
void SigmaProcess::store2Kin(double x1in, double x2in, double sHin, double tHin, double m3in, double m4in, double runBW3in, double runBW4in)
{
  // Default ordering of particles 3 and 4.
  // TODO: so whats TU got to do with those particles?
  swapTU   = false;

  // Incoming parton momentum fractions.
  x1Save   = x1in;
  x2Save   = x2in;

  // Incoming masses and their squares.
  bool masslessKin = (id3Mass == 0) && (id4Mass == 0);
  if (masslessKin) 
  {
    m3     = 0.;
    m4     = 0.;
  } 
  else 
  {
    m3     = m3in;
    m4     = m4in;
  }
  mSave[3] = m3;
  mSave[4] = m4;
  s3       = m3 * m3;
  s4       = m4 * m4;

  // Standard Mandelstam variables and their squares.
  sH       = sHin;
  tH       = tHin;
  uH       = (masslessKin) ? -(sH + tH) : s3 + s4 - (sH + tH);
  mH       = sqrt(sH);
  sH2      = sH * sH;
  tH2      = tH * tH;
  uH2      = uH * uH;

  // The nominal Breit-Wigner factors with running width.
  runBW3   = runBW3in;
  runBW4   = runBW4in;

  // Calculate squared transverse momentum.
  pT2 = (masslessKin) ?  tH * uH / sH : (tH * uH - s3 * s4) / sH;

  // Special case: pick scale as if 2 -> 1 process in disguise.
  if (isSChannel) {

    // Different options for renormalization scale, but normally sHat.
    Q2RenSave                        = renormMultFac * sH;
    if (renormScale1 == 2) Q2RenSave = renormFixScale;

    // Different options for factorization scale, but normally sHat.
    Q2FacSave                        = factorMultFac * sH;
    if (factorScale1 == 2) Q2FacSave = factorFixScale;

  // Normal case with "true" 2 -> 2.
  } else {

    // Different options for renormalization scale.
    if (masslessKin)            Q2RenSave = (renormScale2 < 4) ? pT2 : sH;
    else if (renormScale2 == 1) Q2RenSave = pT2 + min(s3, s4);
    else if (renormScale2 == 2) Q2RenSave = sqrt((pT2 + s3) * (pT2 + s4));
    else if (renormScale2 == 3) Q2RenSave = pT2 + 0.5 * (s3 + s4);
    else                        Q2RenSave = sH;
    Q2RenSave                            *= renormMultFac;
    if      (renormScale2 == 5) Q2RenSave = renormFixScale;

    // Different options for factorization scale.
    if (masslessKin)            Q2FacSave = (factorScale2 < 4) ? pT2 : sH;
    else if (factorScale2 == 1) Q2FacSave = pT2 + min(s3, s4);
    else if (factorScale2 == 2) Q2FacSave = sqrt((pT2 + s3) * (pT2 + s4));
    else if (factorScale2 == 3) Q2FacSave = pT2 + 0.5 * (s3 + s4);
    else                        Q2FacSave = sH;
    Q2FacSave                            *= factorMultFac;
    if      (factorScale2 == 5) Q2FacSave = factorFixScale;
  }

  // Evaluate alpha_strong and alpha_EM.
  alpS  = pState->couplings->alphaS(Q2RenSave);
  alpEM = pState->couplings->alphaEM(Q2RenSave);
}

// (helper)
void SigmaProcess::store2KinMPI(double x1in, double x2in, double sHin, double tHin, double uHin, double alpSin, double alpEMin, bool needMasses, double m3in, double m4in)
{
  // Default ordering of particles 3 and 4.
  swapTU    = false;

  // Incoming x values.
  x1Save    = x1in;
  x2Save    = x2in;

  // Standard Mandelstam variables and their squares.
  sH        = sHin;
  tH        = tHin;
  uH        = uHin;
  mH        = sqrt(sH);
  sH2       = sH * sH;
  tH2       = tH * tH;
  uH2       = uH * uH;

  // Strong and electroweak couplings.
  alpS      = alpSin;
  alpEM     = alpEMin;

  // Assume vanishing masses. (Will be modified in final kinematics.)
  m3        = 0.;
  s3        = 0.;
  m4        = 0.;
  s4        = 0.;
  sHBeta    = sH;

  // Scattering angle.
  cosTheta  = (tH - uH) / sH;
  sinTheta  = 2. * sqrtpos( tH * uH ) / sH;

  // In some cases must use masses and redefine meaning of tHat and uHat.
  if (needMasses) {
    m3      = m3in;
    s3      = m3 * m3;
    m4      = m4in;
    s4      = m4 * m4;
    sHMass  = sH - s3 - s4;
    sHBeta  = sqrtpos(sHMass*sHMass - 4. * s3 * s4);
    tH      = -0.5 * (sHMass - sHBeta * cosTheta);
    uH      = -0.5 * (sHMass + sHBeta * cosTheta);
    tH2     = tH * tH;
    uH2     = uH * uH;
  }

  // pT2 with masses (at this stage) included.
  pT2Mass   = 0.25 * sHBeta * pow2(sinTheta);
}

// (helper)
void SigmaProcess::store3Kin(double x1in, double x2in, double sHin, Vec4ref p3cmIn, Vec4ref p4cmIn, Vec4ref p5cmIn, double m3in, double m4in, double m5in, double runBW3in, double runBW4in, double runBW5in)
{
  // Input and complement kinematics for resolved 2 -> 3 process.

  // Default ordering of particles 3 and 4 - not relevant here.
  swapTU   = false;

  // Incoming parton momentum fractions.
  x1Save   = x1in;
  x2Save   = x2in;

  // Incoming masses and their squares.
  if (id3Mass == 0 && id4Mass == 0 && id5Mass == 0) 
  {
    m3     = 0.;
    m4     = 0.;
    m5     = 0.;
  } 
  else 
  {
    m3     = m3in;
    m4     = m4in;
    m5     = m5in;
  }
  mSave[3] = m3;
  mSave[4] = m4;
  mSave[5] = m5;
  s3       = m3 * m3;
  s4       = m4 * m4;
  s5       = m5 * m5;

  // Standard Mandelstam variables and four-momenta in rest frame.
  sH       = sHin;
  mH       = sqrt(sH);
  sH2      = sH * sH;
  p3cm     = p3cmIn;
  p4cm     = p4cmIn;
  p5cm     = p5cmIn;

  // The nominal Breit-Wigner factors with running width.
  runBW3   = runBW3in;
  runBW4   = runBW4in;
  runBW5   = runBW5in;

  // Special case: pick scale as if 2 -> 1 process in disguise.
  if (isSChannel) 
  {
    // Different options for renormalization scale, but normally sHat.
    Q2RenSave = renormMultFac * sH;
    if (renormScale1 == 2) Q2RenSave = renormFixScale;

    // Different options for factorization scale, but normally sHat.
    Q2FacSave = factorMultFac * sH;
    if (factorScale1 == 2) Q2RenSave = factorFixScale;

  // "Normal" 2 -> 3 processes, i.e. not vector boson fusion.
  } 
  else if ( idTchan1 != 23 && idTchan1 != 24 && idTchan2 != 23 && idTchan2 != 24 ) 
  {
    double mT3S = s3 + p3cm.pT2();
    double mT4S = s4 + p4cm.pT2();
    double mT5S = s5 + p5cm.pT2();

    // Different options for renormalization scale.
    if      (renormScale3 == 1) Q2RenSave = min( mT3S, min(mT4S, mT5S) );
    else if (renormScale3 == 2) Q2RenSave = sqrt( mT3S * mT4S * mT5S
      / max( mT3S, max(mT4S, mT5S) ) );
    else if (renormScale3 == 3) Q2RenSave = pow( mT3S * mT4S * mT5S,
                                            1./3. );
    else if (renormScale3 == 4) Q2RenSave = (mT3S + mT4S + mT5S) / 3.;
    else                        Q2RenSave = sH;
    Q2RenSave                            *= renormMultFac;
    if      (renormScale3 == 6) Q2RenSave = renormFixScale;

    // Different options for factorization scale.
    if      (factorScale3 == 1) Q2FacSave = min( mT3S, min(mT4S, mT5S) );
    else if (factorScale3 == 2) Q2FacSave = sqrt( mT3S * mT4S * mT5S
      / max( mT3S, max(mT4S, mT5S) ) );
    else if (factorScale3 == 3) Q2FacSave = pow( mT3S * mT4S * mT5S,
                                            1./3. );
    else if (factorScale3 == 4) Q2FacSave = (mT3S + mT4S + mT5S) / 3.;
    else                        Q2FacSave = sH;
    Q2FacSave                            *= factorMultFac;
    if      (factorScale3 == 6) Q2FacSave = factorFixScale;

  // Vector boson fusion 2 -> 3 processes; recoils in positions 4 and 5.
  } 
  else 
  {
    double sV4   = pow2( pState->particleData.m0(idTchan1) );
    double sV5   = pow2( pState->particleData.m0(idTchan2) );
    double mT3S  = s3  + p3cm.pT2();
    double mTV4S = sV4 + p4cm.pT2();
    double mTV5S = sV5 + p5cm.pT2();

    // Different options for renormalization scale.
    if      (renormScale3VV == 1) Q2RenSave = max( sV4, sV5);
    else if (renormScale3VV == 2) Q2RenSave = sqrt( mTV4S * mTV5S );
    else if (renormScale3VV == 3) Q2RenSave = pow( mT3S * mTV4S * mTV5S,
                                              1./3. );
    else if (renormScale3VV == 4) Q2RenSave = (mT3S * mTV4S * mTV5S) / 3.;
    else                          Q2RenSave = sH;
    Q2RenSave                              *= renormMultFac;
    if      (renormScale3VV == 6) Q2RenSave = renormFixScale;

    // Different options for factorization scale.
    if      (factorScale3VV == 1) Q2FacSave = max( sV4, sV5);
    else if (factorScale3VV == 2) Q2FacSave = sqrt( mTV4S * mTV5S );
    else if (factorScale3VV == 3) Q2FacSave = pow( mT3S * mTV4S * mTV5S,
                                              1./3. );
    else if (factorScale3VV == 4) Q2FacSave = (mT3S * mTV4S * mTV5S) / 3.;
    else                          Q2FacSave = sH;
    Q2FacSave                              *= factorMultFac;
    if      (factorScale3VV == 6) Q2FacSave = factorFixScale;
  }

  // Evaluate alpha_strong and alpha_EM.
  alpS  = pState->couplings->alphaS(Q2RenSave);
  alpEM = pState->couplings->alphaEM(Q2RenSave);

}

// stores the kinematic quantities; calls sigmaKin
void SigmaProcess::set1Kin(double x1in, double x2in, double sHin)
{
  store1Kin(x1in, x2in, sHin);
  sigmaKin();
}

// stores the kinematic quantities; calls sigmaKin
void SigmaProcess::set2Kin(double x1in, double x2in, double sHin, double tHin, double m3in, double m4in, double runBW3in, double runBW4in)
{
  store2Kin( x1in, x2in, sHin, tHin, m3in, m4in, runBW3in, runBW4in);
  sigmaKin();
}

// stores the kinematic quantities; calls sigmaKin
void SigmaProcess::set2KinMPI(double x1in, double x2in, double sHin, double tHin,double uHin, double alpSin, double alpEMin, bool needMasses, double m3in, double m4in)
{
  store2KinMPI( x1in, x2in, sHin, tHin, uHin, alpSin, alpEMin, needMasses, m3in, m4in);
  sigmaKin();
}

// stores the kinematic quantities; calls sigmaKin
void SigmaProcess::set3Kin(double x1in, double x2in, double sHin, Vec4ref p3prel, Vec4ref p4prel, Vec4ref p5prel, double m3in, double m4in, double m5in, double runBW3in, double runBW4in, double runBW5in)
{
  store3Kin( x1in, x2in, sHin, p3prel, p4prel, p5prel, m3in, m4in, m5in, runBW3in, runBW4in, runBW5in);
  sigmaKin();
}

// sets the incoming parton indices; calculates sigmaHat; converts units to mb
double SigmaProcess::sigmaHatWrap(int id1in, int id2in)
{
  id1 = id1in; 
  id2 = id2in;
  double result = sigmaHat();

  if (type == ProcessType::P2to1)
  {
    if (convertM2) 
    {
      result /= 2. * sH; // TODO: what is sH?
      // Convert 2 * pi * delta(p^2 - m^2) to Breit-Wigner with same area.
      // TODO: why do we assume resonance just based on units?
      double mTmp   = pState->particleData.m0(resonanceA);
      double GamTmp = pState->particleData.mWidth(resonanceA);
      result *= 2. * mTmp * GamTmp / ( pow2(sH - mTmp * mTmp) + pow2(mTmp * GamTmp) );
    }
  }

  if (type == ProcessType::P2to2)
  {
    if (convertM2) result /= 16. * M_PI * sH2;
  }

  // SigmaProcess
  if (convert2mb) result *= CONVERT2MB;

  return result;
}

// loop over beam particles; get parton densities; calc hadronic cross-section (summed for all beam candidates)
// (only used in PhaseSpace.cc)
double SigmaProcess::sigmaPDF()
{
  // Since no PDF's there is no difference
  if (type == ProcessType::P2to0) return sigmaHat();

  // K factor, multiplying resolved processes. (But not here for MPI.)
  double Kfactor         = pState->settings.get(Param::SigmaProcess_Kfactor);

  // SigmaProcess

  // Evaluate and store the required parton densities.
  // TODO: do these actually change between different calls?
  //       maybe we should evaluate them in sigmaHatWrap
  // TODO: what are x1Save, Q2FacSave?

  // x1Save is the momentum fraction
  for (int j = 0; j < inBeamA.size(); ++j)
    inBeamA[j].pdf = pState->beamA->xfHard( inBeamA[j].id, x1Save, Q2FacSave);
  for (int j = 0; j < inBeamB.size(); ++j)
    inBeamB[j].pdf = pState->beamB->xfHard( inBeamB[j].id, x2Save, Q2FacSave);

  // store total ??Hadronic?? cross-section here; summed for all parton channels
  sigmaSumSave = 0.;

  // Loop over the allowed pairs incoming parton for this process
  for (int i = 0; i < inPair.size(); ++i) 
  {
    // Evaluate hard-scattering cross section. Include K factor.
    // WARNING: Kfactor is an alias
    // Note: this is where sigmaHatWrap inputs setup; they are just the parton flavours!
    inPair[i].pdfSigma = Kfactor * sigmaHatWrap(inPair[i].A.id, inPair[i].B.id);

    // Multiply by respective parton densities.
    // @OVERHEAD: we are just sifting out the appropriate multiplier
    for (int j = 0; j < inBeamA.size(); ++j)
    {
      if (inPair[i].A.id == inBeamA[j].id) 
      {
        inPair[i].A.pdf      = inBeamA[j].pdf;
        inPair[i].pdfSigma *= inBeamA[j].pdf;
        break;
      }
    }
    for (int j = 0; j < inBeamB.size(); ++j)
    {
      if (inPair[i].B.id == inBeamB[j].id) 
      {
        inPair[i].B.pdf      = inBeamB[j].pdf;
        inPair[i].pdfSigma *= inBeamB[j].pdf;
        break;
      }
    }

    // Sum for all channels.
    sigmaSumSave += inPair[i].pdfSigma;
  }

    // Done.
  return sigmaSumSave;
}

// picks the parton pair at random; with the probability weighted by the pdfSigma
void SigmaProcess::pickInState(int id1in, int id2in)
{
  // Note: so this function just selects one of the parton channels

  // Multiparton interactions: partons already selected.
  // Note: so the inputs are used to select the parton channel in this case
  if (id1in != 0 && id2in != 0) 
  {
    id1 = id1in;
    id2 = id2in;
    return;
  }

  // Pick channel. Extract channel flavours and pdf's.
  // Note: so we pick a channel at random, with the probability weighted by the pdfSigma
  double sigmaRand =  sigmaSumSave * pState->rndm.flat();
  for (int i = 0; i < inPair.size(); ++i) 
  {
    sigmaRand -= inPair[i].pdfSigma;
    if (sigmaRand <= 0.) 
    {
      id1      = inPair[i].A.id;
      id2      = inPair[i].B.id;
      pdf1Save = inPair[i].A.pdf;
      pdf2Save = inPair[i].B.pdf;
      break;
    }
  }

}

// sets phi randomly and calculates pTFin, sets parton[0-3] and p1Res/p2Res
// still not clear what these actually are
bool SigmaProcess::final2KinMPI(int i1Res, int i2Res, Vec4ref p1Res, Vec4ref p2Res, double m1Res, double m2Res)
{
  // Perform kinematics for a multiparton interaction, including a rescattering.

  // TODO: what do you mean by 'perform kinematics'?
  // TODO: why do we need a final kinematics only for MPI?


  // Have to set flavours and colours.
  setIdColAcol();

  // Check that masses of outgoing particles not too big.
  m3           = pState->particleData.m0(idSave[3]);
  m4           = pState->particleData.m0(idSave[4]);
  mH           = sqrt(sH);
  if (m3 + m4 + MASSMARGIN > mH) return false;
  s3           = m3 * m3;
  s4           = m4 * m4;

  // Do kinematics of the production; without or with masses.
  double e1In  = 0.5 * mH;
  double e2In  = e1In;
  double pzIn  = e1In;
  if (i1Res > 0 || i2Res > 0) 
  {
    double s1  = m1Res * m1Res;
    double s2  = m2Res * m2Res;
    e1In       = 0.5 * (sH + s1 - s2) / mH;
    e2In       = 0.5 * (sH + s2 - s1) / mH;
    pzIn       = sqrtpos( e1In*e1In - s1 );
  }

  // Do kinematics of the decay.
  // Note: so we seem to be selecting the azmuth here?
  double e3    = 0.5 * (sH + s3 - s4) / mH;
  double e4    = 0.5 * (sH + s4 - s3) / mH;
  double pAbs  = sqrtpos( e3*e3 - s3 );
  phi          = 2. * M_PI * pState->rndm.flat();
  double pZ    = pAbs * cosTheta;
  pTFin        = pAbs * sinTheta;
  double pX    = pTFin * sin(phi);
  double pY    = pTFin * cos(phi);
  double scale = 0.5 * mH * sinTheta;

  // Fill particle info.
  // TODO: check the particle class. Is this usually used for the process catalog?
  int status1  = (i1Res == 0) ? -31 : -34;
  int status2  = (i2Res == 0) ? -31 : -34;
  parton[1]    = Particle( idSave[1], status1, 0, 0, 3, 4,
    colSave[1], acolSave[1],  0.,  0.,  pzIn, e1In, m1Res, scale);
  parton[2]    = Particle( idSave[2], status2, 0, 0, 3, 4,
    colSave[2], acolSave[2],  0.,  0., -pzIn, e2In, m2Res, scale);
  parton[3]    = Particle( idSave[3],      33, 1, 2, 0, 0,
    colSave[3], acolSave[3],  pX,  pY,    pZ,   e3,    m3, scale);
  parton[4]    = Particle( idSave[4],      33, 1, 2, 0, 0,
    colSave[4], acolSave[4], -pX, -pY,   -pZ,   e4,    m4, scale);

  // Boost particles from subprocess rest frame to event rest frame.
  // Normal multiparton interaction: only longitudinal boost.
  if (i1Res == 0 && i2Res == 0) {
    double betaZ = (x1Save - x2Save) / (x1Save + x2Save);
    for (int i = 1; i <= 4; ++i) parton[i].bst(0., 0., betaZ);
  // Rescattering: generic rotation and boost required.
  } else {
    RotBstMatrix M;
    M.fromCMframe( p1Res, p2Res);
    for (int i = 1; i <= 4; ++i) parton[i].rotbst(M);
  }

  // Done.
  return true;
}

// sets Q2RenSave and Q2FacSave for LHA processes
// (only used in PhaseSpace.h)
void SigmaProcess::setScale()
{
  // If scale has not been set, then to set.
  double scaleLHA = pState->lhaUp->scale();
  if (scaleLHA < 0.) 
  {
    // Final-state partons and their invariant mass.
    vector<int> iFin;
    Vec4 pFinSum;
    for (int i = 3; i < pState->lhaUp->sizePart(); ++i)
    if (pState->lhaUp->mother1(i) == 1) {
      iFin.push_back(i);
      pFinSum += Vec4( pState->lhaUp->px(i), pState->lhaUp->py(i),
        pState->lhaUp->pz(i), pState->lhaUp->e(i) );
    }
    int nFin = iFin.size();
    sH       = pFinSum * pFinSum;
    mH       = sqrt(sH);
    sH2      = sH * sH;

    // If 1 final-state particle then use Sigma1Process logic.
    if (nFin == 1) {
      Q2RenSave                             = renormMultFac * sH;
      if (renormScale1 == 2) Q2RenSave      = renormFixScale;
      Q2FacSave                             = factorMultFac * sH;
      if (factorScale1 == 2) Q2FacSave      = factorFixScale;

    // If 2 final-state particles then use Sigma2Process logic.
    } else if (nFin == 2) {
      double s3  = pow2(pState->lhaUp->m(iFin[0]));
      double s4  = pow2(pState->lhaUp->m(iFin[1]));
      double pT2 = pow2(pState->lhaUp->px(iFin[0])) + pow2(pState->lhaUp->py(iFin[0]));
      if      (renormScale2 == 1) Q2RenSave = pT2 + min(s3, s4);
      else if (renormScale2 == 2) Q2RenSave = sqrt((pT2 + s3) * (pT2 + s4));
      else if (renormScale2 == 3) Q2RenSave = pT2 + 0.5 * (s3 + s4);
      else                        Q2RenSave = sH;
      Q2RenSave                            *= renormMultFac;
      if      (renormScale2 == 5) Q2RenSave = renormFixScale;
      if      (factorScale2 == 1) Q2FacSave = pT2 + min(s3, s4);
      else if (factorScale2 == 2) Q2FacSave = sqrt((pT2 + s3) * (pT2 + s4));
      else if (factorScale2 == 3) Q2FacSave = pT2 + 0.5 * (s3 + s4);
      else                        Q2FacSave = sH;
      Q2FacSave                            *= factorMultFac;
      if      (factorScale2 == 5) Q2FacSave = factorFixScale;

    // If 3 or more final-state particles then use Sigma3Process logic.
    } else {
      double mTSlow  = sH;
      double mTSmed  = sH;
      double mTSprod = 1.;
      double mTSsum  = 0.;
      for (int i = 0; i < nFin; ++i) {
        double mTSnow = pow2(pState->lhaUp->m(iFin[i]))
          + pow2(pState->lhaUp->px(iFin[i])) + pow2(pState->lhaUp->py(iFin[i]));
        if      (mTSnow < mTSlow) {mTSmed = mTSlow; mTSlow = mTSnow;}
        else if (mTSnow < mTSmed) mTSmed = mTSnow;
        mTSprod *= mTSnow;
        mTSsum  += mTSnow;
      }
      if      (renormScale3 == 1) Q2RenSave = mTSlow;
      else if (renormScale3 == 2) Q2RenSave = sqrt(mTSlow * mTSmed);
      else if (renormScale3 == 3) Q2RenSave = pow(mTSprod, 1. / nFin);
      else if (renormScale3 == 4) Q2RenSave = mTSsum / nFin;
      else                        Q2RenSave = sH;
      Q2RenSave                            *= renormMultFac;
      if      (renormScale3 == 6) Q2RenSave = renormFixScale;
      if      (factorScale3 == 1) Q2FacSave = mTSlow;
      else if (factorScale3 == 2) Q2FacSave = sqrt(mTSlow * mTSmed);
      else if (factorScale3 == 3) Q2FacSave = pow(mTSprod, 1. / nFin);
      else if (factorScale3 == 4) Q2FacSave = mTSsum / nFin;
      else                        Q2FacSave = sH;
      Q2FacSave                            *= factorMultFac;
      if      (factorScale3 == 6) Q2FacSave = factorFixScale;
    }
  }

  // If alpha_strong and alpha_EM have not been set, then set them.
  if (pState->lhaUp->alphaQCD() < 0.001) {
    double Q2RenNow = (scaleLHA < 0.) ? Q2RenSave : pow2(scaleLHA);
    alpS = pState->couplings->alphaS(Q2RenNow);
  }
  if (pState->lhaUp->alphaQED() < 0.001) {
    double Q2RenNow = (scaleLHA < 0.) ? Q2RenSave : pow2(scaleLHA);
    alpEM = pState->couplings->alphaEM(Q2RenNow);
  }

}

// saves a copy of some kinematic variables (already in the class)
// (only used in MultiPartonInteractions.h)
void SigmaProcess::saveKin()
{
    for (int i = 0; i < 12; i++) 
    { 
      partonT[i] = parton[i];
      mSaveT[i] = mSave[i]; 
    }
    pTFinT = pTFin; phiT = phi; cosThetaT = cosTheta; sinThetaT = sinTheta;
}

// loads the kinematic variables that were saved
// (only used in MultiPartonInteractions.h)
void SigmaProcess::loadKin()
{
  for (int i = 0; i < 12; i++) 
  { 
    parton[i] = partonT[i];
    mSave[i] = mSaveT[i]; 
  }
  pTFin = pTFinT; cosTheta = cosThetaT; sinTheta = sinThetaT; phi = phiT;
}

// similar to load
// (only used in MultiPartonInteractions.h)
void SigmaProcess::swapKin()
{
  for (int i = 0; i < 12; i++)
  { 
    swap(parton[i], partonT[i]);
    swap(mSave[i], mSaveT[i]);
  }
  swap(pTFin, pTFinT); 
  swap(cosTheta, cosThetaT);
  swap(sinTheta, sinThetaT); 
  swap(phi, phiT);
}

// (unused)
// sets allowME, mME, pME. Not entirely sure what it is. Usually 'matrix elements' are presented as an alternative to 'parton showers' somehow
bool SigmaProcess::setupForMatrixElement()
{
  if (type == ProcessType::P2to3)
  {
    // Common initial-state handling.
      bool allowME = setupForMatrixElementInitial();

      // Correct outgoing c, b, mu and tau to be massive or not.
      mME[2] = m3; // Note: mass for matrix element??
      int id3Tmp = abs(id3Mass);
      if (id3Tmp ==  4) mME[2] = mcME;
      if (id3Tmp ==  5) mME[2] = mbME;
      if (id3Tmp == 13) mME[2] = mmuME;
      if (id3Tmp == 15) mME[2] = mtauME;
      mME[3] = m4;
      int id4Tmp = abs(id4Mass);
      if (id4Tmp ==  4) mME[3] = mcME;
      if (id4Tmp ==  5) mME[3] = mbME;
      if (id4Tmp == 13) mME[3] = mmuME;
      if (id4Tmp == 15) mME[3] = mtauME;
      mME[4] = m5;
      int id5Tmp = abs(id5Mass);
      if (id5Tmp ==  4) mME[4] = mcME;
      if (id5Tmp ==  5) mME[4] = mbME;
      if (id5Tmp == 13) mME[4] = mmuME;
      if (id5Tmp == 15) mME[4] = mtauME;

      // If kinematically impossible turn to massless case, but set error.
      if (mME[2] + mME[3] + mME[4] >= mH) {
        mME[2] = 0.;
        mME[3] = 0.;
        mME[4] = 0.;
        allowME = false;
      }

      // Form new average masses if identical particles.
      if (id3Tmp != 0 && id4Tmp == id3Tmp && id5Tmp == id3Tmp) {
        double mAvg = (mME[2] + mME[3] + mME[4]) / 3.;
        mME[2] = mAvg;
        mME[3] = mAvg;
        mME[4] = mAvg;
      } else if (id3Tmp != 0 && id4Tmp == id3Tmp) {
        mME[2] = sqrtpos(0.5 * (pow2(mME[2]) + pow2(mME[3]))
              - 0.25 * pow2(pow2(mME[2]) - pow2(mME[3])) / sH);
        mME[3] = mME[2];
      } else if (id3Tmp != 0 && id5Tmp == id3Tmp) {
        mME[2] = sqrtpos(0.5 * (pow2(mME[2]) + pow2(mME[4]))
              - 0.25 * pow2(pow2(mME[2]) - pow2(mME[4])) / sH);
        mME[4] = mME[2];
      } else if (id4Tmp != 0 && id5Tmp == id4Tmp) {
        mME[3] = sqrtpos(0.5 * (pow2(mME[3]) + pow2(mME[4]))
              - 0.25 * pow2(pow2(mME[3]) - pow2(mME[4])) / sH);
        mME[4] = mME[2];
      }

      // Iterate rescaled three-momenta until convergence.
      double m2ME3 = pow2(mME[2]);
      double m2ME4 = pow2(mME[3]);
      double m2ME5 = pow2(mME[4]);
      double p2ME3 = p3cm.pAbs2();
      double p2ME4 = p4cm.pAbs2();
      double p2ME5 = p5cm.pAbs2();
      double p2sum = p2ME3 + p2ME4 + p2ME5;
      double eME3  = sqrt(m2ME3 + p2ME3);
      double eME4  = sqrt(m2ME4 + p2ME4);
      double eME5  = sqrt(m2ME5 + p2ME5);
      double esum  = eME3 + eME4 + eME5;
      double p2rat = p2ME3 / eME3 + p2ME4 / eME4 + p2ME5 / eME5;
      int iStep = 0;
      while ( abs(esum - mH) > COMPRELERR * mH && iStep < NCOMPSTEP ) 
      {
        ++iStep;
        double compFac = 1. + 2. * (mH - esum) / p2rat;
        p2ME3 *= compFac;
        p2ME4 *= compFac;
        p2ME5 *= compFac;
        eME3   = sqrt(m2ME3 + p2ME3);
        eME4   = sqrt(m2ME4 + p2ME4);
        eME5   = sqrt(m2ME5 + p2ME5);
        esum   = eME3 + eME4 + eME5;
        p2rat  = p2ME3 / eME3 + p2ME4 / eME4 + p2ME5 / eME5;
      }

      // If failed convergence set error flag.
      if (abs(esum - mH) > COMPRELERR * mH) allowME = false;

      // Set up accepted kinematics.
      double totFac = sqrt( (p2ME3 + p2ME4 + p2ME5) / p2sum);
      pME[2] = totFac * p3cm;
      pME[2].e( eME3);
      pME[3] = totFac * p4cm;
      pME[3].e( eME4);
      pME[4] = totFac * p5cm;
      pME[4].e( eME5);

      // Done.
      return allowME;
  }

  if (type == ProcessType::P2to2)
  {
    // Common initial-state handling.
    bool allowME = setupForMatrixElementInitial();

    // Correct outgoing c, b, mu and tau to be massive or not.
    mME[2] = m3;
    int id3Tmp = abs(id3Mass);
    if (id3Tmp ==  4) mME[2] = mcME;
    if (id3Tmp ==  5) mME[2] = mbME;
    if (id3Tmp == 13) mME[2] = mmuME;
    if (id3Tmp == 15) mME[2] = mtauME;
    mME[3] = m4;
    int id4Tmp = abs(id4Mass);
    if (id4Tmp ==  4) mME[3] = mcME;
    if (id4Tmp ==  5) mME[3] = mbME;
    if (id4Tmp == 13) mME[3] = mmuME;
    if (id4Tmp == 15) mME[3] = mtauME;

    // If kinematically impossible turn to massless case, but set error.
    if (mME[2] + mME[3] >= mH) 
    {
      mME[2] = 0.;
      mME[3] = 0.;
      allowME = false;
    }

    // Calculate scattering angle in subsystem rest frame.
    double sH34 = sqrtpos( pow2(sH - s3 - s4) - 4. * s3 * s4);
    double cThe = (tH - uH) / sH34;
    double sThe = sqrtpos(1. - cThe * cThe);

    // Setup massive kinematics with preserved scattering angle.
    double s3ME   = pow2(mME[2]);
    double s4ME   = pow2(mME[3]);
    double sH34ME = sqrtpos( pow2(sH - s3ME - s4ME) - 4. * s3ME * s4ME);
    double pAbsME = 0.5 * sH34ME / mH;

    // Normally allowed with unequal (or vanishing) masses.
    if (id3Tmp == 0 || id3Tmp != id4Tmp) 
    {
      pME[2] = Vec4(  pAbsME * sThe, 0.,  pAbsME * cThe,
              0.5 * (sH + s3ME - s4ME) / mH);
      pME[3] = Vec4( -pAbsME * sThe, 0., -pAbsME * cThe,
              0.5 * (sH + s4ME - s3ME) / mH);

    // For equal (anti)particles (e.g. W+ W-) use averaged mass.
    } else {
      mME[2] = sqrtpos(0.5 * (s3ME + s4ME) - 0.25 * pow2(s3ME - s4ME) / sH);
      mME[3] = mME[2];
      pME[2] = Vec4(  pAbsME * sThe, 0.,  pAbsME * cThe, 0.5 * mH);
      pME[3] = Vec4( -pAbsME * sThe, 0., -pAbsME * cThe, 0.5 * mH);
    }

    // Done.
    return allowME;
  }

  // Note: change ME to MatrixElements
  // Note: change to setupMatrixElementFinalState
  if (type == ProcessType::P2to1)
  {
    // Common initial-state handling.
    bool allowME = setupForMatrixElementInitial();
    
    // Note: so we are just setting the final state 4-mom
    // Final state trivial here.
    mME[2] = mH;
    pME[2] = Vec4( 0., 0., 0., mH);

    // Done.
    return allowME;
  }

  // SigmaProcess
  return true;
}

// (helper function)
bool SigmaProcess::setupForMatrixElementInitial()
{
  // Initially assume it will work out to set up modified kinematics.
  bool allowME = true;

  // Note: id1, id2 are the incoming partons
  // Correct incoming c, b, mu and tau to be massive or not.
  mME[0] = 0.;
  int id1Tmp = abs(id1);
  if (id1Tmp ==  4) mME[0] = mcME; // Note: these are aliases from particleData
  if (id1Tmp ==  5) mME[0] = mbME;
  if (id1Tmp == 13) mME[0] = mmuME;
  if (id1Tmp == 15) mME[0] = mtauME;
  mME[1] = 0.;
  int id2Tmp = abs(id2);
  if (id2Tmp ==  4) mME[1] = mcME;
  if (id2Tmp ==  5) mME[1] = mbME;
  if (id2Tmp == 13) mME[1] = mmuME;
  if (id2Tmp == 15) mME[1] = mtauME;

  // If kinematically impossible return to massless case, but set error.
  if (mME[0] + mME[1] >= mH) // TODO: what is mH?
  {
    mME[0] = 0.;
    mME[1] = 0.;
    allowME = false;
  }

  // Do incoming two-body kinematics for massless or massive cases.
  // TODO: what is it? the 4-momentum of the incoming partons?
  if (mME[0] == 0. && mME[1] == 0.) 
  {
    pME[0] = 0.5 * mH * Vec4( 0., 0.,  1., 1.);
    pME[1] = 0.5 * mH * Vec4( 0., 0., -1., 1.);
  } 
  else 
  {
    double e0   = 0.5 * (mH * mH + mME[0] * mME[0] - mME[1] * mME[1]) / mH;
    double pz0  = sqrtpos(e0 * e0 - mME[0] * mME[0]);
    pME[0] = Vec4( 0., 0.,  pz0, e0);
    pME[1] = Vec4( 0., 0., -pz0, mH - e0);
  }

  // Done.
  return allowME;

}

void SigmaProcess::setId(int id1in, int id2in, int id3in, int id4in, int id5in)
{
  idSave[1] = id1in; idSave[2] = id2in; idSave[3] = id3in;
  idSave[4] = id4in; idSave[5] = id5in;
}

// set colors of particles 1 to 5. This should be deleted.
void SigmaProcess::setColAcol(int col1, int acol1, int col2, int acol2, int col3, int acol3, int col4, int acol4, int col5, int acol5)
{
  colSave[1] = col1; acolSave[1] = acol1; colSave[2] = col2;
  acolSave[2] = acol2; colSave[3] = col3; acolSave[3] = acol3;
  colSave[4] = col4; acolSave[4] = acol4; colSave[5] = col5;
  acolSave[5] = acol5;
}

// function for swapping colors/anticolors. Probably should be deleted
void SigmaProcess::swapColAcol()
{
  swap(colSave[1], acolSave[1]);
  swap(colSave[2], acolSave[2]); 
  swap(colSave[3], acolSave[3]);
  swap(colSave[4], acolSave[4]); 
  swap(colSave[5], acolSave[5]);
}

// function for swapping colors/anticolors. Probably should be deleted
void SigmaProcess::swapCol1234()
{
  swap(colSave[1], colSave[2]);
  swap(colSave[3], colSave[4]); 
  swap(acolSave[1], acolSave[2]);
  swap(acolSave[3], acolSave[4]);
}

// function for swapping colors/anticolors. Probably should be deleted
void SigmaProcess::swapCol12()
{
  swap(colSave[1], colSave[2]);
  swap(acolSave[1], acolSave[2]);
}

// function for swapping colors/anticolors. Probably should be deleted
void SigmaProcess::swapCol34()
{
  swap(colSave[3], colSave[4]);
  swap(acolSave[3], acolSave[4]);
}

// why do we need special weight functions for top/higgs?
double SigmaProcess::weightTopDecay(Event& process, int iResBeg, int iResEnd)
{
  // Note: iResBeg to iResEnd is a range of parton ids to test decays of

  // TODO: shouldn't there be another input parton along with t?
  // Evaluate weight for W decay distribution in t -> W b -> f fbar b.

  // If not pair W d/s/b and mother t then return unit weight.
  // Note: so it looks like this function expects 2 of them
  if (iResEnd - iResBeg != 1) return 1.;

  int iW1  = iResBeg; // id of W (W+ or W-)
  int iB2  = iResBeg + 1; // id of quark

  // Note: so these are not ids by process indices; we get the ids here
  int idW1 = process[iW1].idAbs();
  int idB2 = process[iB2].idAbs();

  // Note: make sure the ids are what we expect
  if (idW1 != 24)
  {
    swap(iW1, iB2);
    swap(idW1, idB2);
  }

  if (idW1 != 24 || (idB2 != 1 && idB2 != 3 && idB2 != 5)) return 1.;
  
  // Note: so the id of the W mother must be a top quark
  int iT   = process[iW1].mother1();
  if (iT <= 0 || process[iT].idAbs() != 6) return 1.;

  // Find sign-matched order of W decay products.
  int iF    = process[iW1].daughter1();
  int iFbar = process[iW1].daughter2();

  // Note: so the daughters of the W must be a ffbar pair that differ by 1 generation
  if (iFbar - iF != 1) return 1.;
  if (process[iT].id() * process[iF].id() < 0) swap(iF, iFbar);

  // Note: so the weight depends on the 4-mom of the W daughter, the t mother and the b
  // TODO: what is this weight and why is it set this way?
  // Weight and maximum weight.
  double wt    = (process[iT].p() * process[iFbar].p()) * (process[iF].p() * process[iB2].p());
  double wtMax = ( pow4(process[iT].m()) - pow4(process[iW1].m()) ) / 8.;

  // Done.
  return wt / wtMax;

}

double SigmaProcess::weightHiggsDecay(Event& process, int iResBeg, int iResEnd)
{
  // Note: see weightTopDecay as this is very similar

  // If not pair Z0 Z0, W+ W- or gamma Z0 then return unit weight.
  if (iResEnd - iResBeg != 1) return 1.;
  int iZW1  = iResBeg;
  int iZW2  = iResBeg + 1;
  int idZW1 = process[iZW1].id();
  int idZW2 = process[iZW2].id();
  if (idZW1 < 0 || idZW2 == 22) {
    swap(iZW1, iZW2);
    swap(idZW1, idZW2);
  }
  if ( (idZW1 != 23 || idZW2 != 23) && (idZW1 != 24 || idZW2 != -24)
    && (idZW1 != 22 || idZW2 != 23) ) return 1.;

  // If mother is not Higgs then return unit weight.
  int iH  = process[iZW1].mother1();
  if (iH <= 0) return 1.;
  int idH = process[iH].id();
  if (idH != 25 && idH != 35 && idH !=36) return 1.;

  // H -> gamma Z0 -> gamma f fbar is 1 + cos^2(theta) in Z rest frame.
  if (idZW1 == 22) {
    int i5 = process[iZW2].daughter1();
    int i6 = process[iZW2].daughter2();
    double pgmZ = process[iZW1].p() * process[iZW2].p();
    double pgm5 = process[iZW1].p() * process[i5].p();
    double pgm6 = process[iZW1].p() * process[i6].p();
    return (pow2(pgm5) + pow2(pgm6)) / pow2(pgmZ);
  }


  // CP violation parameters for the BSM Higgs sector.
  int higgsH1parity   = pState->settings.get(Mode::HiggsH1_parity);
  double higgsH1eta      = pState->settings.get(Param::HiggsH1_etaParity);
  double higgsH1phi      = pState->settings.get(Param::HiggsH1_phiParity);
  int higgsH2parity   = pState->settings.get(Mode::HiggsH2_parity);
  double higgsH2eta      = pState->settings.get(Param::HiggsH2_etaParity);
  double higgsH2phi      = pState->settings.get(Param::HiggsH2_phiParity);
  int higgsA3parity   = pState->settings.get(Mode::HiggsA3_parity);
  double higgsA3eta      = pState->settings.get(Param::HiggsA3_etaParity);
  double higgsA3phi      = pState->settings.get(Param::HiggsA3_phiParity);

  // If BSM not switched on then H1 should have SM properties.
  if (!pState->settings.get(Flag::Higgs_useBSM))
  {
    higgsH1parity = 1;
    higgsH1eta    = 0.;
    higgsH1phi    = M_PI / 2.;
  }

  // Parameters depend on Higgs type: H0(H_1), H^0(H_2) or A^0(H_3).
  int    higgsParity = higgsH1parity;
  double higgsEta    = higgsH1eta;
  if (idH == 35) {
    higgsParity      = higgsH2parity;
    higgsEta         = higgsH2eta;
  } else if (idH == 36) {
    higgsParity      = higgsA3parity;
    higgsEta         = higgsA3eta;
  }

  // Option with isotropic decays (also for pseudoscalar fermion couplings).
  if (higgsParity == 0 || higgsParity > 3) return 1.;

  // Maximum and initial weight.
  double wtMax = pow4(process[iH].m());
  double wt    = wtMax;

  // Find sign-matched order of Z0/W+- decay products.
  int i3 = process[iZW1].daughter1();
  int i4 = process[iZW1].daughter2();
  if (process[i3].id() < 0) swap( i3, i4);
  int i5 = process[iZW2].daughter1();
  int i6 = process[iZW2].daughter2();
  if (process[i5].id() < 0) swap( i5, i6);

  // Evaluate four-vector products and find masses..
  double p35  = 2. * process[i3].p() * process[i5].p();
  double p36  = 2. * process[i3].p() * process[i6].p();
  double p45  = 2. * process[i4].p() * process[i5].p();
  double p46  = 2. * process[i4].p() * process[i6].p();
  double p34  = 2. * process[i3].p() * process[i4].p();
  double p56  = 2. * process[i5].p() * process[i6].p();
  double mZW1 = process[iZW1].m();
  double mZW2 = process[iZW2].m();

  // For mixed CP states need epsilon product and gauge boson masses.
  double epsilonProd = 0.;
  if (higgsParity == 3) {
    double p[4][4];
    for (int i = 0; i < 4; ++i) {
      int         ii = i3;
      if (i == 1) ii = i4;
      if (i == 2) ii = i5;
      if (i == 3) ii = i6;
      p[i][0] = process[ii].e();
      p[i][1] = process[ii].px();
      p[i][2] = process[ii].py();
      p[i][3] = process[ii].pz();
    }
    epsilonProd
      = p[0][0]*p[1][1]*p[2][2]*p[3][3] - p[0][0]*p[1][1]*p[2][3]*p[3][2]
      - p[0][0]*p[1][2]*p[2][1]*p[3][3] + p[0][0]*p[1][2]*p[2][3]*p[3][1]
      + p[0][0]*p[1][3]*p[2][1]*p[3][2] - p[0][0]*p[1][3]*p[2][2]*p[3][1]
      - p[0][1]*p[1][0]*p[2][2]*p[3][3] + p[0][1]*p[1][0]*p[2][3]*p[3][2]
      + p[0][1]*p[1][2]*p[2][0]*p[3][3] - p[0][1]*p[1][2]*p[2][3]*p[3][0]
      - p[0][1]*p[1][3]*p[2][0]*p[3][2] + p[0][1]*p[1][3]*p[2][2]*p[3][0]
      + p[0][2]*p[1][0]*p[2][1]*p[3][3] - p[0][2]*p[1][0]*p[2][3]*p[3][1]
      - p[0][2]*p[1][1]*p[2][0]*p[3][3] + p[0][2]*p[1][1]*p[2][3]*p[3][0]
      + p[0][2]*p[1][3]*p[2][0]*p[3][1] - p[0][2]*p[1][3]*p[2][1]*p[3][0]
      - p[0][3]*p[1][0]*p[2][1]*p[3][2] + p[0][3]*p[1][0]*p[2][2]*p[3][1]
      + p[0][3]*p[1][1]*p[2][0]*p[3][2] - p[0][3]*p[1][1]*p[2][2]*p[3][0]
      - p[0][3]*p[1][2]*p[2][0]*p[3][1] + p[0][3]*p[1][2]*p[2][1]*p[3][0];
  }

  // Z0 Z0 decay: vector and axial couplings of two fermion pairs.
  if (idZW1 == 23) {
    double vf1 = pState->couplings->vf(process[i3].idAbs());
    double af1 = pState->couplings->af(process[i3].idAbs());
    double vf2 = pState->couplings->vf(process[i5].idAbs());
    double af2 = pState->couplings->af(process[i5].idAbs());
    double va12asym = 4. * vf1 * af1 * vf2 * af2
      / ( (vf1*vf1 + af1*af1) * (vf2*vf2 + af2*af2) );
    double vh = 1;
    double ah = higgsEta / pow2( pState->particleData.m0(23) );

    // Normal CP-even decay.
    if (higgsParity == 1) wt = 8. * (1. + va12asym) * p35 * p46
      + 8. * (1. - va12asym) * p36 * p45;

    // CP-odd decay (normal for A0(H_3)).
    else if (higgsParity == 2) wt = ( pow2(p35 + p46)
      + pow2(p36 + p45) - 2. * p34 * p56
      - 2. * pow2(p35 * p46 - p36 * p45) / (p34 * p56)
      + va12asym * (p35 + p36 - p45 - p46) * (p35 + p45 - p36 - p46) )
      / (1. +  va12asym);

    // Mixed CP states.
    else wt = 32. * ( 0.25 * pow2(vh) * ( (1. + va12asym) * p35 * p46
      + (1. - va12asym) * p36 * p45 ) - 0.5 * vh * ah * epsilonProd
      * ( (1. + va12asym) * (p35 + p46) - (1. - va12asym) * (p36 + p45) )
      + 0.0625 * pow2(ah) * (-2. * pow2(p34 * p56)
      - 2. * pow2(p35 * p46 - p36 * p45)
      + p34 * p56 * (pow2(p35 + p46) + pow2(p36 + p45))
      + va12asym * p34 * p56 * (p35 + p36 - p45 - p46)
      * (p35 + p45 - p36 - p46) ) )
      / ( pow2(vh) + 2. * abs(vh * ah) * mZW1 * mZW2
      + 2. * pow2(ah * mZW1 * mZW2) * (1. + va12asym) );

  // W+ W- decay.
  } else if (idZW1 == 24) {
    double vh = 1;
    double ah = higgsEta / pow2( pState->particleData.m0(24) );

    // Normal CP-even decay.
    if (higgsParity == 1) wt = 16. * p35 * p46;

    // CP-odd decay (normal for A0(H_3)).
    else if (higgsParity == 2) wt = 0.5 * ( pow2(p35 + p46)
      + pow2(p36 + p45) - 2. * p34 * p56
      - 2. * pow2(p35 * p46 - p36 * p45) / (p34 * p56)
      + (p35 + p36 - p45 - p46) * (p35 + p45 - p36 - p46) );

    // Mixed CP states.
    else wt = 32. * ( 0.25 * pow2(vh) * 2. * p35 * p46
      - 0.5 * vh * ah * epsilonProd * 2. * (p35 + p46)
      + 0.0625 * pow2(ah) * (-2. * pow2(p34 * p56)
      - 2. * pow2(p35 * p46 - p36 * p45)
      + p34 * p56 * (pow2(p35 + p46) + pow2(p36 + p45))
      + p34 * p56 * (p35 + p36 - p45 - p46) * (p35 + p45 - p36 - p46) ) )
      / ( pow2(vh) + 2. * abs(vh * ah) * mZW1 * mZW2
      + 2. * pow2(ah * mZW1 * mZW2) );
  }

  // Done.
  return wt / wtMax;

}

// --- virtual functions ---

// used for setting process id, name, resonance ids, etc
void SigmaProcess::initProc()
{
  // (empty)
}

void SigmaProcess::sigmaKin()
{
  // (empty)
}

// calculates the partonic cross-section for the process
double SigmaProcess::sigmaHat()
{
  // SigmaProcess
  return 0.;
}

// sets the colors dynamically based on the particles
// not sure why this needs to be specialized; couldn't we just use a look-up-table
// (only used in ProcessContainer::constructState)
void SigmaProcess::setIdColAcol()
{
  // (empty)
}

double SigmaProcess::weightDecayFlav(Event& process)
{
  // SigmaProcess
  return 1.0;
}

// what exactly is this?
double SigmaProcess::weightDecay(Event& process, int iResBeg, int iResEnd)
{
  if (type == ProcessType::LHA)
  {
    // Do nothing if decays present already at input.
    if (iResBeg < process.savedSizeValue()) return 1.;

    // Identity of mother of decaying reseonance(s).
    int idMother = process[process[iResBeg].mother1()].idAbs();

    // For Higgs decay hand over to standard routine.
    if (idMother == 25 || idMother == 35 || idMother == 36)
      return weightHiggsDecay( process, iResBeg, iResEnd);

    // For top decay hand over to standard routine.
    if (idMother == 6)
      return weightTopDecay( process, iResBeg, iResEnd);
  }

  // SigmaProcess
  return 1.;
}

} // end namespace Pythia