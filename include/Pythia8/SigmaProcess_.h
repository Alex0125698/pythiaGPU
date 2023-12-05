// SigmaProcess.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for hard-process differential cross sections.
// SigmaProcess: base class for cross sections.
// Sigma0Process: base class for unresolved processes, derived from above.
// Sigma1Process: base class for 2 -> 1 processes, derived from above.
// Sigma2Process: base class for 2 -> 2 processes, derived from above.
// Sigma3Process: base class for 2 -> 3 processes, derived from above.
// SigmaLHAProcess: wrapper class for Les Houches Accord external input.
// Actual physics processes are found in separate files:
// SigmaQCD for QCD processes;
// SigmaEW for electroweak processes (including photon production);
// SigmaOnia for charmonium and bottomonium processes;
// SigmaHiggs for Higgs processes;
// SigmaSUSY for supersymmetric production;
// SigmaLeftRightSym for  processes in left-right-symmetric scenarios;
// SigmaLeptoquark for leptoquark production.
// SigmaExtraDim for processes in scenarios for extra dimensions;
// SigmaGeneric for generic scalar/fermion/vector pair production.

#ifndef Pythia8_SigmaProcess_H
#define Pythia8_SigmaProcess_H

#include "Pythia8/Basics.h"
#include "Pythia8/BeamParticle.h"
#include "Pythia8/Event.h"
#include "Pythia8/Info.h"
#include "Pythia8/LesHouches.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PartonDistributions.h"
#include "Pythia8/PythiaComplex.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/ResonanceWidths.h"
#include "Pythia8/Settings.h"
#include "Pythia8/SigmaTotal.h"
#include "Pythia8/StandardModel.h"
#include "Pythia8/SLHAinterface.h"
#include "Pythia8/SusyLesHouches.h"

namespace Pythia8 {

// replacing the flux type strings with an enum
// ?? incoming parton pair type ??
enum class FluxType { NONE, GG, QG, QQ, QQBAR, QQBARSAME, FF, FFBAR,  
  FFBARSAME, FFBARCHG, FGM, GGAMMA, GAMMAGAMMA };

// the standard PDG codes
// https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
enum PDG
{
  sbar = -3,
  ubar = -2,
  dbar = -1,
  // skip 0
  d = 1,
  u = 2,
  s = 3,
  c = 4,
  b = 5,
  t = 6,
  bPrime = 7,
  tPrime = 8,
  // 9 previously usef for gluon
  eMinus = 11,
  nuE = 12,
  muMinus = 13,
  nuMu = 14,
  tauMinus = 15,
  nuTau = 16,
  tauMinusPrime = 17,
  nuTauPrime = 18,
  // 19,29 skipped
  g = 21,
  gamma = 22,
  Z = 23,
  Wplus = 24,
  h = 25,
  // skip 26 to 31
  Zprime = 32, Z2 = 32,
  ZprimePrime = 33, Z3 = 33,
  Wprime = 34, Wplus2 = 34,
  H2 = 35,
  A = 36, H3 = 36,
  Hplus = 37,
};

// identified each particle in the theory
enum class PID
{
  d, // down quark; PDG 1
  u, // up-quark; PDG 2
  s,
  c,
  b,
  t,
  bPrime,
  tPrime,
  eMinus,
  nuE,
  muMinus,
  nuMu,
  TauMinus,
  MuTau,
  g,
  gamma,
  Z,
  Wplus,
  h,
};

// add this to Pythia

struct PythiaState
{
  // true singletons
  Settings settings;
  ParticleData particleData;
  // duplicated per pythia event generator
  Info info;
  Rndm rndm;
  BeamParticle beamA;
  BeamParticle beamB;
  Couplings* couplings;
  SigmaTotal sigmaTotal;
  SLHAinterface slhaInterface;
  LHAup* lhaUp;
};


// Need the following globals

// Info, Settings, ParticleData, Rndm, Couplings, SigmaTotal, SLHAinterface

// a parton in the beam and its PDF
struct InBeam { PDG id; double pdf; };

// a parton pair in the beam and their pdfSigma
struct InPair { InBeam A; InBeam B; double pdfSigma; };

// base class for all processes; derived classes represent a particular cross-section
class SigmaProcess
{

public: // data

  string name; // process name
  FluxType fluxType         = FluxType::NONE;
  PythiaState* pState       = nullptr; // pointer to global Pythia state
  bool   convert2mb         = true;
  bool   convertM2          = false;
  bool   isLHA              = false; // is this an LHA process?
  bool   isNonDiff          = false; // 
  bool   isResolved         = true;
  bool   isDiffA            = false;
  bool   isDiffB            = false;
  bool   isDiffC            = false;
  bool   isSUSY             = false;
  bool   allowNegativeSigma = false;
  bool   isSChannel         = false;
  bool   isQCD3body         = false;
  bool   useMirrorWeight    = false;
  int    code               = 0;
  int    nFinal             = 2;
  int    gmZmode            = -1;
  int    id3Mass            = 0;
  int    id4Mass            = 0;
  int    id5Mass            = 0;
  int    resonanceA         = 0;
  int    resonanceB         = 0;
  int    idSChannel         = 0;
  int    idTchan1           = 0;
  int    idTchan2           = 0;
  double tChanFracPow1      = 0.3;
  double tChanFracPow2      = 0.3;
  
  // WARNING: skipping shorthand beam data
  // WARNING: skipping shorthand settings data
  // WARNING: need to override higgsH1parity, higgsH1eta, higgsH1phi for SM
  // WARNING: mass for mcME, mbME, mmuME, mtauME only considered if appropriate setting is on
  // WARNING: cacheched copy of beam properties may be 0 if beams don't exist
  // WARNING: nFinal seems to be 2 for Sigma0Process


  // TODO: sort out what these actually are

  // Constants: could only be changed in the code itself.
  static const double CONVERT2MB, MASSMARGIN, COMPRELERR;
  static const int    NCOMPSTEP;

  // Initialization data, normally only set once.
  int    nQuarkIn, renormScale1, renormScale2, renormScale3, renormScale3VV,
         factorScale1, factorScale2, factorScale3, factorScale3VV;
  double Kfactor, mcME, mbME, mmuME, mtauME, renormMultFac, renormFixScale,
         factorMultFac, factorFixScale;

  // CP violation parameters for Higgs sector, normally only set once.
  int    higgsH1parity, higgsH2parity, higgsA3parity;
  double higgsH1eta, higgsH2eta, higgsA3eta, higgsH1phi, higgsH2phi,
         higgsA3phi;

  // Information on incoming beams.
  int    idA, idB;
  double mA, mB;
  bool   isLeptonA, isLeptonB, hasLeptonBeams;

  // Partons in beams, with PDF's.
  vector<InBeam> inBeamA;
  vector<InBeam> inBeamB;

  // Allowed colliding parton pairs, with pdf's.
  vector<InPair> inPair;

  // Store common subprocess kinematics quantities.
  double mH, sH, sH2;

  // Store Q2 renormalization and factorization scales, and related values.
  double x1Save, x2Save; // incoming parton momentum fractions
  double alpEM, alpS; // electromagnetic and strong coupling constants
  double Q2RenSave, Q2FacSave, pdf1Save, pdf2Save, sigmaSumSave;

  // Store flavour, colour, anticolour, mass, angles and the whole particle.
  int      id1, id2, id3, id4, id5;
  int      idSave[12], colSave[12], acolSave[12];
  double   mSave[12], cosTheta, sinTheta, phi, sHMass, sHBeta, pT2Mass, pTFin;
  Particle parton[12];

  // Minimal set of saved kinematics for trial interactions when
  // using the x-dependent matter profile of multiparton interactions.
  Particle partonT[12];
  double   mSaveT[12], pTFinT, cosThetaT, sinThetaT, phiT;

  // Calculate and store all modified masses and four-vectors
  // intended for matrix elements. Return false if failed.
  double   mME[12];
  Vec4     pME[12];

  // Store whether tHat and uHat are swapped (= same as swap 3 and 4).
  bool swapTU;


public: // functions

  // set the PythiaState pointer
  void init(PythiaState* pState_)
  {
    pState = pState_;
  }

  // init process and allowed flux???
  bool setup()
  {
    // set the incoming particles list based on the flux type

    auto add = [&](PDG idInA, PDG idInB)
    {
      inBeamA.push_back({idInA,0.});
      inBeamA.push_back({idInB,0.});
      inPair.push_back({{idInA,0.},{idInB,0.},0.});
    };

    switch (fluxType)
    {
      case FluxType::NONE:
        pState->info.errorMsg("Error in SigmaProcess::initFlux; fluxType was not specified");
        return false;  

      case FluxType::GG:
        add(PDG::g, PDG::g);
        break;

      case FluxType::QG:
        
        break;

      case FluxType::QQ:
        
        break;

      case FluxType::QQBAR:
        
        break;

      case FluxType::QQBARSAME:
        
        break;

      case FluxType::FF:
        
        break;

      case FluxType::FFBAR:
        
        break;

      case FluxType::FFBARSAME:
        
        break;

      case FluxType::FFBARCHG:
        
        break;

      case FluxType::FGM:
        
        break;

      case FluxType::GGAMMA:
        
        break;

      case FluxType::GAMMAGAMMA:
        
        break;
    };

  }

  // TODO: does this actually change for each derived class? ->NO
  // Calculate modified masses and four-vectors for matrix elements.
  bool setupForMatrixElement();

  // TODO: does this actually change for each derived class? -> NO
  // Store kinematics and set scales for resolved 2 -> 1 process.
  void store1Kin( double x1in, double x2in, double sHin)
  {
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
    alpS   = couplingsPtr->alphaS(Q2RenSave);
    alpEM  = couplingsPtr->alphaEM(Q2RenSave);
  }

  // Store kinematics and set scales for resolved 2 -> 2 process.
  void store2Kin( double x1in, double x2in, double sHin,
    double tHin, double m3in, double m4in, double runBW3in,
    double runBW4in);

  void store2KinMPI( double x1in, double x2in, double sHin,
    double tHin, double uHin, double alpSin, double alpEMin,
    bool needMasses, double m3in, double m4in);

  // Store kinematics and set scales for resolved 2 -> 3 process.
  void store3Kin( double x1in, double x2in, double sHin,
    Vec4ref p3cmIn, Vec4ref p4cmIn, Vec4ref p5cmIn, double m3in, double m4in,
    double m5in, double runBW3in, double runBW4in, double runBW5in)
  {
    // Default ordering of particles 3 and 4 - not relevant here.
    swapTU   = false;

    // Incoming parton momentum fractions.
    x1Save   = x1in;
    x2Save   = x2in;

    // Incoming masses and their squares.
    if (id3Mass() == 0 && id4Mass() == 0 && id5Mass() == 0) {
      m3     = 0.;
      m4     = 0.;
      m5     = 0.;
    } else {
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
    if (isSChannel()) {

      // Different options for renormalization scale, but normally sHat.
      Q2RenSave = renormMultFac * sH;
      if (renormScale1 == 2) Q2RenSave = renormFixScale;

      // Different options for factorization scale, but normally sHat.
      Q2FacSave = factorMultFac * sH;
      if (factorScale1 == 2) Q2RenSave = factorFixScale;

    // "Normal" 2 -> 3 processes, i.e. not vector boson fusion.
    } else if ( idTchan1() != 23 && idTchan1() != 24 && idTchan2() != 23
      && idTchan2() != 24 ) {
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
    } else {
      double sV4   = pow2( particleDataPtr->m0(idTchan1()) );
      double sV5   = pow2( particleDataPtr->m0(idTchan2()) );
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
    alpS  = couplingsPtr->alphaS(Q2RenSave);
    alpEM = couplingsPtr->alphaEM(Q2RenSave);

  }


  // set kinematics for a 2 -> 1 process
  void set1Kin(double x1in, double x2in, double sHin)
  {
    store1Kin( x1in, x2in, sHin);
    sigmaKin();    
  }

  // set kinematics for a 2 -> 2 process
  void set2Kin(double x1in, double x2in, double sHin, double tHin, double m3in, double m4in, double runBW3in, double runBW4in)
  {
    store2Kin( x1in, x2in, sHin, tHin, m3in, m4in, runBW3in, runBW4in); 
    sigmaKin();
  }

  // set kinematics for a 2 -> 2 process (Multiparton Interaction)
  void set2KinMPI(double x1in, double x2in, double sHin, double tHin,double uHin, double 
                      alpSin, double alpEMin, bool needMasses, double m3in, double m4in)
  {

  }

  // set kinematics for a 2 -> 3 process
  void set3Kin( double x1in, double x2in, double sHin, Vec4 p3prel, Vec4 p4prel, Vec4 p5prel,
            double m3in, double m4in, double m5in, double runBW3in, double runBW4in, double runBW5in)
  {
    store3Kin( x1in, x2in, sHin, p3cmIn, p4cmIn, p5cmIn, m3in, m4in, m5in, runBW3in, runBW4in, runBW5in);
    sigmaKin();
  }

  // get the partonic cross-section
  virtual double sigmaHat()
  {
    // TODO: need to incorporate this somehow
    // Maybe have a constant that changes for each base class?
    if (convert2mb) result *= CONVERT2MB;
  }

  // get the hadronic cross-section
  virtual double sigmaPDF()
  {
    // 0 particles: Since no PDF's there is no difference from above.
    if (nFinal == 0) return sigmaHat();
  }

};


} // end namespace Pythia8

#endif // Pythia8_SigmaProcess_H
