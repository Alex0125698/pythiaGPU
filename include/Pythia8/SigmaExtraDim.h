// SigmaExtraDim.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Stefan Ask for the *LED* routines.
// Header file for extra-dimensional-process differential cross sections.
// Contains classes derived from SigmaProcess via Sigma(1/2)Process.

#ifndef Pythia8_SigmaExtraDim_H
#define Pythia8_SigmaExtraDim_H

#include "Pythia8/SigmaProcess.h"

namespace Pythia8 {

//==========================================================================

// A derived class for g g -> G^* (excited graviton state).

class Sigma1gg2GravitonStar : public SigmaProcess {

public:

  // Constructor.
  Sigma1gg2GravitonStar() : SigmaProcess(ProcessType::P2to1)
  {
    name = "g g -> G*";
    code = 5001;
    fluxType = FluxType::GG;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for G* decay angle.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Parameters set at initialization or for current kinematics.
  bool   eDsmbulk, eDvlvl;
  int    idGstar;
  double mRes, GammaRes, m2Res, GamMRat, kappaMG, sigma;

  // Couplings between graviton and SM (indexed by particle id).
  double eDcoupling[27];

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* gStarPtr;

};

//==========================================================================

// A derived class for f fbar -> G^* (excited graviton state).

class Sigma1ffbar2GravitonStar : public SigmaProcess {

public:

  // Constructor.
  Sigma1ffbar2GravitonStar() : SigmaProcess(ProcessType::P2to1)
  {
    name = "f fbar -> G*";
    code = 5002;
    fluxType = FluxType::FFBARSAME;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for G* decay angle.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Parameters set at initialization or for current kinematics.
  bool   eDsmbulk, eDvlvl;
  int    idGstar;
  double mRes, GammaRes, m2Res, GamMRat, kappaMG, sigma0;

  // Couplings between graviton and SM (indexed by particle id).
  double eDcoupling[27];

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* gStarPtr;

};

//==========================================================================

// A derived class for q qbar -> g^*/KK-gluon^* (excited kk-gluon state).

class Sigma1qqbar2KKgluonStar : public SigmaProcess {

public:

  // Constructor.
  Sigma1qqbar2KKgluonStar() : SigmaProcess(ProcessType::P2to1)
  {
    name = "q qbar -> g*/KK-gluon*";
    code = 5006;
    fluxType = FluxType::QQBARSAME;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for g* decay angle.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Parameters set at initialization or for current kinematics.
  int    idKKgluon;
  double mRes, GammaRes, m2Res, GamMRat;
  double sumSM, sumInt, sumKK, sigSM, sigInt, sigKK;

  // Couplings between kk gluon and SM (indexed by particle id).
  // Helicity dependent couplings. Use vector/axial-vector
  // couplings internally, gv/ga = 0.5 * (gL +/- gR).
  double eDgv[10], eDga[10];

  // Interference parameter.
  int interfMode;

  // Pointer to properties of the particle species, to access decay
  // channels.
  ParticleDataEntry* gStarPtr;

};

//==========================================================================

// A derived class for g g -> G^* g (excited graviton state).

class Sigma2gg2GravitonStarg : public SigmaProcess {

public:

  // Constructor.
  Sigma2gg2GravitonStarg() : SigmaProcess(ProcessType::P2to2)
  {
    name = "g g -> G* g";
    code = 5003;
    fluxType = FluxType::GG;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight: currently isotropic (except secondary top decay)..
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Parameters set at initialization or for current kinematics.
  int    idGstar;
  double mRes, GammaRes, m2Res, GamMRat, kappaMG, openFrac, sigma;

};

//==========================================================================

// A derived class for q g -> G^* q (excited graviton state).

class Sigma2qg2GravitonStarq : public SigmaProcess {

public:

  // Constructor.
  Sigma2qg2GravitonStarq() : SigmaProcess(ProcessType::P2to2)
  {
    name = "q g -> G* q";
    code = 5004;
    fluxType = FluxType::QG;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight: currently isotropic (except secondary top decay).
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Parameters set at initialization or for current kinematics.
  int    idGstar;
  double mRes, GammaRes, m2Res, GamMRat, kappaMG, openFrac, sigma;

};

//==========================================================================

// A derived class for q qbar -> G^* g (excited graviton state).

class Sigma2qqbar2GravitonStarg : public SigmaProcess {

public:

  // Constructor.
  Sigma2qqbar2GravitonStarg() : SigmaProcess(ProcessType::P2to2)
  {
    name = "q qbar -> G* g";
    code = 5005;
    fluxType = FluxType::QQBARSAME;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight: currently isotropic (except secondary top decay).
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Parameters set at initialization or for current kinematics.
  int    idGstar;
  double mRes, GammaRes, m2Res, GamMRat, kappaMG, openFrac, sigma;

};

//==========================================================================

// NOAM: A derived class for, f fbar -> (gamma/Z)_KKTower -> F Fbar,
// for one heavy F.
// Process provided by N. Hod et al. and is described in arXiv:XXXX.YYYY

class Sigma2ffbar2TEVffbar : public SigmaProcess {

public:

  // Constructor.
  Sigma2ffbar2TEVffbar(int idIn, int codeIn) : SigmaProcess(ProcessType::P2to2), idNew(idIn)
  {
    code = codeIn;
    fluxType = FluxType::FFBARSAME;
    isSChannel = true;
    idSChannel = 5000023;
    resonanceA = 23;
    resonanceB = 5000023;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for W decay angles in top decay (else inactive).
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Add phase-space sampling also around the Z_KK resonance.
  // TODO: surely this must be a bug...
  virtual int    phaseResonanceA();
  virtual int    phaseResonanceB();

private:

  // Values stored for process type.
  string  nameSave;
  int     idNew, gmZmode, nexcitationmax;
  bool    isPhysical;
  double  gPlusf, gMinusf, gPlusF, gMinusF, gPlusTop, gMinusTop, gf, gF;
  double  mRes, m2Res, mStar, mTop, m2Top, mZKKn, m2ZKKn, m2gmKKn, mgmKKn,
          alphaemfixed;
  double  helicityME2, coefTot, coefAngular;
  double  mr, betaf, cosThe, openFracPair;
  double  wgmKKFactor, wgmKKn, wZKKn,
          wZ0, ttbarwZKKn, ttbarwgmKKn,
          ttbarwFactorA, ttbarwFactorB;
  double  phaseSpacemHatMin, phaseSpacemHatMax;
  complex gammaProp, resProp, gmPropKK, ZPropKK, totalProp;
  complex mI;
};

//==========================================================================

// A derived class for g g -> U/G g (real graviton emission in
// large extra dimensions or unparticle emission).

class Sigma2gg2LEDUnparticleg : public SigmaProcess {

public:

  // Constructor: bool Graviton  = true, to use LED graviton settings.
  Sigma2gg2LEDUnparticleg( bool Graviton ) : SigmaProcess(ProcessType::P2to2), eDgraviton(Graviton)
  {
    name = eDgraviton ? "g g -> G g" : "g g -> U g";
    code = eDgraviton ? 5021 : 5045;
    fluxType = FluxType::GG;
    id3Mass = 5000039;
    id4Mass = 21;

  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section;
  // first step when inflavours unknown.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat); second step for given inflavours.
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  bool   eDgraviton;
  int    eDspin, eDnGrav, eDidG, eDcutoff;
  double mG, mGS, eDsigma0, eDdU, eDLambdaU, eDlambda, eDconstantTerm,
         eDtff, eDcf;

};

//==========================================================================

// A derived class for q g -> U/G q (real graviton emission in
// large extra dimensions or unparticle emission).

class Sigma2qg2LEDUnparticleq : public SigmaProcess {

public:

  // Constructor: bool Graviton  = true, to use LED graviton settings.
  Sigma2qg2LEDUnparticleq( bool Graviton) : SigmaProcess(ProcessType::P2to2), eDgraviton(Graviton)
  {
    name = eDgraviton ? "q g -> G q" : "q g -> U q";
    code = eDgraviton ? 5022 : 5046;
    fluxType = FluxType::QG;
    id3Mass = 5000039;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section;
  // first step when inflavours unknown.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat); second step for given inflavours.
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  bool   eDgraviton;
  int    eDspin, eDnGrav, eDidG, eDcutoff;
  double mG, mGS, eDsigma0, eDdU, eDLambdaU, eDlambda, eDconstantTerm,
         eDtff, eDgf, eDcf;

};

//==========================================================================

// A derived class for q qbar -> U/G g (real graviton emission in
// large extra dimensions or unparticle emission).

class Sigma2qqbar2LEDUnparticleg : public SigmaProcess {

public:

  // Constructor: bool Graviton  = true, to use LED graviton settings.
  Sigma2qqbar2LEDUnparticleg( bool Graviton) : SigmaProcess(ProcessType::P2to2), eDgraviton(Graviton)
  {
    name = eDgraviton ? "q qbar -> G g" : "q qbar -> U g";
    code = eDgraviton ? 5023 : 5047;
    fluxType = FluxType::QQBARSAME;
    id3Mass = 5000039;
    id4Mass = 21;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section;
  // first step when inflavours unknown.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat); second step for given inflavours.
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  bool   eDgraviton;
  int    eDspin, eDnGrav, eDidG, eDcutoff;
  double mG, mGS, eDsigma0, eDdU, eDLambdaU, eDlambda, eDconstantTerm,
         eDtff, eDgf, eDcf;

};

//==========================================================================

// A derived class for f fbar -> U/G Z (real graviton emission in
// large extra dimensions or unparticle emission).

class Sigma2ffbar2LEDUnparticleZ : public SigmaProcess {

public:

  // Constructor: bool Graviton  = true, to use LED graviton settings.
  Sigma2ffbar2LEDUnparticleZ( bool Graviton) : SigmaProcess(ProcessType::P2to2), eDgraviton(Graviton)
  {
    name = eDgraviton ? "f fbar -> G Z" : "f fbar -> U Z";
    code = eDgraviton ? 5024 : 5041;
    fluxType = FluxType::FFBARSAME;
    id3Mass = 5000039;
    id4Mass = 23;
    resonanceA = 23;
    gmZmode = 2;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section;
  // first step when inflavours unknown.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat); second step for given inflavours.
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Constants: could only be changed in the code itself.
  static const double FIXRATIO;

  int    eDspin, eDnGrav, eDcutoff, eDidG;
  bool   eDgraviton;
  double eDdU, eDLambdaU, eDlambda, eDratio, eDlambdaPrime,
         eDtff, eDconstantTerm;
  double sHS, tHS, uHS, tHC, uHC, tHQ, uHQ, tHuH, mU, mUS, mZ, widZ,
         mZS, mwZS, eDsigma0;

};

//==========================================================================

// A derived class for f fbar -> U/G gamma (real graviton emission in
// large extra dimensions or unparticle emission).

class Sigma2ffbar2LEDUnparticlegamma : public SigmaProcess {

public:

  // Constructor: bool Graviton  = true, to use LED graviton settings.
  Sigma2ffbar2LEDUnparticlegamma( bool Graviton) : SigmaProcess(ProcessType::P2to2), eDgraviton(Graviton)
  {
    name = eDgraviton ? "f fbar -> G gamma" : "f fbar -> U gamma";
    code = eDgraviton ? 5025 : 5042;
    fluxType = FluxType::FFBARSAME;
    id3Mass = 5000039;
    id4Mass = 22;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section;
  // first step when inflavours unknown.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat); second step for given inflavours.
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Constants: could only be changed in the code itself.
  static const double FIXRATIO;

  int    eDspin, eDnGrav, eDcutoff, eDidG;
  bool   eDgraviton;
  double eDdU, eDLambdaU, eDlambda, eDratio, eDlambdaPrime,
         eDtff, eDconstantTerm;
  double sHS, tHS, uHS, tHC, uHC, tHQ, uHQ, tHuH, mU, mUS, mZ,
         mZS, eDsigma0;

};

//==========================================================================

// A derived class for f fbar -> (LED G*/U*) -> gamma gamma
// (virtual graviton/unparticle exchange).

class Sigma2ffbar2LEDgammagamma : public SigmaProcess {

public:

  // Constructor: bool Graviton  = true, to use LED graviton settings.
  Sigma2ffbar2LEDgammagamma( bool Graviton) : SigmaProcess(ProcessType::P2to2), eDgraviton(Graviton)
  {
    name = eDgraviton ? "f fbar -> (LED G*) -> gamma gamma" : "f fbar -> (U*) -> gamma gamma";
    code = eDgraviton ? 5026 : 5043;
    fluxType = FluxType::FFBARSAME;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section;
  // first step when inflavours unknown.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat); second step for given inflavours.
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  int    eDspin, eDcutoff, eDnGrav, eDnegInt;
  bool   eDgraviton;
  double eDdU, eDLambdaU, eDlambda, eDlambda2chi,
         eDterm1, eDterm2, eDterm3, eDtff;

};

//==========================================================================

// A derived class for g g -> (LED G*/U*) -> gamma gamma
// (virtual graviton/unparticle exchange).

class Sigma2gg2LEDgammagamma : public SigmaProcess {

public:

  // Constructor: bool Graviton  = true, to use LED graviton settings.
  Sigma2gg2LEDgammagamma( bool Graviton) : SigmaProcess(ProcessType::P2to2), eDgraviton(Graviton)
  {
    name = eDgraviton ? "g g -> (LED G*) -> gamma gamma" : "g g -> (U*) -> gamma gamma";
    code = eDgraviton ? 5027 : 5044;
    fluxType = FluxType::GG;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section;
  // first step when inflavours unknown.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat); second step for given inflavours.
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  int    eDspin, eDcutoff, eDnGrav;
  bool   eDgraviton;
  double eDdU, eDLambdaU, eDlambda, eDlambda2chi, eDsigma0, eDtff;

};

//==========================================================================

// A derived class for f fbar -> (LED G*/U*) -> l lbar
// (virtual graviton/unparticle exchange).
// Does not include t-channel contributions relevant for e^+e^- to e^+e^-

class Sigma2ffbar2LEDllbar : public SigmaProcess {

public:

  // Constructor: bool Graviton  = true, to use LED graviton settings.
  Sigma2ffbar2LEDllbar( bool Graviton) : SigmaProcess(ProcessType::P2to2), eDgraviton(Graviton)
  {
    name = eDgraviton ? "f fbar -> (LED G*) -> l l" : "f fbar -> (U*) -> l l";
    code = eDgraviton ? 5028 : 5048;
    fluxType = FluxType::FFBARSAME;
    isSChannel = true;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section;
  // first step when inflavours unknown.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat); second step for given inflavours.
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  int    eDspin, eDcutoff, eDnGrav,eDnxx, eDnxy, eDnegInt;
  bool   eDgraviton;
  double eDdU, eDLambdaU, eDlambda, eDlambda2chi, eDtff,
         eDmZ, eDmZS, eDGZ, eDGZS, eDabsMeU, eDdenomPropZ, eDrePropGamma,
         eDrePropZ, eDimPropZ, eDabsAS, eDreA, eDreABW, eDpoly1, eDpoly2,
         eDpoly3;

};

//==========================================================================

// A derived class for g g -> (LED G*/U*) -> l lbar
// (virtual graviton/unparticle exchange).

class Sigma2gg2LEDllbar : public SigmaProcess {

public:

  // Constructor: bool Graviton  = true, to use LED graviton settings.
  Sigma2gg2LEDllbar( bool Graviton) : SigmaProcess(ProcessType::P2to2), eDgraviton(Graviton)
  {
    name = eDgraviton ? "g g -> (LED G*) -> l l" : "g g -> (U*) -> l l";
    code = eDgraviton ? 5029 : 5049;
    fluxType = FluxType::GG;
  }


  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section;
  // first step when inflavours unknown.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat); second step for given inflavours.
  virtual double sigmaHat() {return eDsigma0;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  int    eDspin, eDcutoff, eDnGrav;
  bool   eDgraviton;
  double eDdU, eDLambdaU, eDlambda, eDlambda2chi, eDsigma0, eDtff;

};

//==========================================================================

// A derived class for g g -> (LED G*) -> g g.

class Sigma2gg2LEDgg : public SigmaProcess {

public:

  // Constructor.
  Sigma2gg2LEDgg() : SigmaProcess(ProcessType::P2to2)
  {
    name = "g g -> (LED G*) -> g g";
    code = 5030;
    fluxType = FluxType::GG;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Values stored for colour flow selection.
  double sigTS, sigUS, sigTU, sigSum, sigma;

  // Model parameters.
  int eDopMode, eDnGrav, eDcutoff, eDnegInt;
  double eDMD, eDLambdaT, eDtff;

};

//==========================================================================

// A derived class for g g -> (LED G*) -> q qbar.

class Sigma2gg2LEDqqbar : public SigmaProcess {

public:

  // Constructor.
  Sigma2gg2LEDqqbar() : SigmaProcess(ProcessType::P2to2)
  {
    name = "g g -> (LED G*) -> q qbar (uds)";
    code = 5031;
    fluxType = FluxType::GG;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Number of quarks to be considered in massless approximation.
  int    nQuarkNew;

  // Values stored for colour flow selection.
  int    idNew;
  double mNew, m2New, sigTS, sigUS, sigSum, sigma;

  // Model parameters.
  int eDopMode, eDnGrav, eDcutoff, eDnegInt;
  double eDMD, eDLambdaT, eDtff;

};

//==========================================================================

// A derived class for q g -> (LED G*) -> q g.
// Use massless approximation also for Q since no alternative.

class Sigma2qg2LEDqg : public SigmaProcess {

public:

  // Constructor.
  Sigma2qg2LEDqg() : SigmaProcess(ProcessType::P2to2)
  {
    name = "q g -> (LED G*) -> q g";
    code = 5032;
    fluxType = FluxType::QG;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Values stored for colour flow selection.
  double sigTS, sigTU, sigSum, sigma;

  // Model parameters.
  int eDopMode, eDnGrav, eDcutoff, eDnegInt;
  double eDMD, eDLambdaT, eDtff;

};

//==========================================================================

// A derived class for q q(bar)' -> (LED G*) -> q q(bar)'.

class Sigma2qq2LEDqq : public SigmaProcess {

public:

  // Constructor.
  Sigma2qq2LEDqq() : SigmaProcess(ProcessType::P2to2)
  {
    name = "q q(bar)' -> (LED G*) -> q q(bar)'";
    code = 5033;
    fluxType = FluxType::QQ;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

 private:

  // Values stored for colour flow selection.
  double sigT, sigU, sigTU, sigST, sigSum;
  double sigGrT1, sigGrT2, sigGrU, sigGrTU, sigGrST;

  // Model parameters.
  int eDopMode, eDnGrav, eDcutoff, eDnegInt;
  double eDMD, eDLambdaT, eDtff;

};

//==========================================================================

// A derived class for q qbar -> (LED G*) -> g g.

class Sigma2qqbar2LEDgg : public SigmaProcess {

public:

  // Constructor.
  Sigma2qqbar2LEDgg() : SigmaProcess(ProcessType::P2to2)
  {
    name = "q qbar -> (LED G*) -> g g";
    code = 5034;
    fluxType = FluxType::QQBARSAME;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

 private:

  // Values stored for colour flow selection.
  double sigTS, sigUS, sigSum, sigma;

  // Model parameters.
  int eDopMode, eDnGrav, eDcutoff, eDnegInt;
  double eDMD, eDLambdaT, eDtff;

};

//==========================================================================

// A derived class for q qbar -> (LED G*) -> q' qbar'.

class Sigma2qqbar2LEDqqbarNew : public SigmaProcess {

public:

  // Constructor.
  Sigma2qqbar2LEDqqbarNew() : SigmaProcess(ProcessType::P2to2)
  {
    name = "q qbar -> (LED G*) -> q' qbar' (uds)";
    code = 5035;
    fluxType = FluxType::QQBARSAME;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

 private:

  // Number of quarks to be considered in massless approximation.
  int    nQuarkNew;

  // Values stored for colour flow selection.
  int    idNew;
  double mNew, m2New, sigS, sigma;

  // Model parameters.
  int eDopMode, eDnGrav, eDcutoff;
  double eDMD, eDLambdaT, eDtff;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SigmaExtraDim_H
