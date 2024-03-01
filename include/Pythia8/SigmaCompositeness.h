// SigmaCompositeness.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for compositiness-process differential cross sections.
// Contains classes derived from SigmaProcess via Sigma(1/2)Process.

#ifndef Pythia8_SigmaCompositeness_H
#define Pythia8_SigmaCompositeness_H

#include "Pythia8/SigmaProcess.h"

namespace Pythia8 {

//==========================================================================

// A derived class for q g -> q^* (excited quark state).

class Sigma1qg2qStar : public Sigma1Process {

public:

  // Constructor.
  Sigma1qg2qStar(int idqIn)
  {
    idq = idqIn;
    fluxType = FluxType::QG;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for q* decay angles (else inactive).
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Parameters set at initialization or for current kinematics.
  int    idq;
  double mRes, GammaRes, m2Res, GamMRat, Lambda, coupFcol, widthIn, sigBW;

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* qStarPtr;

};

//==========================================================================

// A derived class for l gamma -> l^* (excited lepton state).

class Sigma1lgm2lStar : public Sigma1Process {

public:

  // Constructor.
  Sigma1lgm2lStar(int idlIn) : idl(idlIn)
  {
    fluxType = FluxType::FGAMMA;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for l* decay angles (else inactive).
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Parameters set at initialization or for current kinematics.
  int    idl;
  double mRes, GammaRes, m2Res, GamMRat, Lambda, coupChg, widthIn, sigBW;

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* qStarPtr;

};

//==========================================================================

// A derived class for q q' -> q^* q' (excited quark state).

class Sigma2qq2qStarq : public Sigma2Process {

public:

  // Constructor.
  Sigma2qq2qStarq(int idqIn) : idq(idqIn)
  {
    fluxType = FluxType::QQ;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for q* decay angles (else inactive).
  virtual double weightDecay(Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.

private:

  // Parameters set at initialization or for current kinematics.
  int    idq;
  double Lambda, preFac, openFracPos, openFracNeg, sigmaA, sigmaB;

};

//==========================================================================

// A derived class for q qbar -> l^* lbar (excited lepton state).

class Sigma2qqbar2lStarlbar : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2lStarlbar(int idlIn) : idl(idlIn)
  {
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

  // Evaluate weight for l* decay angles (else inactive).
  virtual double weightDecay(Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.

private:

  // Parameters set at initialization or for current kinematics.
  int    idl;
  double Lambda, preFac, openFracPos, openFracNeg, sigma;

};

//==========================================================================

// A derived class for q qbar -> lStar lStarBar.
// Code contributed by Olga Igonkina.

class Sigma2qqbar2lStarlStarBar: public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2lStarlStarBar(int idlIn) : idl(idlIn)
  {
    fluxType = FluxType::QQBARSAME;
  }

  // Initialize process.
  void initProc();

  // Calculate flavour-independent parts of cross section.
  void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  void setIdColAcol();

  // Evaluate weight for l* decay angles (else inactive).
  virtual double weightDecay(Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.

private:

  // Parameters set at initialization or for current kinematics.
  int    idl, idRes;
  double Lambda, preFac, openFracPos, openFracNeg, sigma;

};

//==========================================================================

// A derived class for q q -> q q (quark contact interactions).
// Based on, Sigma2qq2qq (QCD).

class Sigma2QCqq2qq : public Sigma2Process {

public:

  // Constructor.
  Sigma2QCqq2qq()
  {
    fluxType = FluxType::QQ;
    name = "q q(bar)' -> (QC) -> q q(bar)'";
    code = 4201;
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
  double sigT, sigU, sigTU, sigST, sigSum, sigQCSTU, sigQCUTS;

  // Compositeness parameters.
  double qCLambda2;
  int    qCetaLL, qCetaRR, qCetaLR;

};

//==========================================================================

// A derived class for q qbar -> q' qbar' (quark contact interactions).
// Based on, Sigma2qqbar2qqbarNew(QCD).
// Note: This process give the same contributions for q == q' and q != q'.

class Sigma2QCqqbar2qqbar : public Sigma2Process {

public:

  // Constructor.
  Sigma2QCqqbar2qqbar() 
  {
    fluxType = FluxType::QQBARSAME;
    name = "q qbar -> (QC) -> q' qbar' (uds)";
    code = 4202;
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

  // Number of outgoing quark flavours to be considered, given that
  // matrix elements are calculated in the massless approximation.
  int    qCnQuarkNew;

  // Values stored for colour flow selection.
  int    idNew;
  double mNew, m2New, sigS, sigma;

  // Compositeness parameters.
  double qCLambda2;
  int    qCetaLL, qCetaRR, qCetaLR;

};

//==========================================================================

// A derived class for f fbar -> l lbar
// (contact interactions).
// Does not include t-channel contributions relevant for e^+e^- to e^+e^-

class Sigma2QCffbar2llbar : public Sigma2Process {

public:

  // Constructor: bool Graviton  = true, to use LED graviton settings.
  Sigma2QCffbar2llbar (int idIn, int codeIn) : idNew(idIn)
  {
    code = codeIn;
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

  virtual bool   isSChannel() const {return true;}

private:

  // Process values.
  int    idNew;
  double qCmNew, qCmNew2, qCmZ, qCmZ2, qCGZ, qCGZ2, sigma0;

  // Compositeness parameters.
  double qCLambda2;
  int    qCetaLL, qCetaRR, qCetaLR;
  double qCPropGm, qCrePropZ, qCimPropZ;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SigmaCompositeness_H
