// SigmaEW.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for electroweak process differential cross sections.
// Contains classes derived from SigmaProcess via Sigma(1/2)Process.

#ifndef Pythia8_SigmaEW_H
#define Pythia8_SigmaEW_H

#include "Pythia8/PythiaComplex.h"
#include "Pythia8/SigmaProcess.h"

namespace Pythia8 {


//==========================================================================

// A derived class for q g -> q gamma (q = u, d, s, c, b).
// Use massless approximation also for Q since no alternative.

class Sigma2qg2qgamma : public Sigma2Process {

public:

  // Constructor.
  Sigma2qg2qgamma()
  {
    name = "q g -> q gamma (udscb)";
    code = 201;
    fluxType = FluxType::QG;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Values stored for later use.
  double sigUS, sigma0;

};

//==========================================================================

// A derived class for q qbar -> g gamma.

class Sigma2qqbar2ggamma : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2ggamma()
  {
    name = "q qbar -> g gamma";
    code = 202;
    fluxType = FluxType::QQBARSAME;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Values stored for later use.
  double sigma0;

};

//==========================================================================

// A derived class for g g -> g gamma.

class Sigma2gg2ggamma : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2ggamma()
  {
    name = "g g -> g gamma";
    code = 203;
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

  // Values stored for later use.
  double chargeSum, sigma;

};

//==========================================================================

// A derived class for f fbar -> gamma gamma.

class Sigma2ffbar2gammagamma : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2gammagamma()
  {
    name = "f fbar -> gamma gamma";
    code = 204;
    fluxType = FluxType::FFBARSAME;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Values stored for later use.
  double sigTU, sigma0;

};

//==========================================================================

// A derived class for g g -> gamma gamma.

class Sigma2gg2gammagamma : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2gammagamma()
  {
    name = "g g -> gamma gamma";
    code = 205;
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

  double charge2Sum, sigma;

};

//==========================================================================

// A derived class for f f' -> f f' via t-channel gamma*/Z0 exchange.

class Sigma2ff2fftgmZ : public Sigma2Process {

public:

  // Constructor.
  Sigma2ff2fftgmZ()
  {
    name = "f f' -> f f' (t-channel gamma*/Z0)";
    code = 211;
    fluxType = FluxType::FF;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  //  Z parameters for propagator.
  int    gmZmode;
  double mZ, mZS, thetaWRat, sigmagmgm, sigmagmZ, sigmaZZ;

};

//==========================================================================

// A derived class for f_1 f_2 -> f_3 f_4 via t-channel W+- exchange.

class Sigma2ff2fftW : public Sigma2Process {

public:

  // Constructor.
  Sigma2ff2fftW()
  {
    name = "f_1 f_2 -> f_3 f_4 (t-channel W+-)";
    code = 212;
    fluxType = FluxType::FF;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  //  W parameters for propagator.
  double mW, mWS, thetaWRat, sigma0;

};

//==========================================================================

// A derived class for q q' -> Q q" via t-channel W+- exchange.
// Related to Sigma2ff2fftW class, but with massive matrix elements.

class Sigma2qq2QqtW : public Sigma2Process {

public:

  // Constructor.
  Sigma2qq2QqtW(int idIn, int codeIn) : idNew(idIn)
  {
    fluxType = FluxType::FF;
    resonanceA = idIn;
    code = codeIn;
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

private:

  // Values stored for process type. W parameters for propagator.
  int    idNew;
  double mW, mWS, thetaWRat, sigma0, openFracPos, openFracNeg;

};

//==========================================================================

// A derived class for f fbar -> gamma*/Z0.

class Sigma1ffbar2gmZ : public Sigma1Process {

public:

  // Constructor.
  Sigma1ffbar2gmZ()
  {
    name = "f fbar -> gamma*/Z0";
    code = 221;
    fluxType = FluxType::FFBARSAME;
    resonanceA = 23;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for Z decay angle.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Parameters set at initialization or for each new event.
  int    gmZmode;
  double mRes, GammaRes, m2Res, GamMRat, thetaWRat,
         gamSum, intSum, resSum, gamProp, intProp, resProp;

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* particlePtr;

};

//==========================================================================

// A derived class for f fbar' -> W+-.

class Sigma1ffbar2W : public Sigma1Process {

public:

  // Constructor.
  Sigma1ffbar2W()
  {
    name = "f fbar' -> W+-";
    code = 222;
    fluxType = FluxType::FFBARCHG;
    resonanceA = 24;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for W decay angle.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Parameters set at initialization.
  double mRes, GammaRes, m2Res, GamMRat, thetaWRat, sigma0Pos, sigma0Neg;

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* particlePtr;

};

//==========================================================================

// A derived class for f fbar -> gamma* -> f' fbar', summed over light f'.
// Allows pT-ordered evolution for multiparton interactions.

class Sigma2ffbar2ffbarsgm : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2ffbarsgm()
  {
    name = "f fbar -> f' fbar' (s-channel gamma*)";
    code = 223;
    fluxType = FluxType::FFBARSAME;
    isSChannel = true;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Values stored for later use.
  int    idNew;
  double sigma0;

};

//==========================================================================

// A derived class for f fbar -> gamma*/Z0 -> f' fbar', summed over light f.

class Sigma2ffbar2ffbarsgmZ : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2ffbarsgmZ()
  {
    name = "f fbar -> f' fbar' (s-channel gamma*/Z0)";
    code = 224;
    fluxType = FluxType::FFBARSAME;
    isSChannel = true;
    idSChannel = 23;
    resonanceA = 23;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Parameters set at initialization or for each new event.
  int    gmZmode;
  double mRes, GammaRes, m2Res, GamMRat, thetaWRat, colQ,
         gamSumT, gamSumL, intSumT, intSumL, intSumA, resSumT, resSumL,
         resSumA, gamProp, intProp, resProp, cThe;
  vector<int> idVec;
  vector<double> gamT, gamL, intT, intL, intA, resT, resL, resA, sigTLA;

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* particlePtr;

};

//==========================================================================

// A derived class for f_1 fbar_2 -> W+- -> f_3 fbar_4, summed over light f.

class Sigma2ffbar2ffbarsW : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2ffbarsW()
  {
    name = "f_1 fbar_2 -> f_3 fbar_4 (s-channel W+-)";
    code = 225;
    fluxType = FluxType::FFBARCHG;
    isSChannel = true;
    idSChannel = 24;
    resonanceA = 24;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Parameters set at initialization or stored for later use.
  int    id3New, id4New;
  double mRes, GammaRes, m2Res, GamMRat, thetaWRat, sigma0;

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* particlePtr;

};

//==========================================================================

// A derived class for f fbar -> gamma*/Z0 -> F Fbar, for one heavy F.
// Allows pT cuts as for other 2 -> 2 processes.

class Sigma2ffbar2FFbarsgmZ : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2FFbarsgmZ(int idIn, int codeIn) : idNew(idIn)
  {
    code = codeIn;
    fluxType = FluxType::FFBARSAME;
    isSChannel = true;
    id3Mass = idIn;
    id4Mass = idIn;
    resonanceA = 23;
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

private:

  // Values stored for process type. Z parameters for propagator.
  int    idNew, gmZmode;
  bool   isPhysical;
  double ef, vf, af, mRes, GammaRes, m2Res, GamMRat, thetaWRat,
         mr, betaf, cosThe, gamProp, intProp, resProp, openFracPair;

};

//==========================================================================

// A derived class for f fbar' -> W+- -> F fbar", for one or two heavy F.
// Allows pT cuts as for other 2 -> 2 processes.

class Sigma2ffbar2FfbarsW : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2FfbarsW(int idIn, int idIn2, int codeIn) : 
    idNew2(idIn2)
  {
    fluxType = FluxType::FFBARCHG;
    isSChannel = true;
    code = codeIn;
    resonanceA = 24;
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

private:

  // Values stored for process type. W parameters for propagator.
  int    idNew, idNew2, idPartner;
  bool   isPhysical;
  double V2New, mRes, GammaRes, m2Res, GamMRat, thetaWRat, sigma0,
         openFracPos, openFracNeg;

};

//==========================================================================

// An intermediate class for f fbar -> gamma*/Z0/W+- gamma*/Z0/W-+.

class Sigma2ffbargmZWgmZW : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbargmZWgmZW() {}

protected:

  // Internal products.
  Vec4    pRot[7];
  complex hA[7][7];
  complex hC[7][7];

  // Calculate and store internal products.
  void setupProd( Event& process, int i1, int i2, int i3, int i4,
    int i5, int i6);

  // Evaluate the F function of Gunion and Kunszt.
  complex fGK(int i1, int i2, int i3, int i4, int i5, int i6);

  // Evaluate the Xi function of Gunion and Kunszt.
  double xiGK( double tHnow, double uHnow);

  // Evaluate the Xj function of Gunion and Kunszt.
  double xjGK( double tHnow, double uHnow);

private:

};

//==========================================================================

// A derived class for f fbar -> gamma*/Z0 gamma*/Z0.

class Sigma2ffbar2gmZgmZ : public Sigma2ffbargmZWgmZW {

public:

  // Constructor.
  Sigma2ffbar2gmZgmZ()
  {
    name = "f fbar -> gamma*/Z0 gamma*/Z0";
    code = 231;
    fluxType = FluxType::FFBARSAME;
    id3Mass = 23;
    id4Mass = 23;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for simultaneous flavour choices.
  virtual double weightDecayFlav( Event& process);

  // Evaluate weight for decay angles of the two gamma*/Z0.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Parameters set at initialization or for each new event.
  int    gmZmode, i1, i2, i3, i4, i5, i6;
  double mRes, GammaRes, m2Res, GamMRat, thetaWRat, sigma0,
         gamSum3, intSum3, resSum3, gamProp3, intProp3, resProp3,
         gamSum4, intSum4, resSum4, gamProp4, intProp4, resProp4,
         c3LL, c3LR, c3RL, c3RR, c4LL, c4LR, c4RL, c4RR, flavWt;

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* particlePtr;

};

//==========================================================================

// A derived class for f fbar' -> Z0 W+-. (Here pure Z0, unfortunately.)

class Sigma2ffbar2ZW : public Sigma2ffbargmZWgmZW {

public:

  // Constructor.
  Sigma2ffbar2ZW()
  {
    name = "f fbar' -> Z0 W+- (no gamma*!)";
    code = 232;
    fluxType = FluxType::FFBARCHG;
    id3Mass = 23;
    id4Mass = 24;
    resonanceA = 24;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for Z0 and W+- decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Store W+- mass and width, and couplings.
  double mW, widW, mWS, mwWS, sin2thetaW, cos2thetaW, thetaWRat, cotT,
         thetaWpt, thetaWmm, lun, lde, sigma0, openFracPos, openFracNeg;

};

//==========================================================================

// A derived class for f fbar -> W+ W-.

class Sigma2ffbar2WW : public Sigma2ffbargmZWgmZW {

public:

  // Constructor.
  Sigma2ffbar2WW()
  {
    name = "f fbar -> W+ W-";
    code = 233;
    fluxType = FluxType::FFBARSAME;
    id3Mass = 24;
    id4Mass = -24;
    resonanceA = 23;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for W+ and W- decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Store Z0 mass and width.
  double mZ, widZ, mZS, mwZS, thetaWRat, sigma0, cgg, cgZ, cZZ, cfg,
    cfZ, cff, gSS, gTT, gST, gUU, gSU, openFracPair;

};

//==========================================================================

// An intermediate class for f fbar -> gamma*/Z0 g/gamma and permutations.

class Sigma2ffbargmZggm : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbargmZggm() {}

  // Initialize process.
  virtual void initProc();

  // Evaluate weight for gamma&/Z0 decay angle.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

protected:

  // Parameters set at initialization or for each new event.
  int    gmZmode;
  double mRes, GammaRes, m2Res, GamMRat, thetaWRat,
         gamSum, intSum, resSum, gamProp, intProp, resProp;

  // Evaluate current sum of flavour couplings times phase space.
  void flavSum();

  // Evaluate current propagator terms of cross section.
  void propTerm();

private:

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* particlePtr;

};

//==========================================================================

// A derived class for q qbar -> gamma*/Z0 g.

class Sigma2qqbar2gmZg : public Sigma2ffbargmZggm {

public:

  // Constructor.
  Sigma2qqbar2gmZg()
  {
    name = "q qbar -> gamma*/Z0 g";
    code = 241;
    fluxType = FluxType::QQBARSAME;
    id3Mass = 23;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Values stored for later use.
  double sigma0;

};

//==========================================================================

// A derived class for q g -> gamma*/Z0 q.

class Sigma2qg2gmZq : public Sigma2ffbargmZggm {

public:

  // Constructor.
  Sigma2qg2gmZq()
  {
    name = "q g-> gamma*/Z0 q";
    code = 242;
    fluxType = FluxType::QG;
    id3Mass = 23;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Values stored for later use.
  double sigma0;

};

//==========================================================================

// A derived class for f fbar' -> gamma*/Z0 gamma.

class Sigma2ffbar2gmZgm : public Sigma2ffbargmZggm {

public:

  // Constructor.
  Sigma2ffbar2gmZgm()
  {
    name = "f fbar -> gamma*/Z0 gamma";
    code = 243;
    fluxType = FluxType::FFBARSAME;
    id3Mass = 23;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Values stored for later use.
  double sigma0;

};

//==========================================================================

// A derived class for f gamma -> gamma*/Z0 f.

class Sigma2fgm2gmZf : public Sigma2ffbargmZggm {

public:

  // Constructor.
  Sigma2fgm2gmZf()
  {
    name = "f gamma -> gamma*/Z0 f";
    code = 244;
    fluxType = FluxType::FGAMMA;
    id3Mass = 23;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Values stored for later use.
  double sigma0;

};

//==========================================================================

// An intermediate class for f fbar -> W+- g/gamma and permutations.

class Sigma2ffbarWggm : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbarWggm() {}

  // Evaluate weight for gamma&/Z0 decay angle.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

};

//==========================================================================

// A derived class for q qbar' -> W+- g.

class Sigma2qqbar2Wg : public Sigma2ffbarWggm {

public:

  // Constructor.
  Sigma2qqbar2Wg()
  {
    name = "q qbar' -> W+- g";
    code = 251;
    fluxType = FluxType::FFBARCHG;
    id3Mass = 24;
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

  // Values stored for later use.
  double sigma0, openFracPos, openFracNeg;

};

//==========================================================================

// A derived class for q g -> W+- q'.

class Sigma2qg2Wq : public Sigma2ffbarWggm {

public:

  // Constructor.
  Sigma2qg2Wq()
  {
    name = "q g-> W+- q'";
    code = 252;
    fluxType = FluxType::QG;
    id3Mass = 24;
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

  // Values stored for later use.
  double sigma0, openFracPos, openFracNeg;

};

//==========================================================================

// A derived class for f fbar' -> W+- gamma.

class Sigma2ffbar2Wgm : public Sigma2ffbarWggm {

public:

  // Constructor.
  Sigma2ffbar2Wgm()
  {
    name = "f fbar' -> W+- gamma";
    code = 253;
    fluxType = FluxType::FFBARCHG;
    id3Mass = 24;
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

  // Values stored for later use.
  double sigma0, openFracPos, openFracNeg;

};

//==========================================================================

// A derived class for f gamma -> W+- f'.

class Sigma2fgm2Wf : public Sigma2ffbarWggm {

public:

  // Constructor.
  Sigma2fgm2Wf()
  {
    name = "f gamma -> W+- f'";
    code = 254;
    fluxType = FluxType::FGAMMA;
    id3Mass = 24;
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

  // Values stored for later use.
  double sigma0, openFracPos, openFracNeg;

};
//==========================================================================

// A derived class for gamma gamma -> f fbar.

class Sigma2gmgm2ffbar : public Sigma2Process {

public:

  // Constructor.
  Sigma2gmgm2ffbar(int idIn, int codeIn) : idNew(idIn)
  {
    code = codeIn;
    fluxType = FluxType::GAMMAGAMMA;
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

  // Member variables.
  int    idNew, idMass, idNow;
  double ef4, s34Avg, sigTU, sigma, openFracPair;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SigmaEW_H
