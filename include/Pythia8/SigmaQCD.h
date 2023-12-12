// SigmaQCD.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for QCD process differential cross sections.
// Contains classes derived from SigmaProcess via Sigma(0/2)Process.

#ifndef Pythia8_SigmaQCD_H
#define Pythia8_SigmaQCD_H

#include "Pythia8/SigmaProcess.h"

namespace Pythia8 {

//==========================================================================

// A derived class for minimum-bias (inelastic, nondiffractive) events.

class Sigma0nonDiffractive : public SigmaProcess {

public:

  // Constructor.
  Sigma0nonDiffractive() : SigmaProcess(ProcessType::P2to0)
  {
    name = "non-diffractive";
    code = 101;
    isNonDiff = true;
  }

  // Evaluate sigma.
  virtual double sigmaHat() {return pState->sigmaTotal.sigmaND();}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol() {}

private:

};

//==========================================================================

// A derived class for elastic scattering A B -> A B.

class Sigma0AB2AB : public SigmaProcess {

public:

  // Constructor.
  Sigma0AB2AB() : SigmaProcess(ProcessType::P2to0)
  {
    name = "A B -> A B elastic";
    code = 102;
    isResolved = false;
  }

  // Evaluate sigma.
  virtual double sigmaHat() {return pState->sigmaTotal.sigmaEl();}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

};

//==========================================================================

// A derived class for single diffractive scattering A B -> X B.

class Sigma0AB2XB : public SigmaProcess {

public:

  // Constructor.
  Sigma0AB2XB() : SigmaProcess(ProcessType::P2to0)
  {
    name = "A B -> X B single diffractive";
    code = 103;
    isResolved = false;
    isDiffA = true;
  }

  // Evaluate sigma.
  virtual double sigmaHat() {return pState->sigmaTotal.sigmaXB();}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

};

//==========================================================================

// A derived class for single diffractive scattering A B -> A X.

class Sigma0AB2AX : public SigmaProcess {

public:

  // Constructor.
  Sigma0AB2AX() : SigmaProcess(ProcessType::P2to0)
  {
    name = "A B -> A X single diffractive";
    code = 104;
    isResolved = false;
    isDiffB = true;

  }

  // Evaluate sigma.
  virtual double sigmaHat() {return pState->sigmaTotal.sigmaAX();}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

};

//==========================================================================

// A derived class for double diffractive scattering A B -> X X.

class Sigma0AB2XX : public SigmaProcess {

public:

  // Constructor.
  Sigma0AB2XX() : SigmaProcess(ProcessType::P2to0)
  {
    name = "A B -> X X double diffractive";
    code = 105;
    isResolved = false;
    isDiffA = true;
    isDiffB = true;
  }

  // Evaluate sigma.
  virtual double sigmaHat() {return pState->sigmaTotal.sigmaXX();}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

};

//==========================================================================

// A derived class for central diffractive scattering A B -> A X B.

class Sigma0AB2AXB : public SigmaProcess {

public:

  // Constructor.
  Sigma0AB2AXB() : SigmaProcess(ProcessType::P2to0)
  {
    name = "A B -> A X B central diffractive";
    code = 106;
    nFinal = 3;
    isResolved = false;
    isDiffC = true;
  }

  // Evaluate sigma.
  virtual double sigmaHat() {return pState->sigmaTotal.sigmaAXB();}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

};

//==========================================================================

// A derived class for g g -> g g.

class Sigma2gg2gg : public SigmaProcess {

public:

  // Constructor.
  Sigma2gg2gg() : SigmaProcess(ProcessType::P2to2)
  {
    name = "g g -> g g";
    code = 111;
    fluxType = FluxType::GG;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Values stored for colour flow selection.
  double sigTS, sigUS, sigTU, sigSum, sigma;

};

//==========================================================================

// A derived class for g g -> q qbar (q = u, d, s, i.e. almost massless).

class Sigma2gg2qqbar : public SigmaProcess {

public:

  // Constructor.
  Sigma2gg2qqbar() : SigmaProcess(ProcessType::P2to2)
  {
    name = "g g -> q qbar (uds)";
    code = 112;
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

};

//==========================================================================

// A derived class for q g -> q g (q = u, d, s, c, b).
// Use massless approximation also for Q since no alternative.

class Sigma2qg2qg : public SigmaProcess {

public:

  // Constructor.
  Sigma2qg2qg() : SigmaProcess(ProcessType::P2to2)
  {
    name = "q g -> q g";
    code = 113;
    fluxType = FluxType::QG;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Values stored for colour flow selection.
  double sigTS, sigTU, sigSum, sigma;

};

//==========================================================================

// A derived class for q qbar' -> q qbar' or q q' -> q q'
// (qbar qbar' -> qbar qbar'), q' may be same as q.

class Sigma2qq2qq : public SigmaProcess {

public:

  // Constructor.
  Sigma2qq2qq() : SigmaProcess(ProcessType::P2to2)
  {
    name = "q q(bar)' -> q q(bar)'";
    code = 114;
    fluxType = FluxType::QQ;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

 private:

  // Values stored for colour flow selection.
  double sigT, sigU, sigTU, sigST, sigSum;

};

//==========================================================================

// A derived class for q qbar -> g g.

class Sigma2qqbar2gg : public SigmaProcess {

public:

  // Constructor.
  Sigma2qqbar2gg() : SigmaProcess(ProcessType::P2to2)
  {
    name = "q qbar -> g g";
    code = 115;
    fluxType = FluxType::QQBARSAME;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

 private:

  // Values stored for colour flow selection.
  double sigTS, sigUS, sigSum, sigma;

};

//==========================================================================

// A derived class for q qbar -> q' qbar'.

class Sigma2qqbar2qqbarNew : public SigmaProcess {

public:

  // Constructor.
  Sigma2qqbar2qqbarNew() : SigmaProcess(ProcessType::P2to2)
  {
    name = "q qbar -> q' qbar' (uds)";
    code = 116;
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

};

//==========================================================================

// A derived class for g g -> Q Qbar (Q = c, b or t).

class Sigma2gg2QQbar : public SigmaProcess {

public:

  // Constructor.
  Sigma2gg2QQbar(int idIn, int codeIn) : SigmaProcess(ProcessType::P2to2), idNew(idIn)
  {
    code = codeIn;
    fluxType = FluxType::GG;
    id3Mass = idNew;
    id4Mass = idNew;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for W decay angles in top decay (else inactive).
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

 private:

  // Values stored for process type and colour flow selection.
  int    idNew;
  double sigTS, sigUS, sigSum, sigma, openFracPair;

};

//==========================================================================

// A derived class for q qbar -> Q Qbar (Q = c, b or t).

class Sigma2qqbar2QQbar : public SigmaProcess {

public:

  // Constructor.
  Sigma2qqbar2QQbar(int idIn, int codeIn) : SigmaProcess(ProcessType::P2to2), idNew(idIn)
  {
    code = codeIn;
    fluxType = FluxType::QQBARSAME;
    id3Mass = idNew;
    id4Mass = idNew;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for W decay angles in top decay (else inactive).
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

 private:

  // Values stored for process type.
  int    idNew;
  double sigma, openFracPair;

};

//==========================================================================

// A derived class for g g -> g g g.

class Sigma3gg2ggg : public SigmaProcess {

public:

  // Constructor.
  Sigma3gg2ggg() : SigmaProcess(ProcessType::P2to3)
  {
    name = "g g -> g g g";
    code = 131;
    nFinal = 3;
    fluxType = FluxType::GG;
    isQCD3body = true;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Values stored for colour flow selection.
  double sigma;

  // Intermediate storage and calculation of four-products.
  double pp[6][6];
  double cycle(int i1, int i2, int i3, int i4, int i5) {return
    pp[i1][i2] * pp[i2][i3] * pp[i3][i4] * pp[i4][i5] * pp[i5][i1];}

};

//==========================================================================

// A derived class for q qbar -> g g g.

class Sigma3qqbar2ggg : public SigmaProcess {

public:

  // Constructor.
  Sigma3qqbar2ggg() :  SigmaProcess(ProcessType::P2to3)
  {
    name = "q qbar -> g g g";
    code = 132;
    nFinal = 3;
    fluxType = FluxType::QQBARSAME;
    isQCD3body = true;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

protected:

  // Pick/map a random final state configuration
  int         config;
  inline void pickFinal() { config = int( 6 * pState->rndm.flat() ); }
  inline void mapFinal();

  // |M|^2 calculation
  inline double m2Calc();

  // Four-vectors for |M|^2 calculation
  Vec4 pCM[5];

  // Intermediate storage and calculation of four-products
  double a[3], b[3], pp[3][3], ab[3][3];

  // Values stored for colour flow selection.
  double sigma;

};

//==========================================================================

// A derived class for q g -> q g g
// Derived from Sigma3qqbar2ggg

class Sigma3qg2qgg : public Sigma3qqbar2ggg {

public:

  // Constructor.
  Sigma3qg2qgg()
  {
    name = "q g -> q g g";
    code = 133;
    nFinal = 3;
    fluxType = FluxType::QG;
    isQCD3body = true;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Sigma for (qg) and (gq) incoming
  double sigma[2];

};

//==========================================================================

// A derived class for g g -> q qbar g
// Derived from Sigma3qqbar2ggg

class Sigma3gg2qqbarg : public Sigma3qqbar2ggg {

public:

  // Constructor.
  Sigma3gg2qqbarg()
  {
    name = "g g -> q qbar g";
    code = 138;
    nFinal = 3;
    fluxType = FluxType::GG;
    isQCD3body = true;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Number of quarks to be considered in massless approximation.
  int    nQuarkNew;

};

//==========================================================================

// A derived class for q q' -> q q' g

class Sigma3qq2qqgDiff : public SigmaProcess {

public:

  // Constructor.
  Sigma3qq2qqgDiff() :  SigmaProcess(ProcessType::P2to3)
  {
    name = "q(bar) q(bar)' -> q(bar) q(bar)' g";
    code = 134;
    nFinal = 3;
    fluxType = FluxType::QQ;
    isQCD3body = true;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

protected:

  // Pick/map a random final state configuration
  int         config;
  inline void pickFinal() { config = int( 6 * pState->rndm.flat() ); }
  inline void mapFinal();

  // |M|^2 calculation
  inline double m2Calc();

  // Kinematic configuration
  Vec4 pCM[5];

  // Four-products
  double s, t, u, sp, tp, up;

  // Cross section
  double sigma;

};

//==========================================================================

// A derived class for q qbar -> q' qbar' g
// Derived from Sigma3qq2qqgDiff

class Sigma3qqbar2qqbargDiff : public Sigma3qq2qqgDiff {

public:

  // Constructor.
  Sigma3qqbar2qqbargDiff()
  {
    name = "q qbar -> q' qbar' g";
    code = 136;
    nFinal = 3;
    fluxType = FluxType::QQBARSAME;
    isQCD3body = true;
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

};

//==========================================================================

// A derived class for q g -> q q' qbar'
// Derived from Sigma3qq2qqgDiff

class Sigma3qg2qqqbarDiff : public Sigma3qq2qqgDiff {

public:

  // Constructor.
  Sigma3qg2qqqbarDiff()
  {
    name = "q g -> q q' qbar'";
    code = 139;
    nFinal = 3;
    fluxType = FluxType::QG;
    isQCD3body = true;
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

  // Number of quarks to be considered in massless approximation.
  int    nQuarkNew;

  // gq and qg incoming
  double sigma[2];

};

//==========================================================================

// A derived class for q q -> q q g

class Sigma3qq2qqgSame : public SigmaProcess {

public:

  // Constructor.
  Sigma3qq2qqgSame() :  SigmaProcess(ProcessType::P2to3)
  {
    name = "q(bar) q(bar) -> q(bar) q(bar) g";
    code = 135;
    nFinal = 3;
    fluxType = FluxType::QQ;
    isQCD3body = true;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

protected:

  // Pick/map a random final state configuration
  int         config;
  inline void pickFinal() { config = int( 6 * pState->rndm.flat() ); }
  inline void mapFinal();

  // |M|^2 calculation
  inline double m2Calc();

  // Kinematic configuration
  Vec4 pCM[5];

  // Four-products
  double s, t, u, sp, tp, up;
  double ssp, ttp, uup, s_sp, t_tp, u_up;

  // Cross section
  double sigma;

};

//==========================================================================

// A derived class for q q -> q q g
// Derived from Sigma3qq2qqgSame

class Sigma3qqbar2qqbargSame : public Sigma3qq2qqgSame {

public:

  // Constructor.
  Sigma3qqbar2qqbargSame()
  {
    name = "q qbar -> q qbar g";
    code = 137;
    nFinal = 3;
    fluxType = FluxType::QQBARSAME;
    isQCD3body = true;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

};

//==========================================================================

// A derived class for q g -> q qbar q; same flavour.
// Derived from Sigma3qq2qqgSame

class Sigma3qg2qqqbarSame : public Sigma3qq2qqgSame {

public:

  // Constructor.
  Sigma3qg2qqqbarSame()
  {
    name = "q g -> q q qbar";
    code = 140;
    nFinal = 3;
    fluxType = FluxType::QG;
    isQCD3body = true;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // gq and qg incoming
  double sigma[2];

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SigmaQCD_H
