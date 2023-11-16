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

class Sigma0nonDiffractive : public Sigma0Process {

public:

  // Constructor.
  Sigma0nonDiffractive()
  {
    nameSave2 = "non-diffractive";
    codeSave = 101;
    isNonDiffSave = true;
  }

  // Evaluate sigma.
  virtual double sigmaHat() {return sigmaTotPtr->sigmaND();}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol() {}

  // Info on the subprocess.
  // virtual cstring name()      const {return "non-diffractive";}
  // virtual int    code()      const {return 101;}
  // virtual bool   isNonDiff() const {return true;}

private:

};

//==========================================================================

// A derived class for elastic scattering A B -> A B.

class Sigma0AB2AB : public Sigma0Process {

public:

  // Constructor.
  Sigma0AB2AB()
  {
    nameSave2 = "A B -> A B elastic";
    codeSave = 102;
    isResolvedSave = false;
  }

  // Evaluate sigma.
  virtual double sigmaHat() {return sigmaTotPtr->sigmaEl();}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  // virtual cstring name()       const {return "A B -> A B elastic";}
  // virtual int    code()       const {return 102;}
  // virtual bool   isResolved() const {return false;}

private:

};

//==========================================================================

// A derived class for single diffractive scattering A B -> X B.

class Sigma0AB2XB : public Sigma0Process {

public:

  // Constructor.
  Sigma0AB2XB()
  {
    nameSave2 = "A B -> X B single diffractive";
    codeSave = 103;
    isResolvedSave = false;
    isDiffASave = true;
  }

  // Evaluate sigma.
  virtual double sigmaHat() {return sigmaTotPtr->sigmaXB();}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  // virtual cstring name()       const {return "A B -> X B single diffractive";}
  // virtual int    code()       const {return 103;}
  // virtual bool   isResolved() const {return false;}
  // virtual bool   isDiffA()    const {return true;};

private:

};

//==========================================================================

// A derived class for single diffractive scattering A B -> A X.

class Sigma0AB2AX : public Sigma0Process {

public:

  // Constructor.
  Sigma0AB2AX()
  {
    nameSave2 = "A B -> A X single diffractive";
    codeSave = 104;
    isResolvedSave = false;
    isDiffBSave = true;
  }

  // Evaluate sigma.
  virtual double sigmaHat() {return sigmaTotPtr->sigmaAX();}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  // virtual cstring name()       const {return "A B -> A X single diffractive";}
  // virtual int    code()       const {return 104;}
  // virtual bool   isResolved() const {return false;}
  // virtual bool   isDiffB()    const {return true;};

// private:

};

//==========================================================================

// A derived class for double diffractive scattering A B -> X X.

class Sigma0AB2XX : public Sigma0Process {

public:

  // Constructor.
  Sigma0AB2XX()
  {
    nameSave2 = "A B -> X X double diffractive";
    codeSave = 105;
    isResolvedSave = false;
    isDiffASave = true;
    isDiffBSave = true;
  }

  // Evaluate sigma.
  virtual double sigmaHat() {return sigmaTotPtr->sigmaXX();}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  // virtual cstring name()       const {return "A B -> X X double diffractive";}
  // virtual int    code()       const {return 105;}
  // virtual bool   isResolved() const {return false;}
  // virtual bool   isDiffA()    const {return true;};
  // virtual bool   isDiffB()    const {return true;};

// private:

};

//==========================================================================

// A derived class for central diffractive scattering A B -> A X B.

class Sigma0AB2AXB : public Sigma0Process {

public:

  // Constructor.
  Sigma0AB2AXB()
  {
    nameSave2 = "A B -> A X B central diffractive";
    codeSave = 106;
    isResolvedSave = false;
    isDiffCSave = true;
    nFinalSave = 3;
  }

  // Evaluate sigma.
  virtual double sigmaHat() {return sigmaTotPtr->sigmaAXB();}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  // virtual cstring name()      const {return "A B -> A X B central diffractive";}
  // virtual int    code()       const {return 106;}
  // virtual int    nFinal()     const {return 3;}
  // virtual bool   isResolved() const {return false;}
  // virtual bool   isDiffC()    const {return true;};

// private:

};

//==========================================================================

// A derived class for g g -> g g.

class Sigma2gg2gg : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2gg()
  {
    nameSave2 = "g g -> g g";
    codeSave = 111;
    inFluxSave = "gg";
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  // virtual cstring name()   const {return "g g -> g g";}
  // virtual int    code()   const {return 111;}
  // virtual cstring inFlux() const {return "gg";}

private:

  // Values stored for colour flow selection.
  double sigTS, sigUS, sigTU, sigSum, sigma;

};

//==========================================================================

// A derived class for g g -> q qbar (q = u, d, s, i.e. almost massless).

class Sigma2gg2qqbar : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2qqbar()
  {
    nameSave2 = "g g -> q qbar (uds)";
    codeSave = 112;
    inFluxSave = "gg";
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  // virtual cstring name()   const {return "g g -> q qbar (uds)";}
  // virtual int    code()   const {return 112;}
  // virtual cstring inFlux() const {return "gg";}

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

class Sigma2qg2qg : public Sigma2Process {

public:

  // Constructor.
  Sigma2qg2qg()
  {
    nameSave2 = "q g -> q g";
    codeSave = 113;
    inFluxSave = "qg";
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  // virtual cstring name()   const {return "q g -> q g";}
  // virtual int    code()   const {return 113;}
  // virtual cstring inFlux() const {return "qg";}

private:

  // Values stored for colour flow selection.
  double sigTS, sigTU, sigSum, sigma;

};

//==========================================================================

// A derived class for q qbar' -> q qbar' or q q' -> q q'
// (qbar qbar' -> qbar qbar'), q' may be same as q.

class Sigma2qq2qq : public Sigma2Process {

public:

  // Constructor.
  Sigma2qq2qq()
  {
    nameSave2 = "q q(bar)' -> q q(bar)'";
    codeSave = 114;
    inFluxSave = "qq";
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  // virtual cstring name()   const {return "q q(bar)' -> q q(bar)'";}
  // virtual int    code()   const {return 114;}
  // virtual cstring inFlux() const {return "qq";}

 private:

  // Values stored for colour flow selection.
  double sigT, sigU, sigTU, sigST, sigSum;

};

//==========================================================================

// A derived class for q qbar -> g g.

class Sigma2qqbar2gg : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2gg()
  {
    nameSave2 = "q qbar -> g g";
    codeSave = 115;
    inFluxSave = "qqbarSame";
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  // virtual cstring name()   const {return "q qbar -> g g";}
  // virtual int    code()   const {return 115;}
  // virtual cstring inFlux() const {return "qqbarSame";}

 private:

  // Values stored for colour flow selection.
  double sigTS, sigUS, sigSum, sigma;

};

//==========================================================================

// A derived class for q qbar -> q' qbar'.

class Sigma2qqbar2qqbarNew : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2qqbarNew()
  {
    nameSave2 = "q qbar -> q' qbar' (uds)";
    codeSave = 116;
    inFluxSave = "qqbarSame";
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  // virtual cstring name()   const {return "q qbar -> q' qbar' (uds)";}
  // virtual int    code()   const {return 116;}
  // virtual cstring inFlux() const {return "qqbarSame";}

 private:

  // Number of quarks to be considered in massless approximation.
  int    nQuarkNew;

  // Values stored for colour flow selection.
  int    idNew;
  double mNew, m2New, sigS, sigma;

};

//==========================================================================

// A derived class for g g -> Q Qbar (Q = c, b or t).

class Sigma2gg2QQbar : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2QQbar(int idIn, int codeIn) : idNew(idIn)
  {
    inFluxSave = "gg";
    id3MassSave = idNew;
    id4MassSave = idNew;
    codeSave = codeIn;
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

  // Info on the subprocess.
  // virtual cstring name()    const {return nameSave;}
  // virtual int    code()    const {return codeSave;}
  // virtual cstring inFlux()  const {return "gg";}
  // virtual int    id3Mass() const {return idNew;}
  // virtual int    id4Mass() const {return idNew;}

 private:

  // Values stored for process type and colour flow selection.
  int    idNew;
  // stringbuf nameSave;
  double sigTS, sigUS, sigSum, sigma, openFracPair;

};

//==========================================================================

// A derived class for q qbar -> Q Qbar (Q = c, b or t).

class Sigma2qqbar2QQbar : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2QQbar(int idIn, int codeIn) : idNew(idIn)
  {
    inFluxSave = "qqbarSame";
    id3MassSave = idNew;
    id4MassSave = idNew;
    codeSave = codeIn;
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

  // Info on the subprocess.
  // virtual cstring name()    const {return nameSave;}
  // virtual int    code()    const {return codeSave;}
  // virtual cstring inFlux()  const {return "qqbarSame";}
  // virtual int    id3Mass() const {return idNew;}
  // virtual int    id4Mass() const {return idNew;}

 private:

  // Values stored for process type.
  int    idNew;
  // stringbuf nameSave;
  double sigma, openFracPair;

};

//==========================================================================

// A derived class for g g -> g g g.

class Sigma3gg2ggg : public Sigma3Process {

public:

  // Constructor.
  Sigma3gg2ggg()
  {
    nameSave2 = "g g -> g g g";
    codeSave = 131;
    nFinalSave = 3;
    inFluxSave = "gg";
    isQCD3bodySave = true;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  // virtual cstring name()       const {return "g g -> g g g";}
  // virtual int    code()       const {return 131;}
  // virtual int    nFinal()     const {return 3;}
  // virtual cstring inFlux()     const {return "gg";}
  // virtual bool   isQCD3body() const {return true;}

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

class Sigma3qqbar2ggg : public Sigma3Process {

public:

  // Constructor.
  Sigma3qqbar2ggg()
  {
    nameSave2 = "q qbar -> g g g";
    codeSave = 132;
    nFinalSave = 3;
    inFluxSave = "qqbarSame";
    isQCD3bodySave = true;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  // virtual cstring name()       const {return "q qbar -> g g g";}
  // virtual int    code()       const {return 132;}
  // virtual int    nFinal()     const {return 3;}
  // virtual cstring inFlux()     const {return "qqbarSame";}
  // virtual bool   isQCD3body() const {return true;}

protected:

  // Pick/map a random final state configuration
  int         config;
  inline void pickFinal() { config = int( 6 * rndmPtr->flat() ); }
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
    nameSave2 = "q g -> q g g";
    codeSave = 133;
    nFinalSave = 3;
    inFluxSave = "qg";
    isQCD3bodySave = true;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  // virtual cstring name()       const {return "q g -> q g g";}
  // virtual int    code()       const {return 133;}
  // virtual int    nFinal()     const {return 3;}
  // virtual cstring inFlux()     const {return "qg";}
  // virtual bool   isQCD3body() const {return true;}

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
    nameSave2 = "g g -> q qbar g";
    codeSave = 138;
    nFinalSave = 3;
    inFluxSave = "gg";
    isQCD3bodySave = true;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  // virtual cstring name()       const {return "g g -> q qbar g";}
  // virtual int    code()       const {return 138;}
  // virtual int    nFinal()     const {return 3;}
  // virtual cstring inFlux()     const {return "gg";}
  // virtual bool   isQCD3body() const {return true;}

private:

  // Number of quarks to be considered in massless approximation.
  int    nQuarkNew;

};

//==========================================================================

// A derived class for q q' -> q q' g

class Sigma3qq2qqgDiff : public Sigma3Process {

public:

  // Constructor.
  Sigma3qq2qqgDiff()
  {
    nameSave2 = "q(bar) q(bar)' -> q(bar) q(bar)' g";
    codeSave = 134;
    nFinalSave = 3;
    inFluxSave = "qq";
    isQCD3bodySave = true;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  // virtual cstring name()       const
    // {return "q(bar) q(bar)' -> q(bar) q(bar)' g";}
  // virtual int    code()       const {return 134;}
  // virtual int    nFinal()     const {return 3;}
  // virtual cstring inFlux()     const {return "qq";}
  // virtual bool   isQCD3body() const {return true;}

protected:

  // Pick/map a random final state configuration
  int         config;
  inline void pickFinal() { config = int( 6 * rndmPtr->flat() ); }
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
    nameSave2 = "q qbar -> q' qbar' g";
    codeSave = 136;
    nFinalSave = 3;
    inFluxSave = "qqbarSame";
    isQCD3bodySave = true;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  // virtual cstring name()       const {return "q qbar -> q' qbar' g";}
  // virtual int    code()       const {return 136;}
  // virtual int    nFinal()     const {return 3;}
  // virtual cstring inFlux()     const {return "qqbarSame";}
  // virtual bool   isQCD3body() const {return true;}

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
    nameSave2 = "q g -> q q' qbar'";
    codeSave = 139;
    nFinalSave = 3;
    inFluxSave = "qg";
    isQCD3bodySave = true;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  // virtual cstring name()       const {return "q g -> q q' qbar'";}
  // virtual int    code()       const {return 139;}
  // virtual int    nFinal()     const {return 3;}
  // virtual cstring inFlux()     const {return "qg";}
  // virtual bool   isQCD3body() const {return true;}

private:

  // Number of quarks to be considered in massless approximation.
  int    nQuarkNew;

  // gq and qg incoming
  double sigma[2];

};

//==========================================================================

// A derived class for q q -> q q g

class Sigma3qq2qqgSame : public Sigma3Process {

public:

  // Constructor.
  Sigma3qq2qqgSame()
  {
    nameSave2 = "q(bar) q(bar) -> q(bar) q(bar) g";
    codeSave = 135;
    nFinalSave = 3;
    inFluxSave = "qq";
    isQCD3bodySave = true;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  // virtual cstring name()       const
    // {return "q(bar) q(bar) -> q(bar) q(bar) g";}
  // virtual int    code()       const {return 135;}
  // virtual int    nFinal()     const {return 3;}
  // virtual cstring inFlux()     const {return "qq";}
  // virtual bool   isQCD3body() const {return true;}

protected:

  // Pick/map a random final state configuration
  int         config;
  inline void pickFinal() { config = int( 6 * rndmPtr->flat() ); }
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
    nameSave2 = "q qbar -> q qbar g";
    codeSave = 137;
    nFinalSave = 3;
    inFluxSave = "qqbarSame";
    isQCD3bodySave = true;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  // virtual cstring name()       const {return "q qbar -> q qbar g";}
  // virtual int    code()       const {return 137;}
  // virtual int    nFinal()     const {return 3;}
  // virtual cstring inFlux()     const {return "qqbarSame";}
  // virtual bool   isQCD3body() const {return true;}

// private:

};

//==========================================================================

// A derived class for q g -> q qbar q; same flavour.
// Derived from Sigma3qq2qqgSame

class Sigma3qg2qqqbarSame : public Sigma3qq2qqgSame {

public:

  // Constructor.
  Sigma3qg2qqqbarSame()
  {
    nameSave2 = "q g -> q q qbar";
    codeSave = 140;
    nFinalSave = 3;
    inFluxSave = "qg";
    isQCD3bodySave = true;
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  // virtual cstring name()       const {return "q g -> q q qbar";}
  // virtual int    code()       const {return 140;}
  // virtual int    nFinal()     const {return 3;}
  // virtual cstring inFlux()     const {return "qg";}
  // virtual bool   isQCD3body() const {return true;}

private:

  // gq and qg incoming
  double sigma[2];

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SigmaQCD_H
