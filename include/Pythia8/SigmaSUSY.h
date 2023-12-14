// SigmaSUSY.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// Main authors of this file: N. Desai, P. Skands
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for Supersymmetric process differential cross sections.
// Contains classes derived from SigmaProcess via Sigma2Process.

#ifndef Pythia8_SigmaSUSY_H
#define Pythia8_SigmaSUSY_H

#include "Pythia8/PhaseSpace.h"
#include "Pythia8/PythiaComplex.h"
#include "Pythia8/SigmaProcess.h"
#include "Pythia8/SusyCouplings.h"

namespace Pythia8 {

//==========================================================================

// An intermediate class for SUSY 2 -> 2 with nontrivial decay angles.

class Sigma2SUSY : public SigmaProcess {

public:

  // Constructor.
  Sigma2SUSY() : SigmaProcess(ProcessType::P2to2) { };

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

 };

//==========================================================================

// A derived class for q qbar -> neutralino_i neutralino_j.

class Sigma2qqbar2chi0chi0 : public Sigma2SUSY {

public:

  // Constructor.
  Sigma2qqbar2chi0chi0() {};

  // Constructor.
  Sigma2qqbar2chi0chi0(int id3chiIn, int id4chiIn, int codeIn) {

    // Save ordering indices and process code
    id3chi   = id3chiIn;
    id4chi   = id4chiIn;
    code = codeIn;
    fluxType = FluxType::FF;
    resonanceA = 23;
    isSUSY = true;


    // Construct id codes from ordering indices.
    id3                  = 1000022;
    if (id3chi == 2) id3 = 1000023;
    if (id3chi == 3) id3 = 1000025;
    if (id3chi == 4) id3 = 1000035;
    if (id3chi == 5) id3 = 1000045;
    id4                  = 1000022;
    if (id4chi == 2) id4 = 1000023;
    if (id4chi == 3) id4 = 1000025;
    if (id4chi == 4) id4 = 1000035;
    if (id4chi == 5) id4 = 1000045;

    id3Mass = abs(id3);
    id4Mass = abs(id4);

  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  //  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual double getSigma0() const {return sigma0;}

 protected:

  // Basic process information
  int     id3chi, id4chi;

  // Values stored for later use
  double  sigma0, ui, uj, ti, tj, openFracPair;
  complex propZ;

  CoupSUSY* coupSUSYPtr;

};

//==========================================================================

// A derived class for q qbar -> neutralino_i chargino_j.

class Sigma2qqbar2charchi0 : public Sigma2qqbar2chi0chi0 {

public:

  // Constructor.
  Sigma2qqbar2charchi0(int id3chiIn, int id4chiIn, int codeIn) {

    // Save ordering indices and process code
    id3chi   = id3chiIn;
    id4chi   = id4chiIn;
    code = codeIn;
    resonanceA = 24;

    // Construct id codes from ordering indices.
    id3 = (abs(id3chi) == 2) ? 1000037 : 1000024;
    if (id3chi < 0)  id3 = -id3;

    id4                  = 1000022;
    if (id4chi == 2) id4 = 1000023;
    if (id4chi == 3) id4 = 1000025;
    if (id4chi == 4) id4 = 1000035;
    if (id4chi == 5) id4 = 1000045;

  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

protected :

  complex propW;

};

//==========================================================================

// A derived class for q qbar -> chargino+_i chargino-_j.

class Sigma2qqbar2charchar : public Sigma2qqbar2chi0chi0 {

public:

  // Constructor.
  Sigma2qqbar2charchar(int id3chiIn, int id4chiIn, int codeIn) {

    // Save ordering indices and process code
    id3chi   = id3chiIn;
    id4chi   = id4chiIn;
    code = codeIn;

    // Construct id codes from ordering indices.
    id3 = (abs(id3chi) == 2) ?  1000037 :  1000024;
    id4 = (abs(id4chi) == 2) ? -1000037 : -1000024;

  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

};

//==========================================================================

// A derived class for q g -> neutralino_i squark_j (and cc)

class Sigma2qg2chi0squark : public SigmaProcess {

public:

  // Constructor.
  Sigma2qg2chi0squark()  : SigmaProcess(ProcessType::P2to2) { };

  // Constructor.
  Sigma2qg2chi0squark(int id3chiIn, int id4sqIn, bool isUp, int codeIn)  : SigmaProcess(ProcessType::P2to2) {

    // Save ordering indices and process code
    id3chi   = id3chiIn;
    id4sq    = id4sqIn;
    code = codeIn;
    fluxType = FluxType::QG;
    isSUSY = true;

    // Construct id codes from ordering indices.
    id3                  = 1000022;
    if (id3chi == 2) id3 = 1000023;
    if (id3chi == 3) id3 = 1000025;
    if (id3chi == 4) id3 = 1000035;
    if (id3chi == 5) id3 = 1000045;
    id4                  = 1000001 + (isUp ? 1 : 0);
    if (id4sq  == 2) id4 = 1000003 + (isUp ? 1 : 0);
    if (id4sq  == 3) id4 = 1000005 + (isUp ? 1 : 0);
    if (id4sq  == 4) id4 = 2000001 + (isUp ? 1 : 0);
    if (id4sq  == 5) id4 = 2000003 + (isUp ? 1 : 0);
    if (id4sq  == 6) id4 = 2000005 + (isUp ? 1 : 0);

    id3Mass = abs(id3);
    id4Mass = abs(id4);

  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

 protected:

  // Basic process information
  int     id3chi, id4sq;

  // Values stored for later use
  double  sigma0, ui, uj, ti, tj, openFracPair;

  //SUSY couplings
  CoupSUSY* coupSUSYPtr;

};

//==========================================================================

// A derived class for q g -> chargino_i squark_j (incl cc)

class Sigma2qg2charsquark : public Sigma2qg2chi0squark {

public:

  // Constructor.
  Sigma2qg2charsquark(int id3chiIn, int id4sqIn, bool isUp, int codeIn) {

    // Save ordering indices and process code
    id3chi   = id3chiIn;
    id4sq    = id4sqIn;
    code = codeIn;
    fluxType = FluxType::QG;
    isSUSY = true;

    // Construct id codes from ordering indices.
    id3Sav                       = 1000024;
    if (abs(id3chi) == 2) id3Sav = 1000037;
    if (isUp)             id3Sav = -id3Sav;
    id4Sav                       = 1000001 + (isUp ? 1 : 0);
    if (id4sq  == 2) id4Sav      = 1000003 + (isUp ? 1 : 0);
    if (id4sq  == 3) id4Sav      = 1000005 + (isUp ? 1 : 0);
    if (id4sq  == 4) id4Sav      = 2000001 + (isUp ? 1 : 0);
    if (id4sq  == 5) id4Sav      = 2000003 + (isUp ? 1 : 0);
    if (id4sq  == 6) id4Sav      = 2000005 + (isUp ? 1 : 0);

    // Initial values, can be swapped to charge conjugates event by event.
    id3 = id3Sav;
    id4 = id4Sav;

    id3Mass = abs(id3);
    id4Mass = abs(id4);

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

  // Basic process information
  int id3Sav, id4Sav;

};

//==========================================================================

// A derived class for q q' -> ~q_i ~q_j

class Sigma2qq2squarksquark : public SigmaProcess {

public:

  // Constructor.
  // Sigma2qq2squarksquark() {}

  // Constructor.
  Sigma2qq2squarksquark(int id3In, int id4In, int codeIn) : SigmaProcess(ProcessType::P2to2) {

    // Save ordering indices and process code
    id3Sav = id3In;
    id4Sav = id4In;
    code = codeIn;
    // Initial values (flipped for c.c.)
    id3    = id3Sav;
    id4    = id4Sav;
    fluxType = FluxType::QQ;
    isSUSY = true;

    id3Mass = abs(id3Sav);
    id4Mass = abs(id4Sav);

  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() { return my_sigmaHat(); }
  double my_sigmaHat();
  virtual double sigmaPDF();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Basic process information
  int     id3Sav, id4Sav, iGen3, iGen4, nNeut;
  bool    isUD, onlyQCD;

  // Storage of mass squares
  double m2Glu;
  vector<double> m2Neut, m2Char;

  // Flavor-independent prefactors.
  double sigmaChar, sigmaNeut, sigmaGlu;
  double sigmaCharNeut, sigmaCharGlu, sigmaNeutGlu;
  double openFracPair;

  // Point-by-point info
  double tGlu, uGlu;
  vector<double> tNeut, uNeut, tChar, uChar;
  double sumCt, sumCu, sumNt, sumNu, sumGt, sumGu, sumInterference;

  //SUSY couplings
  CoupSUSY* coupSUSYPtr;
};

//==========================================================================

// A derived class for q qbar' -> ~q_i ~q*_j

class Sigma2qqbar2squarkantisquark : public SigmaProcess {

public:

  // Constructor.
  Sigma2qqbar2squarkantisquark() : SigmaProcess(ProcessType::P2to2) {}

  // Constructor.
  Sigma2qqbar2squarkantisquark(int id3In, int id4In, int codeIn) : SigmaProcess(ProcessType::P2to2) {

    // Save ordering indices and process code
    // (always store squark first, antisquark second)
    id3Sav = abs(id3In);
    id4Sav = -abs(id4In);
    code = codeIn;
    // Initial values
    id3    = id3Sav;
    id4    = id4Sav;
    fluxType = FluxType::QQ;
    isSUSY = true;

    id3Mass = abs(id3Sav);
    id4Mass = abs(id4Sav);
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

  // Basic process information
  int     id3Sav, id4Sav, iGen3, iGen4, nNeut;
  bool    isUD, isCC, onlyQCD;

  // Storage of mass squares
  double m2Glu;
  vector<double> m2Neut;

  // Flavor-independent prefactors: EW, strong, and interference
  double xW;
  double openFracPair;
  double sigmaEW, sigmaGlu, sigmaEWG;

  // Point-by-point info
  double tGlu, uGlu;
  vector<double> tNeut, uNeut;
  complex propZW;
  double sumColS, sumColT, sumInterference;

  //SUSY couplings
  CoupSUSY* coupSUSYPtr;

};

//==========================================================================

// A derived class for g g -> ~q ~q*

class Sigma2gg2squarkantisquark : public SigmaProcess {

public:

  // Constructor.
  // Sigma2gg2squarkantisquark() {
  // }

  // Constructor.
  Sigma2gg2squarkantisquark(int id34In, int codeIn) : SigmaProcess(ProcessType::P2to2) {

    // Save ordering indices and process code
    // (always store squark first, antisquark second)
    id3Sav = abs(id34In);
    id4Sav = -abs(id34In);
    code = codeIn;
    // Initial values
    id3    = id3Sav;
    id4    = id4Sav;
    fluxType = FluxType::GG;
    isSUSY = true;
    id3Mass = abs(id3Sav);
    id4Mass = abs(id4Sav);

  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() { Benchmark_start(sigmaHat_SUSY_Sigma2gg2squarkantisquark); return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Basic process information
  int     id3Sav, id4Sav;
  double sigma, m2Sq, openFracPair;

  //SUSY couplings
  CoupSUSY* coupSUSYPtr;

};

//==========================================================================

// A derived class for q g -> ~q ~g

class Sigma2qg2squarkgluino : public SigmaProcess {

public:

  // Constructor.
  // Sigma2qg2squarkgluino() {}

  // Constructor.
  Sigma2qg2squarkgluino(int id3In, int codeIn) : SigmaProcess(ProcessType::P2to2) {

    // Save ordering indices and process code
    id3Sav = abs(id3In);
    code = codeIn;
    // Initial values
    id3    = id3Sav;
    id4    = 1000021;
    fluxType = FluxType::QG;
    isSUSY = true;
    id4Mass= 1000021;
    id3Mass = abs(id3Sav);
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

  // Basic process information
  int     id3Sav;
  double sigmaA, sigmaB, comFacHat, m2Glu, m2Sq, openFracPair;

  //SUSY couplings
  CoupSUSY* coupSUSYPtr;

};

//==========================================================================

// A derived class for g g -> gluino gluino.

class Sigma2gg2gluinogluino : public SigmaProcess {

public:

  // Constructor.
  Sigma2gg2gluinogluino() : SigmaProcess(ProcessType::P2to2)
  {
    name = "g g -> gluino gluino";
    code = 1201;
    fluxType = FluxType::GG;
    id3Mass = 1000021;
    id4Mass = 1000021;
    isSUSY = true;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {Benchmark_start(sigmaHat_SUSY_Sigma2gg2gluinogluino); return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

private:

  // Values stored for process type and colour flow selection.
  double sigTS, sigUS, sigTU, sigSum, sigma, openFracPair;

  //SUSY couplings
  CoupSUSY* coupSUSYPtr;

};

//==========================================================================

// A derived class for q qbar -> gluino gluino.

class Sigma2qqbar2gluinogluino : public SigmaProcess {

public:

  // Constructor.
  Sigma2qqbar2gluinogluino() : SigmaProcess(ProcessType::P2to2)
  {
    name = "q qbar -> gluino gluino";
    code = 1202;
    fluxType = FluxType::QQ;
    id3Mass = 1000021;
    id4Mass = 1000021;
    isSUSY = true;
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

  // Values stored for process type and colour flow selection.
  double openFracPair, s34Avg, sigS, tHG, uHG, tHG2, uHG2;

  //SUSY couplings
  CoupSUSY* coupSUSYPtr;

};

//==========================================================================

class Sigma1qq2antisquark : public SigmaProcess {
public:

  // Constructor.
  // Sigma1qq2antisquark() {}


  Sigma1qq2antisquark(int id3In) : SigmaProcess(ProcessType::P2to1) 
  {

    idRes = id3In;
    fluxType = FluxType::QQ;
    isSUSY = true;
    resonanceA = idRes;

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
  virtual bool   isRPV()   const {return true;}

private:

  // Values stored for process type and colour flow selection.
  double mRes, GammaRes, m2Res, sigBW, widthOut;
  int    idRes;

  //SUSY couplings
  CoupSUSY* coupSUSYPtr;

};


//==========================================================================

// A derived class for q qbar -> neutralino_i gluino.

class Sigma2qqbar2chi0gluino : public Sigma2SUSY {

public:

  // Constructor.
  Sigma2qqbar2chi0gluino() {};

  // Constructor.
  Sigma2qqbar2chi0gluino(int id4chiIn, int codeIn) {

    // Save ordering indices and process code
    id3   = 1000021;
    id4chi   = id4chiIn;
    code = codeIn;


    // Construct id codes from ordering indices.
    id4                  = 1000022;
    if (id4chi == 2) id4 = 1000023;
    if (id4chi == 3) id4 = 1000025;
    if (id4chi == 4) id4 = 1000035;
    if (id4chi == 5) id4 = 1000045;

    fluxType = FluxType::FF;
    resonanceA = 23;
    isSUSY = true;

    id3Mass = abs(id3);
    id4Mass = abs(id4);
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  //  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Info on the subprocess.
  virtual double getSigma0() const {return sigma0;}

 protected:

  // Basic process information
  int     id3chi, id4chi;

  // Values stored for later use
  double  sigma0, ui, uj, ti, tj, openFracPair;

  CoupSUSY* coupSUSYPtr;

};

//==========================================================================

// A derived class for q qbar -> neutralino_i chargino_j.

class Sigma2qqbar2chargluino : public Sigma2qqbar2chi0gluino {

public:

  // Constructor.
  Sigma2qqbar2chargluino(int id4chiIn, int codeIn) {

    // Save ordering indices and process code
    id3   = 1000021;
    id4chi   = id4chiIn;
    code = codeIn;

    // Construct id codes from ordering indices.
    id4 = (abs(id4chi) == 2) ? 1000037 : 1000024;
    if (id4chi < 0)  id4 = -id4;

    fluxType = FluxType::FF;
    resonanceA = 24;
    isSUSY = true;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

protected :

  complex propW;

};

//==========================================================================

// A derived class for q qbar' -> ~q_i ~q*_j

class Sigma2qqbar2sleptonantislepton : public Sigma2qqbar2squarkantisquark {

public:

  // Constructor.
  // Sigma2qqbar2sleptonantislepton() {}

  // Constructor.
  Sigma2qqbar2sleptonantislepton(int id3In, int id4In, int codeIn) {

    // Save ordering indices and process code
    // (always store squark first, antisquark second)
    id3Sav = abs(id3In);
    id4Sav = -abs(id4In);
    code = codeIn;
    // Initial values
    id3    = id3Sav;
    id4    = id4Sav;

    fluxType = FluxType::QQ;
    isSUSY = true;

    id3Mass = abs(id3Sav);
    id4Mass = abs(id4Sav);
    
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

  // Basic process information
  int     id3Sav, id4Sav, iGen3, iGen4, nNeut;
  bool    isUD;

  // Storage of mass squares
  vector<double> m2Neut;

  // Flavor-independent prefactors: EW, strong, and interference
  double xW;
  double openFracPair;
  double sigmaEW;

  // Point-by-point info
  vector<double> tNeut, uNeut;
  complex propZW;
  double sumColS, sumColT, sumInterference;

  //SUSY couplings
  CoupSUSY* coupSUSYPtr;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SigmaSUSY_H
