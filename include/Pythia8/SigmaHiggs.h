// SigmaHiggs.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// Part of code written by Marc Montull, CERN summer student 2007.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for Higgs process differential cross sections.
// Contains classes derived from SigmaProcess via Sigma2Process.

#ifndef Pythia8_SigmaHiggs_H
#define Pythia8_SigmaHiggs_H

#include "Pythia8/SigmaProcess.h"

namespace Pythia8 {

//==========================================================================

// A derived class for f fbar -> H0 (SM), H1, H2 or A3 (BSM).

class Sigma1ffbar2H : public SigmaProcess {

public:

  // Constructor.
  Sigma1ffbar2H(int higgsTypeIn) : SigmaProcess(ProcessType::P2to1), higgsType(higgsTypeIn)
  {
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

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // An H0, H1, H2 or A3 resonance object provides coupling
  // and propagator expressions.
  ParticleDataEntry* HResPtr;
  double mRes, GammaRes, m2Res, GamMRat, sigBW, widthOut;
  int    higgsType;
};

//==========================================================================

// A derived class for g g -> H0 (SM), H1, H2 or A3 (BSM).

class Sigma1gg2H : public SigmaProcess {

public:

  // Constructor.
  Sigma1gg2H(int higgsTypeIn) : SigmaProcess(ProcessType::P2to1), higgsType(higgsTypeIn)
  {
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

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // A H0, H1, H2 or A3 resonance object provides coupling
  // and propagator expressions.
  ParticleDataEntry* HResPtr;
  double mRes, GammaRes, m2Res, GamMRat, sigma;
  int    higgsType;
};

//==========================================================================

// A derived class for gamma gamma -> H0 (SM Higgs), H1, H2 or A3 (BSM Higgs).

class Sigma1gmgm2H : public SigmaProcess {

public:

  // Constructor.
  Sigma1gmgm2H(int higgsTypeIn) :SigmaProcess(ProcessType::P2to1), higgsType(higgsTypeIn)
  {
    fluxType = FluxType::GAMMAGAMMA;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // A H0, H1, H2 or A3 resonance object provides coupling
  // and propagator expressions.
  ParticleDataEntry* HResPtr;
  double mRes, GammaRes, m2Res, GamMRat, sigma;
  int    higgsType;
};

//==========================================================================

// A derived class for f fbar -> H Z0.
// (H can be H0 SM or H1, H2, A3 from BSM).
class Sigma2ffbar2HZ : public SigmaProcess {

public:

  // Constructor.
  Sigma2ffbar2HZ(int higgsTypeIn) : SigmaProcess(ProcessType::P2to2), higgsType(higgsTypeIn)
  {
    fluxType = FluxType::FFBARSAME;
    isSChannel = true;
    id3Mass = 23;
    id4Mass = 23;
    gmZmode = 2;
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
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Store Z0 mass and width.
  double mZ, widZ, mZS, mwZS, thetaWRat, sigma0, openFracPair, coup2Z;
  int    higgsType;
};

//==========================================================================

// A derived class for f fbar -> H W+- (Standard Model Higgs).
// (H can be H0 SM or H1, H2, A3 from BSM).

class Sigma2ffbar2HW : public SigmaProcess {

public:

  // Constructor.
  Sigma2ffbar2HW(int higgsTypeIn) : SigmaProcess(ProcessType::P2to2), higgsType(higgsTypeIn)
  {
    fluxType = FluxType::FFBARCHG;
    isSChannel = true;
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

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Store W+- mass and width, and couplings.
  double mW, widW, mWS, mwWS, thetaWRat, sigma0, openFracPairPos,
         openFracPairNeg, coup2W;
  int    higgsType, idRes;
};

//==========================================================================

// A derived class for f f' -> H f f' (Z0 Z0 fusion of SM or BSM Higgs).
// (H can be H0 SM or H1, H2, A3 from BSM).

class Sigma3ff2HfftZZ : public SigmaProcess {

public:

  // Constructor.
  Sigma3ff2HfftZZ(int higgsTypeIn) : SigmaProcess(ProcessType::P2to3), higgsType(higgsTypeIn)
  {
    fluxType = FluxType::FF;
    idTchan1 = 23;
    idTchan2 = 23;
    tChanFracPow1 = 0.05;
    tChanFracPow2 = 0.9;
    useMirrorWeight = true;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Store standard factors.
  double mZS, prefac, sigma1, sigma2, openFrac, coup2Z;
  int    higgsType, codeSave, idRes;
};

//==========================================================================

// A derived class for f_1 f_2 -> H f_3 f_4 (W+ W- fusion of SM or BSM Higgs).
// (H can be H0 SM or H1, H2, A3 from BSM).

class Sigma3ff2HfftWW : public SigmaProcess {

public:

  // Constructor.
  Sigma3ff2HfftWW(int higgsTypeIn) : SigmaProcess(ProcessType::P2to3), higgsType(higgsTypeIn)
  {
    fluxType = FluxType::FF;
    idTchan1 = 24;
    idTchan2 = 24;
    tChanFracPow1 = 0.05;
    tChanFracPow2 = 0.9;
    useMirrorWeight = true;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Store standard prefactor.
  double mWS, prefac, sigma0, openFrac, coup2W;
  int    higgsType, codeSave, idRes;
};

//==========================================================================

// A derived class for g g -> H Q Qbar (Q Qbar fusion of SM or BSM Higgs).
// (H can be H0 SM or H1, H2, A3 from BSM).

class Sigma3gg2HQQbar : public SigmaProcess {

public:

  // Constructor.
  Sigma3gg2HQQbar(int idIn, int higgsTypeIn) : SigmaProcess(ProcessType::P2to3), idNew(idIn),
    higgsType(higgsTypeIn)
    {
      fluxType = FluxType::GG;
      tChanFracPow1 = 0.4;
      tChanFracPow2 = 0.2;
      useMirrorWeight = false;
    }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Store flavour-specific process information and standard prefactor.
  double prefac, sigma, openFracTriplet, coup2Q;
  int    idNew, higgsType, codeSave, idRes;

};

//==========================================================================

// A derived class for q qbar -> H Q Qbar (Q Qbar fusion of SM or BSM Higgs).
// (H can be H0 SM or H1, H2, A3 from BSM).

class Sigma3qqbar2HQQbar : public SigmaProcess {

public:

  // Constructor.
  Sigma3qqbar2HQQbar(int idIn, int higgsTypeIn) : SigmaProcess(ProcessType::P2to3), idNew(idIn),
    higgsType(higgsTypeIn)
  {
    fluxType = FluxType::QQBARSAME;
    tChanFracPow1 = 0.4;
    tChanFracPow2 = 0.2;
    useMirrorWeight = false;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Store flavour-specific process information and standard prefactor.
  double prefac, sigma, openFracTriplet, coup2Q;
  int    idNew, higgsType, codeSave, idRes;

};

//==========================================================================

// A derived class for q g -> H q (SM or BSM Higgs).
// (H can be H0 SM or H1, H2, A3 from BSM).

class Sigma2qg2Hq : public SigmaProcess {

public:

  // Constructor.
  Sigma2qg2Hq(int idIn, int higgsTypeIn) : SigmaProcess(ProcessType::P2to2), idNew(idIn),
    higgsType(higgsTypeIn)
  {
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

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Store flavour-specific process information and standard prefactor.
  double m2W, thetaWRat, sigma, openFrac;
  int    idNew, higgsType, codeSave, idRes;

};

//==========================================================================

// A derived class for g g -> H0 g (SM or BSM Higgs via heavy top loop).
// (H can be H0 SM or H1, H2, A3 from BSM).

class Sigma2gg2Hglt : public SigmaProcess {

public:

  // Constructor.
  Sigma2gg2Hglt(int higgsTypeIn) : SigmaProcess(ProcessType::P2to2), higgsType(higgsTypeIn)
  {
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

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Store standard prefactor.
  double widHgg, sigma, openFrac;
  int    higgsType, codeSave, idRes;
};

//==========================================================================

// A derived class for q g -> H q (SM or BSM Higgs via heavy top loop).
// (H can be H0 SM or H1, H2, A3 from BSM).

class Sigma2qg2Hqlt : public SigmaProcess {

public:

  // Constructor.
  Sigma2qg2Hqlt(int higgsTypeIn) : SigmaProcess(ProcessType::P2to2), higgsType(higgsTypeIn)
  {
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

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Store standard prefactor.
  double widHgg, sigma, openFrac;
  int    higgsType, codeSave, idRes;
};

//==========================================================================

// A derived class for q qbar -> H g (SM or BSM Higgs via heavy top loop).
// (H can be H0 SM or H1, H2, A3 from BSM).

class Sigma2qqbar2Hglt : public SigmaProcess {

public:

  // Constructor.
  Sigma2qqbar2Hglt(int higgsTypeIn) : SigmaProcess(ProcessType::P2to2), higgsType(higgsTypeIn)
  {
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

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Store standard prefactor.
  double widHgg, sigma, openFrac;
  int    higgsType, codeSave, idRes;
};

//==========================================================================

// A derived class for f fbar' -> H+-.

class Sigma1ffbar2Hchg : public SigmaProcess {

public:

  // Constructor.
  Sigma1ffbar2Hchg() : SigmaProcess(ProcessType::P2to1)
  {
    name = "f fbar' -> H+-";
    code = 1061;
    fluxType = FluxType::FFBARCHG;
    resonanceA = 37;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // A H0 resonance object provides coupling and propagator expressions.
  ParticleDataEntry* HResPtr;
  double mRes, GammaRes, m2Res, GamMRat, m2W, thetaWRat, tan2Beta, sigBW,
         widthOutPos, widthOutNeg;

};

//==========================================================================

// A derived class for q g -> H+- q'.

class Sigma2qg2Hchgq : public SigmaProcess {

public:

  // Constructor.
  Sigma2qg2Hchgq(int idIn, int codeIn, stringref nameIn) : SigmaProcess(ProcessType::P2to2), idNew(idIn)
  {
    code = codeIn;
    name = nameIn;
    fluxType = FluxType::QG;
    id3Mass = 37;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Store flavour-specific process information and standard prefactor.
  int    idNew, codeSave, idOld, idUp, idDn;
  double m2W, thetaWRat, tan2Beta, sigma, openFracPos, openFracNeg;

};

//==========================================================================

// A derived class for f fbar -> A0(H_3) h0(H_1) or A0(H_3) H0(H_2).

class Sigma2ffbar2A3H12 : public SigmaProcess {

public:

  // Constructor.
  Sigma2ffbar2A3H12(int higgsTypeIn) : SigmaProcess(ProcessType::P2to2), higgsType(higgsTypeIn)
  {
    fluxType = FluxType::FFBARSAME;
    id3Mass = 36;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Store flavour-specific process information and standard prefactor.
  int    higgsType, higgs12, codeSave;
  double coupZA3H12, m2Z, mGammaZ, thetaWRat, openFrac, sigma0;

};

//==========================================================================

// A derived class for f fbar -> H+- h0(H_1) or H+- H0(H_2).

class Sigma2ffbar2HchgH12 : public SigmaProcess {

public:

  // Constructor.
  Sigma2ffbar2HchgH12(int higgsTypeIn) : SigmaProcess(ProcessType::P2to2), higgsType(higgsTypeIn)
  {
    fluxType = FluxType::FFBARCHG;
    id3Mass = 37;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Store flavour-specific process information and standard prefactor.
  int    higgsType, higgs12, codeSave;
  double coupWHchgH12, m2W, mGammaW, thetaWRat, openFracPos, openFracNeg,
         sigma0;

};

//==========================================================================

// A derived class for f fbar -> H+ H-.

class Sigma2ffbar2HposHneg : public SigmaProcess {

public:

  // Constructor.
  Sigma2ffbar2HposHneg() : SigmaProcess(ProcessType::P2to2)
  {
    name = "f fbar -> H+ H-";
    code = 1085;
    fluxType = FluxType::FFBARSAME;
    id3Mass = 37;
    id4Mass = 37;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

private:

  // Store flavour-specific process information and standard prefactor.
  double m2Z, mGammaZ, thetaWRat, eH, lH, openFrac, gamSig, intSig, resSig;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SigmaHiggs_H
