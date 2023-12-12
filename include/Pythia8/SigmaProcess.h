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
#include "Pythia8/PythiaState.h"

namespace Pythia8 {

  class PythiaState;

// incoming parton pair type
enum class FluxType { 
  NONE, GG, QG, QQ, QQBAR, QQBARSAME, FF, FFBAR,  
  FFBARSAME, FFBARCHG, FGAMMA, GGAMMA, GAMMAGAMMA 
};

enum class ProcessType
{
  NONE,
  P2to0,
  P2to1,
  P2to2,
  P2to3,
  LHA,
};

// a parton in the beam and its PDF
struct InBeam { int id; double pdf; };

// a parton pair in the beam and their pdfSigma
struct InPair { InBeam A; InBeam B; double pdfSigma; };

// base class for all processes; derived classes represent a particular cross-section
class SigmaProcess
{

public: // public data

  // type of process
  ProcessType type = ProcessType::NONE; 

  // pointer to global Pythia state
  PythiaState* pState = nullptr;

  // getter functions (new)
  // Note: these were originally virtual functions
  string name; // name to identify process
  bool   convert2mb         = true;  // do we need to convert from GeV^-2 to mb
  bool   convertM2          = false; // For 2 -> 2 process optional conversion from |M|^2 to d(sigmaHat)/d(tHat).
  bool   isLHA              = false; // is this an LHA process?
  bool   isNonDiff          = false; // is it not a diffractive process?
  bool   isResolved         = true;  // is it a resolved process?
  bool   isDiffA            = false; // is product A diffractive?
  bool   isDiffB            = false; // is product B diffractive?
  bool   isDiffC            = false; // is product C diffractive?
  bool   isSUSY             = false; // is this a SUSY process?
  bool   allowNegativeSigma = false; // are negative cross-sections allowed?
  bool   isSChannel         = false; // 2 -> 2 and 2 -> 3 processes only through s-channel exchange.
  bool   isQCD3body         = false; // QCD 2 -> 3 processes need special phase space selection machinery.
  bool   useMirrorWeight    = false; // (for  2 -> 3 with two massive propagators)
  FluxType fluxType         = FluxType::NONE; // incomping parton type
  int    code               = 0; // code to identify proces
  int    nFinal             = 2; // number of final state particles
  int    gmZmode            = -1; // // Special process-specific gamma*/Z0 choice if >=0 (e.g. f fbar -> H0 Z0).
  int    id3Mass            = 0; // mass of first product for (2->3 processes)
  int    id4Mass            = 0; // mass of second product for (2->3 processes)
  int    id5Mass            = 0; // mass of third product for (2->3 processes)
  int    resonanceA         = 0; // id of first resonance
  int    resonanceB         = 0; // id of second resonance
  int    idSChannel         = 0; // NOAM: Insert an intermediate resonance in 2 -> 1 -> 2 (or 3) listings.
  int    idTchan1           = 0; // (for  2 -> 3 with two massive propagators)
  int    idTchan2           = 0; // (for  2 -> 3 with two massive propagators)
  double tChanFracPow1      = 0.3; // (for  2 -> 3 with two massive propagators)
  double tChanFracPow2      = 0.3; // (for  2 -> 3 with two massive propagators)
  
  // getter functions (original)
  double thetaMPI() const {return atan2( sinTheta, cosTheta);}
  // Note: rest were deleted since they were just accessers

public: // private data

  // Constants: could only be changed in the code itself.
  static const double CONVERT2MB, MASSMARGIN, COMPRELERR;
  static const int    NCOMPSTEP;

  // Initialization data, normally only set once.
  // WARNING: alias of some settings / particleData
  int    nQuarkIn, renormScale1, renormScale2, renormScale3, renormScale3VV,
         factorScale1, factorScale2, factorScale3, factorScale3VV;
  double Kfactor, mcME, mbME, mmuME, mtauME, renormMultFac, renormFixScale,
         factorMultFac, factorFixScale;

  // WARNING: alias of some settings
  // CP violation parameters for Higgs sector, normally only set once.
  int    higgsH1parity, higgsH2parity, higgsA3parity;
  double higgsH1eta, higgsH2eta, higgsA3eta, higgsH1phi, higgsH2phi,
         higgsA3phi;

  // Information on incoming beams.
  // WARNING: alias of beam info
  int    idA, idB;
  double mA, mB;
  bool   isLeptonA, isLeptonB, hasLeptonBeams;

  // Partons in beams, with PDF's.
  vector<InBeam> inBeamA;
  vector<InBeam> inBeamB;

  // Allowed colliding parton pairs, with pdf's.
  vector<InPair> inPair;


  // Store Q2 renormalization and factorization scales, and related values.
  double x1Save=0., x2Save=0.; // incoming parton momentum fractions
  
  // change to alpha -> alphaEMRen, alphaSRen
  double alpEM=0., alpS=0.; // electromagnetic and strong coupling constants
  double Q2RenSave=0., Q2FacSave=0., pdf1Save=0., pdf2Save=0., sigmaSumSave=0.;

  // Store flavour, colour, anticolour, mass, angles and the whole particle.
  int      id1, id2, id3, id4, id5; // i wonder if these are just duplicates of below
  int idSave[12]; // the ids of what?? i assume 1 + 2 -> 3 + 4 + 5 + ... but why would we need 12 of them?
  int colSave[12]; // color of first 12 particles
  int acolSave[12]; // anti-color of first 12 particles
  double   mSave[12] = {0,0,0,0, 0,0,0,0, 0,0,0,0}; // mass of first 12 particles
  double cosTheta, sinTheta, phi, sHMass, sHBeta, pT2Mass, pTFin;
  // phiMPI, sHBetaMPI, pT2MPI, pTMPIFin
  // TODO: is phi the same as phiMPI? -> yes
  // TODO: is pTFin the same as pTMIPFin? -> yes
  // TODO: is sHBeta the same as sHBetaMPI? -> yes
  Particle parton[12]; // i guess this just sotres the complete info. but why do we need the cached versions?

  // Minimal set of saved kinematics for trial interactions when
  // using the x-dependent matter profile of multiparton interactions.
  Particle partonT[12];
  double   mSaveT[12], pTFinT, cosThetaT, sinThetaT, phiT;

  // Calculate and store all modified masses and four-vectors
  // intended for matrix elements. Return false if failed.
  double   mME[12]; // mass for matrix elements; but why do we need a separate mass for them??
  Vec4     pME[12]; // 4-momentum for matrix elements

  // Store whether tHat and uHat are swapped (= same as swap 3 and 4).
  bool swapTU; // change to swappedTU

  // sH, tH, uH are Mandelstam variables; H for Hat ???
  // common subprocess kinematics quantities.
  double mH, sH, sH2;

  // (Sigma2Process) Store subprocess kinematics quantities.
  double tH, uH, tH2, uH2, m3, s3, m4, s4, pT2, runBW3, runBW4;

  // (Sigma3Process) Store subprocess kinematics quantities.
  double m5, s5, runBW5;
  Vec4   p3cm, p4cm, p5cm;
  // double m3, s3, m4, s4, runBW3, runBW4

public: // public functions

  // constructor
  SigmaProcess(ProcessType type);

  // destructor
  virtual ~SigmaProcess();

  // Perform simple initialization and store pointers.
  void init(PythiaState* pState);

  // Set up allowed flux of incoming partons. Default is no flux.
  bool initFlux();

  // Store kinematics and set scales for resolved 2 -> 1 process.
  void store1Kin(double x1in, double x2in, double sHin);

  // Store kinematics and set scales for resolved 2 -> 2 process.
  void store2Kin(double x1in, double x2in, double sHin, double tHin, double m3in, double m4in, double runBW3in, double runBW4in);

  // Store kinematics and set scales for resolved 2 -> 2 process.
  void store2KinMPI(double x1in, double x2in, double sHin, double tHin, double uHin, double alpSin, double alpEMin, bool needMasses, double m3in, double m4in);

  // Store kinematics and set scales for resolved 2 -> 3 process.
  void store3Kin(double x1in, double x2in, double sHin, Vec4ref p3cmIn, Vec4ref p4cmIn, Vec4ref p5cmIn, double m3in, double m4in, double m5in, double runBW3in, double runBW4in, double runBW5in);

  // set kinematics for a 2 -> 1 process
  void set1Kin(double x1in, double x2in, double sHin);

  // set kinematics for a 2 -> 2 process
  void set2Kin(double x1in, double x2in, double sHin, double tHin, double m3in, double m4in, double runBW3in, double runBW4in);

  // set kinematics for a 2 -> 2 process (Multiparton Interaction)
  void set2KinMPI(double x1in, double x2in, double sHin, double tHin,double uHin, double alpSin, double alpEMin, bool needMasses, double m3in, double m4in);
  
  // set kinematics for a 2 -> 3 process
  void set3Kin(double x1in, double x2in, double sHin, Vec4ref p3prel, Vec4ref p4prel, Vec4ref p5prel, double m3in, double m4in, double m5in, double runBW3in, double runBW4in, double runBW5in);

  // Wrapper to sigmaHat, to (a) store current incoming flavours and
  // (b) convert from GeV^-2 to mb where required.
  // For 2 -> 1/2 also (c) convert from from |M|^2 to d(sigmaHat)/d(tHat).
  double sigmaHatWrap(int id1in = 0, int id2in = 0);

  // Convolute above with parton flux and K factor. Sum over open channels.
  double sigmaPDF();

  // Select incoming parton channel and extract parton densities (resolved).
  void pickInState(int id1in = 0, int id2in = 0);

  // Perform kinematics for a Multiparton Interaction, in its rest frame.
  bool final2KinMPI( int = 0, int = 0, Vec4ref = 0., Vec4ref = 0., double = 0., double = 0.);

  // Set scale, when that is missing for an external LHA process.
  void setScale();

  // Save kinematics for trial interactions
  void saveKin();

  // Load kinematics for trial interactions
  void loadKin();

  // switch (active <--> saved) kinematics for trial interactions
  void swapKin();

public: // private functions

  // Calculate and store all modified masses and four-vectors
  // intended for matrix elements. Return false if failed.
  bool setupForMatrixElement();
  bool setupForMatrixElementInitial();

  // helper functions
  void setId(int id1in = 0, int id2in = 0, int id3in = 0, int id4in = 0, int id5in = 0);
  void setColAcol(int col1 = 0, int acol1 = 0, int col2 = 0, int acol2 = 0, int col3 = 0, int acol3 = 0, int col4 = 0, int acol4 = 0, int col5 = 0, int acol5 = 0);
  void swapColAcol();
  void swapCol1234();
  void swapCol12();
  void swapCol34();

  // TODO: why are these here rather than in the relavent process
  // Common code for top and Higgs secondary decay angular weights.
  double weightTopDecay(Event& process, int iResBeg, int iResEnd);
  double weightHiggsDecay(Event& process, int iResBeg, int iResEnd);

public: // virtual functions

  // Initialize process. Only used for some processes.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigma for unresolved, sigmaHat(sHat) for 2 -> 1 processes,
  // d(sigmaHat)/d(tHat) for (resolved) 2 -> 2 processes, and |M|^2 for
  // 2 -> 3 processes. Answer in "native" units, either mb or GeV^-2.
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for simultaneous flavours (only gamma*/Z0 gamma*/Z0).
  // WARNING: only overridden once
  virtual double weightDecayFlav(Event& process);

  // Evaluate weight for decay angular configuration.
  // iResBeg <= i < iResEnd is range of sister partons to test decays of.
  virtual double weightDecay(Event& process, int iResBeg, int iResEnd);

};

} // end namespace Pythia8

#endif // Pythia8_SigmaProcess_H
