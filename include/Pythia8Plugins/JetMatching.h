// JetMatching.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Richard Corke (implementation of MLM matching as
// in Alpgen for Alpgen input)
// and Stephen Mrenna (implementation of MLM-style matching as
// in Madgraph for Alpgen or Madgraph 5 input.)
// and Simon de Visscher, Stefan Prestel (implementation of shower-kT
// MLM-style matching and flavour treatment for Madgraph input)
// and Stefan Prestel (FxFx NLO jet matching with aMC@NLO.)
// This file provides the classes to perform MLM matching of
// Alpgen or MadGraph 5 input.
// Example usage is shown in main32.cc, and further details
// can be found in the 'Jet Matching Style' manual page.

#ifndef Pythia8_JetMatching_H
#define Pythia8_JetMatching_H

// Includes
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/GeneratorInput.h"

namespace Pythia8 {

//==========================================================================

// Declaration of main JetMatching class to perform MLM matching.
// Note that it is defined with virtual inheritance, so that it can
// be combined with other UserHooks classes, see e.g. main33.cc.

class JetMatching : virtual public UserHooks {

public:

  // Constructor and destructor
 JetMatching() : cellJet(NULL), slowJet(NULL), slowJetHard(NULL) {}
  ~JetMatching() {
    if (cellJet) delete cellJet;
    if (slowJet) delete slowJet;
    if (slowJetHard) delete slowJetHard;
  }

  // Initialisation
  virtual bool initAfterBeams() = 0;

  // Process level vetos
  bool canVetoProcessLevel() { return doMerge; }
  bool doVetoProcessLevel(Event& process) {
    eventProcessOrig = process;
    return false;
  }

  // Parton level vetos (before beam remnants and resonance decays)
  bool canVetoPartonLevelEarly() { return doMerge; }
  bool doVetoPartonLevelEarly(const Event& event);

  // Shower step vetoes (after the first emission, for Shower-kT scheme)
  int  numberVetoStep() {return 1;}
  bool canVetoStep() { return false; }
  bool doVetoStep(int,  int, int, const Event& ) { return false; }

protected:

  // Constants to be changed for debug printout or extra checks.
  static const bool MATCHINGDEBUG, MATCHINGCHECK;

  // Different steps of the matching algorithm.
  virtual void sortIncomingProcess(const Event &)=0;
  virtual void jetAlgorithmInput(const Event &, int)=0;
  virtual void runJetAlgorithm()=0;
  virtual bool matchPartonsToJets(int)=0;
  virtual int  matchPartonsToJetsLight()=0;
  virtual int  matchPartonsToJetsHeavy()=0;

  enum vetoStatus { NONE, LESS_JETS, MORE_JETS, HARD_JET, UNMATCHED_PARTON };
  enum partonTypes { ID_CHARM=4, ID_BOT=5, ID_TOP=6, ID_LEPMIN=11,
    ID_LEPMAX=16, ID_GLUON=21, ID_PHOTON=22 };

  // Master switch for merging
  bool   doMerge;
  // Switch for merging in the shower-kT scheme. Needed here because
  // the scheme uses different UserHooks functionality.
  bool   doShowerKt;

  // Maximum and current number of jets
  int    nJetMax, nJet;

  // Jet algorithm parameters
  int    jetAlgorithm;
  double eTjetMin, coneRadius, etaJetMax, etaJetMaxAlgo;

  // Internal jet algorithms
  CellJet* cellJet;
  SlowJet* slowJet;
  SlowJet* slowJetHard;

  // SlowJet specific
  int    slowJetPower;

  // Event records to store original incoming process, final-state of the
  // incoming process and what will be passed to the jet algorithm.
  // Not completely necessary to store all steps, but makes tracking the
  // steps of the algorithm a lot easier.
  Event eventProcessOrig, eventProcess, workEventJet;

  // Sort final-state of incoming process into light/heavy jets and 'other'
  vector<int> typeIdx[3];
  set<int>    typeSet[3];

  // Momenta output of jet algorithm (to provide same output regardless of
  // the selected jet algorithm)
  vector<Vec4> jetMomenta;

  // CellJet specific
  int    nEta, nPhi;
  double eTseed, eTthreshold;

  // Merging procedure parameters
  int    jetAllow, jetMatch, exclusiveMode;
  double coneMatchLight, coneRadiusHeavy, coneMatchHeavy;
  bool   exclusive;

  // Store the minimum eT/pT of matched light jets
  double eTpTlightMin;

};

//==========================================================================

// Declaration of main UserHooks class to perform Alpgen matching.

class JetMatchingAlpgen : virtual public JetMatching {

public:

  // Constructor and destructor
  JetMatchingAlpgen() { }
  ~JetMatchingAlpgen() { }

  // Initialisation
  bool initAfterBeams();

private:

  // Different steps of the matching algorithm.
  void sortIncomingProcess(const Event &);
  void jetAlgorithmInput(const Event &, int);
  void runJetAlgorithm();
  bool matchPartonsToJets(int);
  int  matchPartonsToJetsLight();
  int  matchPartonsToJetsHeavy();

  // Sorting utility
  void sortTypeIdx(vector < int > &vecIn);

  // Constants
  static const double GHOSTENERGY, ZEROTHRESHOLD;

};

//==========================================================================

// Declaration of main UserHooks class to perform Madgraph matching.

class JetMatchingMadgraph : virtual public JetMatching {

public:

  // Constructor and destructor
  JetMatchingMadgraph() { }
  ~JetMatchingMadgraph() { if (slowJetDJR) delete slowJetDJR; }

  // Initialisation
  bool initAfterBeams();

  // Process level vetos
  bool canVetoProcessLevel() { return doMerge; }
  bool doVetoProcessLevel(Event& process);

  // Shower step vetoes (after the first emission, for Shower-kT scheme)
  int  numberVetoStep() {return 1;}
  bool canVetoStep() { return doShowerKt; }
  bool doVetoStep(int,  int, int, const Event& );

  // Jet algorithm to access the jet separations in the cleaned event
  // after showering.
  SlowJet* slowJetDJR;
  // Functions to return the jet clustering scales and number of ME partons.
  // These are useful to investigate the matching systematics.
  vector<double> GetDJR() { return DJR;}
  pair<int,int> nMEpartons() { return nMEpartonsSave;}

protected:

  // Different steps of the matching algorithm.
  void sortIncomingProcess(const Event &);
  void jetAlgorithmInput(const Event &, int);
  void runJetAlgorithm();
  bool matchPartonsToJets(int);
  int  matchPartonsToJetsLight();
  int  matchPartonsToJetsHeavy();
  bool doShowerKtVeto(double pTfirst);

  // Functions to clear and set the jet clustering scales.
  void ClearDJR() { DJR.resize(0);}
  void SetDJR( const Event& event);
  // Functions to clear and set the jet clustering scales.
  void clear_nMEpartons() { nMEpartonsSave.first = nMEpartonsSave.second =-1;}
  void set_nMEpartons( const int nOrig, const int nMatch) {
    clear_nMEpartons();
    nMEpartonsSave.first  = nOrig;
    nMEpartonsSave.second = nMatch;
  };

  // Variables.
  vector<int> origTypeIdx[3];
  int    nQmatch;
  double qCut, qCutSq, clFact;
  bool   doFxFx;
  int    nPartonsNow;
  double qCutME, qCutMESq;

  // Vector to store the jet clustering scales.
  vector<double> DJR;
  // Pair of integers giving the number of ME partons read from LHEF and used
  // in the matching (can be different if some partons should not be matched)
  pair<int,int> nMEpartonsSave;

  // Function to get the current number of partons in the Born state, as
  // read from LHE.
  int npNLO();

};

} // end namespace Pythia8

#endif // end Pythia8_JetMatching_H
