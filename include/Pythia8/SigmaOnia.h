// SigmaOnia.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for charmonia/bottomonia process differential cross sections.
// Contains classes derived from SigmaProcess via Sigma2Process.

#ifndef Pythia8_SigmaOnia_H
#define Pythia8_SigmaOnia_H

#include "Pythia8/SigmaProcess.h"

namespace Pythia8 {

//==========================================================================

// A helper class used to setup the onia processes.

class SigmaOniaSetup {

public:

  // Constructors.
  SigmaOniaSetup() {};
  SigmaOniaSetup(Info* infoPtrIn, Settings* settingsPtrIn,
    ParticleData* particleDataPtrIn, int flavourIn);

  // Initialise the SigmaProcesses for gg, qg, or qqbar production.
  void setupSigma2gg(vector<SigmaProcess*> &procs, bool oniaIn = false);
  void setupSigma2qg(vector<SigmaProcess*> &procs, bool oniaIn = false);
  void setupSigma2qq(vector<SigmaProcess*> &procs, bool oniaIn = false);

private:

  // Intialise and check settings.
  void initStates(stringref wave, const vector<int> &states,
    vector<int> &jnums, bool &valid);
  void initSettings(stringref wave, unsigned int size,
    const vector<string> &names, vector< vector<double> > &pvecs, bool &valid);
  void initSettings(stringref wave, unsigned int size,
    const vector<string> &names, vector< vector<bool> > &fvecs, bool &valid);

  // Stored pointers.
  Info* infoPtr;
  Settings* settingsPtr;
  ParticleData* particleDataPtr;

  // Stored vectors of settings.
  vector<int> states3S1, states3PJ, states3DJ, spins3S1, spins3PJ, spins3DJ;
  vector<string> meNames3S1, meNames3PJ, meNames3DJ;
  vector< vector<double> > mes3S1, mes3PJ, mes3DJ;
  vector<string> ggNames3S1, qgNames3S1, qqNames3S1,
    ggNames3PJ, qgNames3PJ, qqNames3PJ, ggNames3DJ, qgNames3DJ, qqNames3DJ;
  vector< vector<bool> > ggs3S1, qgs3S1, qqs3S1, ggs3PJ, qgs3PJ, qqs3PJ,
    ggs3DJ, qgs3DJ, qqs3DJ;

  // Stored validity and production flags.
  bool onia, onia3S1, onia3PJ, onia3DJ, oniaFlavour;
  bool valid3S1, valid3PJ, valid3DJ;
  int flavour;
  string cat, key;

  // Stored parameters.
  double mSplit;

};

//==========================================================================

// A derived class for g g -> QQbar[3S1(1)] g (Q = c or b).

class Sigma2gg2QQbar3S11g : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2QQbar3S11g(int idHadIn, double oniumMEIn, int codeIn) :
    oniumME(oniumMEIn)
  {
    inFluxSave = "gg";
    codeSave = codeIn;
    id3MassSave = abs(idHadIn);
    idHad = idHadIn;
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

  // Values stored for process type and colour flow selection.
  int    idHad;
  double oniumME, sigma;

};

//==========================================================================

// A derived class for g g -> QQbar[3PJ(1)] g (Q = c or b, J = 0, 1 or 2).

class Sigma2gg2QQbar3PJ1g : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2QQbar3PJ1g(int idHadIn, double oniumMEIn, int jIn, int codeIn) :
    jSave(jIn), oniumME(oniumMEIn)
  {
    inFluxSave = "gg";
    codeSave = codeIn;
    id3MassSave = idHadIn;
    idHad = idHadIn;
  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

protected:

  ccstring namePrefix = "g g";
  ccstring namePostfix = "g";

  // Name pre-, post-, and mid-fix.
  ccstring nameMidfix() const {return (codeSave - codeSave%100)/100
      == 4 ? "ccbar" : "bbbar";}

  // Values stored for process type and colour flow selection.
  int    idHad, jSave;
  double oniumME, sigma;

};

//==========================================================================

// A derived class for q g -> QQbar[3PJ(1)] q (Q = c or b, J = 0, 1 or 2).

class Sigma2qg2QQbar3PJ1q : public Sigma2gg2QQbar3PJ1g {

public:

  // Constructor.
  Sigma2qg2QQbar3PJ1q(int idHadIn, double oniumMEIn, int jIn, int codeIn) :
    Sigma2gg2QQbar3PJ1g(idHadIn, oniumMEIn, jIn, codeIn)
  {
    inFluxSave = "qg";
    namePrefix = "q g";
    namePostfix = "q";
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

};

//==========================================================================

// A derived class for q qbar -> QQbar[3PJ(1)] g (Q = c or b, J = 0, 1 or 2).

class Sigma2qqbar2QQbar3PJ1g : public Sigma2gg2QQbar3PJ1g {

public:

  // Constructor.
  Sigma2qqbar2QQbar3PJ1g(int idHadIn, double oniumMEIn, int jIn, int codeIn) :
    Sigma2gg2QQbar3PJ1g(idHadIn, oniumMEIn, jIn, codeIn)
  {
    inFluxSave = "qqbarSame";
    namePrefix = "q qbar";
    namePostfix = "g";
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();
};

//==========================================================================

// A derived class for g g -> QQbar[3DJ(1)] g (Q = c or b).

class Sigma2gg2QQbar3DJ1g : public Sigma2gg2QQbar3PJ1g {

public:

  // Constructor.
  Sigma2gg2QQbar3DJ1g(int idHadIn, double oniumMEIn, int jIn, int codeIn) :
    Sigma2gg2QQbar3PJ1g(idHadIn, oniumMEIn, jIn, codeIn)
  {

  }

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

};

//==========================================================================

// A derived class for g g -> QQbar[X(8)] g (Q = c or b, X = 3S1, 1S0 or 3PJ).

class Sigma2gg2QQbarX8g : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2QQbarX8g(int idHadIn, double oniumMEIn, int stateIn,
    double mSplitIn, int codeIn) : stateSave(stateIn),
    oniumME(oniumMEIn), mSplit(mSplitIn)
  {
    inFluxSave = "gg";
    namePrefix = "g g";
    namePostfix = "g";
    id3MassSave = idHadIn;
    idHad = idHadIn;
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

protected:

  // Values stored for process type and colour flow selection.
  int    idHad, stateSave;
  double oniumME, sigma, mSplit;
  ccstring namePrefix = "g g";
  ccstring namePostfix = "g";
};

//==========================================================================

// A derived class for q g -> QQbar[X(8)] q (Q = c or b, X = 3S1, 1S0 or 3PJ).

class Sigma2qg2QQbarX8q : public Sigma2gg2QQbarX8g {

public:

  // Constructor.
  Sigma2qg2QQbarX8q(int idHadIn, double oniumMEIn, int stateIn,
    double mSplitIn, int codeIn) :
    Sigma2gg2QQbarX8g(idHadIn, oniumMEIn, stateIn, mSplitIn, codeIn)
  {
    inFluxSave = "qg";
    namePrefix = "q g";
    namePostfix = "q";
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

};

//==========================================================================

// A derived class for q qbar -> QQbar[X(8)] g (Q = c or b,
//   X = 3S1, 1S0 or 3PJ).

class Sigma2qqbar2QQbarX8g : public Sigma2gg2QQbarX8g {

public:

  // Constructor.
  Sigma2qqbar2QQbarX8g(int idHadIn, double oniumMEIn, int stateIn,
    double mSplitIn, int codeIn) :
    Sigma2gg2QQbarX8g(idHadIn, oniumMEIn, stateIn, mSplitIn, codeIn)
  {
    inFluxSave = "qqbarSame";
    namePrefix = "q qbar";
    namePostfix = "g";
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SigmaOnia_H
