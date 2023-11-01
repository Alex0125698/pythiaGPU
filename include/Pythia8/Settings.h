// Settings.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the settings database.
// Flag: helper class with bool flags.
// Mode: helper class with int modes.
// Parm: (short for parameter) helper class with double parameters.
// Word: helper class with string words.
// MVec: vector of Modes (integers).
// PVec: vector of Parms (doubles).
// Settings: maps of flags, modes, parms and words with input/output.

#ifndef Pythia8_Settings_H
#define Pythia8_Settings_H

#include "Pythia8/Info.h"
#include "Pythia8/PythiaStdlib.h"

namespace Pythia8 {

//==========================================================================

// Class for bool flags.

class Flag {

public:

  // Constructor
  Flag(string nameIn = " ", bool defaultIn = false) : name(nameIn),
    valNow(defaultIn) , valDefault(defaultIn) { }

  // Data members.
  string name;
  bool   valNow, valDefault;

};

//==========================================================================

// Class for integer modes.

class Mode {

public:

  // Constructor
  Mode(stringref nameIn = " ", int defaultIn = 0, bool hasMinIn = false,
    bool hasMaxIn = false, int minIn = 0,  int maxIn = 0,
    bool optOnlyIn = false) :  name(nameIn), valNow(defaultIn),
    valDefault(defaultIn), hasMin(hasMinIn), hasMax(hasMaxIn),
    valMin(minIn), valMax(maxIn), optOnly(optOnlyIn)  { }

  // Data members.
  string name;
  int    valNow, valDefault;
  bool   hasMin, hasMax;
  int    valMin, valMax;
  bool   optOnly;

};

//==========================================================================

// Class for double parms (where parm is shorthand for parameter).

class Parm {

public:

  // Constructor
  Parm(stringref nameIn = " ", double defaultIn = 0.,
    bool hasMinIn = false, bool hasMaxIn = false, double minIn = 0.,
    double maxIn = 0.) :  name(nameIn), valNow(defaultIn),
    valDefault(defaultIn), hasMin(hasMinIn), hasMax(hasMaxIn),
    valMin(minIn), valMax(maxIn) { }

  // Data members.
  string name;
  double valNow, valDefault;
  bool   hasMin, hasMax;
  double valMin, valMax;

};

//==========================================================================

// Class for string words.

class Word {

public:

  // Constructor
  Word(stringref nameIn = " ", stringref defaultIn = " ") : name(nameIn),
    valNow(defaultIn) , valDefault(defaultIn) { }

  // Data members.
  string name, valNow, valDefault;

};

//==========================================================================

// Class for vector of bool flags.

class FVec {

public:

  // Constructor
  FVec(stringref nameIn = " ", const vector<bool>& defaultIn = vector<bool>(1, false)) :
    name(nameIn), valNow(defaultIn) , valDefault(defaultIn) { }

  // Data members.
  string name;
  vector<bool> valNow, valDefault;

};

//==========================================================================

// Class for vector of integers.

class MVec {

public:

  // Constructor
  MVec(stringref nameIn = " ", const vector<int>& defaultIn = vector<int>(1, 0),
    bool hasMinIn = false, bool hasMaxIn = false, int minIn = 0,
    int maxIn = 0) :  name(nameIn), valNow(defaultIn),
    valDefault(defaultIn), hasMin(hasMinIn), hasMax(hasMaxIn),
    valMin(minIn), valMax(maxIn) { }

  // Data members.
  string name;
  vector<int> valNow, valDefault;
  bool   hasMin, hasMax;
  int    valMin, valMax;

};

//==========================================================================

// Class for vector of doubles.

class PVec {

public:

  // Constructor
  PVec(stringref nameIn = " ", const vector<double>& defaultIn = vector<double>(1, 0.),
    bool hasMinIn = false, bool hasMaxIn = false, double minIn = 0.,
    double maxIn = 0.) :  name(nameIn), valNow(defaultIn),
    valDefault(defaultIn), hasMin(hasMinIn), hasMax(hasMaxIn),
    valMin(minIn), valMax(maxIn) { }

  // Data members.
  string name;
  vector<double> valNow, valDefault;
  bool   hasMin, hasMax;
  double valMin, valMax;

};

//==========================================================================

// This class holds info on flags (bool), modes (int), parms (double),
// words (string), fvecs (vector of bool), mvecs (vector of int) and pvecs
// (vector of double).

class Settings {

public:

  // Constructor.
  Settings() : isInit(false), readingFailedSave(false) {}

  // Initialize Info pointer.
  void initPtr(Info* infoPtrIn) {infoPtr = infoPtrIn;}

  // Read in database from specific file.
  bool init(stringref startFile = "../xmldoc/Index.xml", bool append = false,
    ostream& os = cout) ;

  // Overwrite existing database by reading from specific file.
  bool reInit(stringref startFile = "../xmldoc/Index.xml", ostream& os = cout) ;

  // Read in one update from a single line.
  bool readString(stringref line, bool warn = true, ostream& os = cout) ;

  // Keep track whether any readings have failed, invalidating run setup.
  bool readingFailed() {return readingFailedSave;}

  // Write updates or everything to user-defined file.
  bool writeFile(stringref toFile, bool writeAll = false) ;
  bool writeFile(ostream& os = cout, bool writeAll = false) ;

  // Print out table of database, either all or only changed ones,
  // or ones containing a given string.
  void listAll(ostream& os = cout) {
    list( true, false, " ", os); }
  void listChanged(ostream& os = cout) {
    list (false, false, " ", os); }
  void list(stringref match, ostream& os = cout) {
    list (false, true, match, os); }

  // Give back current value(s) as a string, whatever the type.
  string output(stringref keyIn, bool fullLine = true);

  // Reset all values to their defaults.
  void resetAll() ;

  // Query existence of an entry.
  bool isFlag(stringref keyIn) {
    return (flags.find(toLower(keyIn)) != flags.end()); }
  bool isMode(stringref keyIn) {
    return (modes.find(toLower(keyIn)) != modes.end()); }
  bool isParm(stringref keyIn) {
    return (parms.find(toLower(keyIn)) != parms.end()); }
  bool isWord(stringref keyIn) {
    return (words.find(toLower(keyIn)) != words.end()); }
  bool isFVec(stringref keyIn) {
    return (fvecs.find(toLower(keyIn)) != fvecs.end()); }
  bool isMVec(stringref keyIn) {
    return (mvecs.find(toLower(keyIn)) != mvecs.end()); }
  bool isPVec(stringref keyIn) {
    return (pvecs.find(toLower(keyIn)) != pvecs.end()); }

  // Add new entry.
  void addFlag(stringref keyIn, bool defaultIn) {
    flags[toLower(keyIn)] = Flag(keyIn, defaultIn); }
  void addMode(stringref keyIn, int defaultIn, bool hasMinIn,
    bool hasMaxIn, int minIn, int maxIn, bool optOnlyIn = false) {
    modes[toLower(keyIn)] = Mode(keyIn, defaultIn, hasMinIn, hasMaxIn,
    minIn, maxIn, optOnlyIn); }
  void addParm(stringref keyIn, double defaultIn, bool hasMinIn,
    bool hasMaxIn, double minIn, double maxIn) { parms[toLower(keyIn)]
    = Parm(keyIn, defaultIn, hasMinIn, hasMaxIn, minIn, maxIn); }
  void addWord(stringref keyIn, stringref defaultIn) {
    words[toLower(keyIn)] = Word(keyIn, defaultIn); }
  void addFVec(stringref keyIn, const vector<bool>& defaultIn) {
    fvecs[toLower(keyIn)] = FVec(keyIn, defaultIn); }
  void addMVec(stringref keyIn, const vector<int>& defaultIn, bool hasMinIn,
    bool hasMaxIn, int minIn, int maxIn) { mvecs[toLower(keyIn)]
    = MVec(keyIn, defaultIn, hasMinIn, hasMaxIn, minIn, maxIn); }
   void addPVec(stringref keyIn, const vector<double>& defaultIn, bool hasMinIn,
    bool hasMaxIn, double minIn, double maxIn) { pvecs[toLower(keyIn)]
    = PVec(keyIn, defaultIn, hasMinIn, hasMaxIn, minIn, maxIn); }

  // Give back current value, with check that key exists.
  bool   flag(stringref keyIn);
  int    mode(stringref keyIn);
  double parm(stringref keyIn);
  string word(stringref keyIn);
  vector<bool>   fvec(stringref keyIn);
  vector<int>    mvec(stringref keyIn);
  vector<double> pvec(stringref keyIn);

  // Give back default value, with check that key exists.
  bool   flagDefault(stringref keyIn);
  int    modeDefault(stringref keyIn);
  double parmDefault(stringref keyIn);
  string wordDefault(stringref keyIn);
  vector<bool>   fvecDefault(stringref keyIn);
  vector<int>    mvecDefault(stringref keyIn);
  vector<double> pvecDefault(stringref keyIn);

  // Give back a map of all entries whose names match the string "match".
  map<string, Flag> getFlagMap(stringref match);
  map<string, Mode> getModeMap(stringref match);
  map<string, Parm> getParmMap(stringref match);
  map<string, Word> getWordMap(stringref match);
  map<string, FVec> getFVecMap(stringref match);
  map<string, MVec> getMVecMap(stringref match);
  map<string, PVec> getPVecMap(stringref match);

  // Change current value, respecting limits.
  void flag(stringref keyIn, bool nowIn);
  bool mode(stringref keyIn, int nowIn);
  void parm(stringref keyIn, double nowIn);
  void word(stringref keyIn, stringref nowIn);
  void fvec(stringref keyIn, const vector<bool>& nowIn);
  void mvec(stringref keyIn, const vector<int>& nowIn);
  void pvec(stringref keyIn, const vector<double>& nowIn);

  // Change current value, disregarding limits.
  void forceMode(stringref keyIn, int nowIn);
  void forceParm(stringref keyIn, double nowIn);
  void forceMVec(stringref keyIn, const vector<int>& nowIn);
  void forcePVec(stringref keyIn, const vector<double>& nowIn);

  // Restore current value to default.
  void resetFlag(stringref keyIn);
  void resetMode(stringref keyIn);
  void resetParm(stringref keyIn);
  void resetWord(stringref keyIn);
  void resetFVec(stringref keyIn);
  void resetMVec(stringref keyIn);
  void resetPVec(stringref keyIn);

private:

  // Pointer to various information on the generation.
  Info* infoPtr;

  // Map for bool flags.
  map<string, Flag> flags;

  // Map for integer modes.
  map<string, Mode> modes;

  // Map for double parms.
  map<string, Parm> parms;

  // Map for string words.
  map<string, Word> words;

  // Map for vectors of bool.
  map<string, FVec> fvecs;

  // Map for vectors of int.
  map<string, MVec> mvecs;

  // Map for vectors of double.
  map<string, PVec> pvecs;

  // Flags that initialization has been performed; whether any failures.
  bool isInit, readingFailedSave;

  // Print out table of database, called from listAll and listChanged.
  void list(bool doListAll, bool doListString, stringref match,
    ostream& os = cout);

  // Master switch for program printout.
  void printQuiet(bool quiet);

  // Restore settings used in tunes to e+e- and pp/ppbar data.
  void resetTuneEE();
  void resetTunePP();

  // Initialize tunes to e+e- and pp/ppbar data.
  void initTuneEE(int eeTune);
  void initTunePP(int ppTune);

  // Useful functions for string handling.
  string toLower(const string& name);
  bool   boolString(stringref tag);
  string attributeValue(stringref line, stringref attribute);
  bool   boolAttributeValue(stringref line, stringref attribute);
  int    intAttributeValue(stringref line, stringref attribute);
  double doubleAttributeValue(stringref line, stringref attribute);
  vector<bool>   boolVectorAttributeValue(stringref line, stringref attribute);
  vector<int>    intVectorAttributeValue(stringref line, stringref attribute);
  vector<double> doubleVectorAttributeValue(stringref line, stringref attribute);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_Settings_H
