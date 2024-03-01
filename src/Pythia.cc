// Pythia.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the Pythia class.

#include "Pythia8/Pythia.h"

// Access time information.
#include <ctime>

// Allow string and character manipulation.
#include <cctype>

namespace Pythia8 {

//==========================================================================

// The Pythia class.

//--------------------------------------------------------------------------

// The current Pythia (sub)version number, to agree with XML version.
const double Pythia::VERSIONNUMBERHEAD = PYTHIA_VERSION;
const double Pythia::VERSIONNUMBERCODE = 8.212;

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Maximum number of tries to produce parton level from given input.
const int Pythia::NTRY          = 10;

// Negative integer to denote that no subrun has been set.
const int Pythia::SUBRUNDEFAULT = -999;

//--------------------------------------------------------------------------

// Constructor.

// done
Pythia::Pythia(stringref xmlDir_, bool printBanner) {

  // !@!@
  thread_local string xmlDir; xmlDir = xmlDir_;

  // Benchmark::init();
  Benchmark_start(PythiaP);
  Benchmark_start(PythiaP_clearData); // trivial

  // Initial values for pointers to PDF's.
  useNewPdfA      = false;
  useNewPdfB      = false;
  useNewPdfHard   = false;
  useNewPdfPomA   = false;
  useNewPdfPomB   = false;
  pdfAPtr         = 0;
  pdfBPtr         = 0;
  pdfHardAPtr     = 0;
  pdfHardBPtr     = 0;
  pdfPomAPtr      = 0;
  pdfPomBPtr      = 0;

  // Initial values for pointers to Les Houches Event objects.
  doLHA           = false;
  useNewLHA       = false;
  pState.lhaUp        = 0;

  //Initial value for couplings pointer
  pState.couplings    = &couplings; // !!!!

  // Initial value for pointer to external decay handler.
  decayHandlePtr  = 0;

  // Initial value for pointer to user hooks.
  userHooksPtr    = 0;

  // Initial value for pointer to merging hooks.
  doMerging          = false;
  hasMergingHooks    = false;
  hasOwnMergingHooks = false;
  mergingHooksPtr    = 0;

  // Initial value for pointer to beam shape.
  useNewBeamShape = false;
  beamShapePtr    = 0;

  // Initial values for pointers to timelike and spacelike showers.
  useNewTimesDec  = false;
  useNewTimes     = false;
  useNewSpace     = false;
  timesDecPtr     = 0;
  timesPtr        = 0;
  spacePtr        = 0;

  Benchmark_stop(PythiaP_clearData);
  Benchmark_start(PythiaP_loadXmlSettings); // trivial

  // Find path to data files, i.e. xmldoc directory location.
  // Environment variable takes precedence, then constructor input,
  // and finally the pre-processor constant XMLDIR.
  xmlPath = "";
  const char* PYTHIA8DATA = "PYTHIA8DATA";
  char* envPath = getenv(PYTHIA8DATA);
  if (envPath != 0 && *envPath != '\0') {
    int i = 0;
    while (*(envPath+i) != '\0') xmlPath += *(envPath+(i++));
  } else {
    if (xmlDir[ xmlDir.length() - 1 ] != '/') xmlDir += "/";
    xmlPath = xmlDir;
    ifstream xmlFile((xmlPath + "Index.xml").c_str());
    if (!xmlFile.good()) xmlPath = XMLDIR;
    xmlFile.close();
  }
  if (xmlPath[ xmlPath.length() - 1 ] != '/') xmlPath += "/";

  // Read in files with all flags, modes, parms and words.
  // pState.settings.initPtr( &pState.info);
  string initFile = xmlPath + "Index.xml";
  isConstructed = pState.settings.init( initFile); // !!!
  if (!isConstructed) {
    pState.info.errorMsg("Abort from Pythia::Pythia: settings unavailable");
    return;
  }

  Benchmark_stop(PythiaP_loadXmlSettings);
  Benchmark_start(PythiaP_checkXmlSettings); // trivial

  // Check that XML version number matches code version number.
  double versionNumberXML = parm("Pythia:versionNumber");
  isConstructed = (abs(versionNumberXML - VERSIONNUMBERCODE) < 0.0005);
  if (!isConstructed) {
    ostringstream errCode;
    errCode << fixed << setprecision(3) << ": in code " << VERSIONNUMBERCODE
            << " but in XML " << versionNumberXML;
    pState.info.errorMsg("Abort from Pythia::Pythia: unmatched version numbers",
      errCode.str());
    return;
  }

  // Check that header version number matches code version number.
  isConstructed = (abs(VERSIONNUMBERHEAD - VERSIONNUMBERCODE) < 0.0005);
  if (!isConstructed) {
    ostringstream errCode;
    errCode << fixed << setprecision(3) << ": in code " << VERSIONNUMBERCODE
            << " but in header " << VERSIONNUMBERHEAD;
    pState.info.errorMsg("Abort from Pythia::Pythia: unmatched version numbers",
      errCode.str());
    return;
  }

  Benchmark_stop(PythiaP_checkXmlSettings);
  Benchmark_start(PythiaP_loadXmlParticleData); // trivial

  // Read in files with all particle data.
  pState.particleData.initPtr( &pState.info, &pState.settings, &pState.rndm, pState.couplings);
  string dataFile = xmlPath + "ParticleData.xml";
  isConstructed = pState.particleData.init( dataFile);
  if (!isConstructed) {
    pState.info.errorMsg("Abort from Pythia::Pythia: particle data unavailable");
    return;
  }

  Benchmark_stop(PythiaP_loadXmlParticleData);
  Benchmark_start(PythiaP_printBanner); // trivial

  // Write the Pythia banner to output.
  if (printBanner) banner();

  // Not initialized until at the end of the init() call.
  isInit = false;
  pState.info.addCounter(0);

}

//--------------------------------------------------------------------------

// Destructor.

// done
Pythia::~Pythia() {

  // Delete the PDF's created with new.
  if (useNewPdfHard && pdfHardAPtr != pdfAPtr) delete pdfHardAPtr;
  if (useNewPdfHard && pdfHardBPtr != pdfBPtr) delete pdfHardBPtr;
  if (useNewPdfA) delete pdfAPtr;
  if (useNewPdfB) delete pdfBPtr;
  if (useNewPdfPomA) delete pdfPomAPtr;
  if (useNewPdfPomB) delete pdfPomBPtr;

  // Delete the Les Houches object created with new.
  if (useNewLHA) delete pState.lhaUp;

  // Delete the MergingHooks object created with new.
  if (hasOwnMergingHooks) delete mergingHooksPtr;

  // Delete the BeamShape object created with new.
  if (useNewBeamShape) delete beamShapePtr;

  // Delete the timelike and spacelike showers created with new.
  if (useNewTimesDec) delete timesDecPtr;
  if (useNewTimes && !useNewTimesDec) delete timesPtr;
  if (useNewSpace) delete spacePtr;

}

//--------------------------------------------------------------------------

// Read in one update for a setting or particle data from a single line.

// done
bool Pythia::readString(stringref line, bool warn) {

  // Check that constructor worked.
  if (!isConstructed) return false;

  // If empty line then done.
  if (line.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos) return true;

  // If first character is not a letter/digit, then taken to be a comment.
  int firstChar = line.find_first_not_of(" \n\t\v\b\r\f\a");
  if (!isalnum(line[firstChar])) return true;

  // Send on particle data to the ParticleData database.
  if (isdigit(line[firstChar])) {
    bool passed = pState.particleData.readString(line, warn);
    if (passed) particleDataBuffer << line << endl;
    return passed;
  }

  // Everything else sent on to pState.Settings.
  return pState.settings.readString(line, warn);

}

//--------------------------------------------------------------------------

// Read in updates for settings or particle data from user-defined file.

// done
bool Pythia::readFile(stringref fileName, bool warn, int subrun) {

  // Check that constructor worked.
  if (!isConstructed) return false;

  // Open file for reading.
  const char* cstring = fileName.c_str();
  ifstream is(cstring);
  if (!is.good()) {
    pState.info.errorMsg("Error in Pythia::readFile: did not find file", fileName);
    return false;
  }

  // Hand over real work to next method.
  return readFile( is, warn, subrun);

}

//--------------------------------------------------------------------------

// Read in updates for settings or particle data
// from user-defined stream (or file).

// done
bool Pythia::readFile(istream& is, bool warn, int subrun) {

  // Check that constructor worked.
  if (!isConstructed) return false;

  // Read in one line at a time.
  string line;
  bool isCommented = false;
  bool accepted = true;
  int subrunNow = SUBRUNDEFAULT;
  while ( getline(is, line) ) {

    // Check whether entering, leaving or inside commented-commands section.
    int commentLine = readCommented( line);
    if      (commentLine == +1)  isCommented = true;
    else if (commentLine == -1)  isCommented = false;
    else if (isCommented) ;

    else {
      // Check whether entered new subrun.
      int subrunLine = readSubrun( line, warn);
      if (subrunLine >= 0) subrunNow = subrunLine;

      // Process the line if in correct subrun.
      if ( (subrunNow == subrun || subrunNow == SUBRUNDEFAULT)
         && !readString( line, warn) ) accepted = false;
    }

  // Reached end of input file.
  };
  return accepted;

}

//--------------------------------------------------------------------------

// Routine to pass in pointers to PDF's. Usage optional.

// done
bool Pythia::setPDFPtr( PDF* pdfAPtrIn, PDF* pdfBPtrIn, PDF* pdfHardAPtrIn,
  PDF* pdfHardBPtrIn, PDF* pdfPomAPtrIn, PDF* pdfPomBPtrIn) {

  // Delete any PDF's created in a previous init call.
  if (useNewPdfHard && pdfHardAPtr != pdfAPtr) delete pdfHardAPtr;
  if (useNewPdfHard && pdfHardBPtr != pdfBPtr) delete pdfHardBPtr;
  if (useNewPdfA) delete pdfAPtr;
  if (useNewPdfB) delete pdfBPtr;
  if (useNewPdfPomA) delete pdfPomAPtr;
  if (useNewPdfPomB) delete pdfPomBPtr;

  // Reset pointers to be empty.
  useNewPdfA    = false;
  useNewPdfB    = false;
  useNewPdfHard = false;
  useNewPdfPomA = false;
  useNewPdfPomB = false;
  pdfAPtr       = 0;
  pdfBPtr       = 0;
  pdfHardAPtr   = 0;
  pdfHardBPtr   = 0;
  pdfPomAPtr    = 0;
  pdfPomBPtr    = 0;

  // Switch off external PDF's by zero as input.
  if (pdfAPtrIn == 0 && pdfBPtrIn == 0) return true;

  // The two PDF objects cannot be one and the same.
  if (pdfAPtrIn == pdfBPtrIn) return false;

  // Save pointers.
  pdfAPtr       = pdfAPtrIn;
  pdfBPtr       = pdfBPtrIn;

  // By default same pointers for hard-process PDF's.
  pdfHardAPtr   = pdfAPtrIn;
  pdfHardBPtr   = pdfBPtrIn;

  // Optionally allow separate pointers for hard process.
  if (pdfHardAPtrIn != 0 && pdfHardBPtrIn != 0) {
    if (pdfHardAPtrIn == pdfHardBPtrIn) return false;
    pdfHardAPtr = pdfHardAPtrIn;
    pdfHardBPtr = pdfHardBPtrIn;
  }

  // Optionally allow pointers for Pomerons in the proton.
  if (pdfPomAPtrIn != 0 && pdfPomBPtrIn != 0) {
    if (pdfPomAPtrIn == pdfPomBPtrIn) return false;
    pdfPomAPtr  = pdfPomAPtrIn;
    pdfPomBPtr  = pdfPomBPtrIn;
  }

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Routine to initialize with the variable values of the Beams kind.

// done
bool Pythia::init() {

  Benchmark_start(Pythia0init);
  Benchmark_start(Pythia0init_readBeamSettings); // trivial

  // Check that constructor worked.
  isInit = false;
  if (!isConstructed) {
    pState.info.errorMsg("Abort from Pythia::init: constructor "
      "initialization failed");
    return false;
  }

  // Early readout, if return false or changed when no beams.
  doProcessLevel = pState.settings.get(Flag::ProcessLevel_all);

  // Check that changes in Settings and ParticleData have worked.
  // if (pState.settings.readingFailed()) {
  //   pState.info.errorMsg("Abort from Pythia::init: some user settings "
  //     "did not make sense");
  //   return false;
  // }
  if (pState.particleData.readingFailed()) {
    pState.info.errorMsg("Abort from Pythia::init: some user particle data "
      "did not make sense");
    return false;
  }

  // Begin initialization. Find which frame type to use.
  pState.info.addCounter(1);
  frameType = mode("Beams:frameType");

  // Initialization with internal processes: read in and set values.
  if (frameType < 4 ) {
    doLHA     = false;
    boostType = frameType;
    idA       = mode("Beams:idA");
    idB       = mode("Beams:idB");
    eCM       = parm("Beams:eCM");
    eA        = parm("Beams:eA");
    eB        = parm("Beams:eB");
    pxA       = parm("Beams:pxA");
    pyA       = parm("Beams:pyA");
    pzA       = parm("Beams:pzA");
    pxB       = parm("Beams:pxB");
    pyB       = parm("Beams:pyB");
    pzB       = parm("Beams:pzB");

   // Initialization with a Les Houches Event File or an LHAup object.
  } else {
    doLHA     = true;
    boostType = 2;
    string lhef        = word("Beams:LHEF");
    string lhefHeader  = word("Beams:LHEFheader");
    bool   readHeaders = flag("Beams:readLHEFheaders");
    bool   setScales   = flag("Beams:setProductionScalesFromLHEF");
    bool   skipInit    = flag("Beams:newLHEFsameInit");
    int    nSkipAtInit = mode("Beams:nSkipLHEFatInit");

    // For file input: renew file stream or (re)new Les Houches object.
    if (frameType == 4) {
      const char* cstring1 = lhef.c_str();
      if (useNewLHA && skipInit) pState.lhaUp->newEventFile(cstring1);
      else {
        if (useNewLHA) delete pState.lhaUp;
        // Header is optional, so use NULL pointer to indicate no value.
        const char* cstring2 = (lhefHeader == "void")
          ? NULL : lhefHeader.c_str();
        pState.lhaUp   = new LHAupLHEF(&pState.info, cstring1, cstring2,
          readHeaders, setScales); // !!!
        useNewLHA  = true;
      }

      // Check that file was properly opened.
      if (!pState.lhaUp->fileFound()) {
        pState.info.errorMsg("Abort from Pythia::init: "
          "Les Houches Event File not found");
        return false;
      }

    // For object input: at least check that not null pointer.
    } else {
      if (pState.lhaUp == 0) {
        pState.info.errorMsg("Abort from Pythia::init: "
          "LHAup object not found");
        return false;
      }

      // LHAup object generic abort using fileFound() routine.
      if (!pState.lhaUp->fileFound()) {
        pState.info.errorMsg("Abort from Pythia::init: "
          "LHAup initialisation error");
        return false;
      }
    }

    // Send in pointer to pState.info. Store or replace LHA pointer in other classes.
    pState.lhaUp->setPtr( &pState.info);
    processLevel.setLHAPtr( pState.lhaUp);

    // If second time around, only with new file, then simplify.
    // Optionally skip ahead a number of events at beginning of file.
    if (skipInit) {
      isInit = true;
      pState.info.addCounter(2);
      if (nSkipAtInit > 0) pState.lhaUp->skipEvent(nSkipAtInit);
      return true;
    }

    // Set up values related to merging hooks.
    if (frameType == 4 || frameType == 5) {

      // Set up values related to CKKW-L merging.
      bool doUserMerging     = pState.settings.get(Flag::Merging_doUserMerging);
      bool doMGMerging       = pState.settings.get(Flag::Merging_doMGMerging);
      bool doKTMerging       = pState.settings.get(Flag::Merging_doKTMerging);
      bool doPTLundMerging   = pState.settings.get(Flag::Merging_doPTLundMerging);
      bool doCutBasedMerging = pState.settings.get(Flag::Merging_doCutBasedMerging);
      // Set up values related to unitarised CKKW merging
      bool doUMEPSTree       = pState.settings.get(Flag::Merging_doUMEPSTree);
      bool doUMEPSSubt       = pState.settings.get(Flag::Merging_doUMEPSSubt);
      // Set up values related to NL3 NLO merging
      bool doNL3Tree         = pState.settings.get(Flag::Merging_doNL3Tree);
      bool doNL3Loop         = pState.settings.get(Flag::Merging_doNL3Loop);
      bool doNL3Subt         = pState.settings.get(Flag::Merging_doNL3Subt);
      // Set up values related to unitarised NLO merging
      bool doUNLOPSTree      = pState.settings.get(Flag::Merging_doUNLOPSTree);
      bool doUNLOPSLoop      = pState.settings.get(Flag::Merging_doUNLOPSLoop);
      bool doUNLOPSSubt      = pState.settings.get(Flag::Merging_doUNLOPSSubt);
      bool doUNLOPSSubtNLO   = pState.settings.get(Flag::Merging_doUNLOPSSubtNLO);
      bool doXSectionEst     = pState.settings.get(Flag::Merging_doXSectionEstimate);
      doMerging = doUserMerging || doMGMerging || doKTMerging
        || doPTLundMerging || doCutBasedMerging || doUMEPSTree || doUMEPSSubt
        || doNL3Tree || doNL3Loop || doNL3Subt || doUNLOPSTree
        || doUNLOPSLoop || doUNLOPSSubt || doUNLOPSSubtNLO || doXSectionEst;

      // Set up MergingHooks object.
      bool inputMergingHooks = (mergingHooksPtr != 0);
      if (doMerging && !inputMergingHooks){
        if (hasOwnMergingHooks && mergingHooksPtr) delete mergingHooksPtr;
        mergingHooksPtr = new MergingHooks();
        hasOwnMergingHooks = true;
      }

      hasMergingHooks  = (mergingHooksPtr != 0);
      // Merging hooks required for merging. If no merging hooks pointer is
      // available, exit.
      if (doMerging && !hasMergingHooks) {
        pState.info.errorMsg("Abort from Pythia::init: "
          "no merging hooks object has been provided");
        return false;
      } else if (doMerging) {
        string lhefIn = (frameType == 4) ? lhef : "";
        mergingHooksPtr->setLHEInputFile( lhefIn);
      }

      // Initialise counting of Les Houches Events significantly above the
      // merging scale.
      pState.info.setCounter(41,0);
    }

    // Set LHAinit information (in some external program).
    if ( !pState.lhaUp->setInit()) {
      pState.info.errorMsg("Abort from Pythia::init: "
        "Les Houches initialization failed");
      return false;
    }

    // Extract beams from values set in an LHAinit object.
    idA = pState.lhaUp->idBeamA();
    idB = pState.lhaUp->idBeamB();
    int idRenameBeams = pState.settings.get(Mode::LesHouches_idRenameBeams);
    if (abs(idA) == idRenameBeams) idA = 16;
    if (abs(idB) == idRenameBeams) idB = -16;
    if (idA == 0 || idB == 0) doProcessLevel = false;
    eA  = pState.lhaUp->eBeamA();
    eB  = pState.lhaUp->eBeamB();

    // Optionally skip ahead a number of events at beginning of file.
    if (nSkipAtInit > 0) pState.lhaUp->skipEvent(nSkipAtInit);
  }

  Benchmark_stop(Pythia0init_readBeamSettings);
  Benchmark_start(Pythia0init_readGenSettings); // trivial

  // Set up values related to user hooks.
  hasUserHooks     = (userHooksPtr != 0);
  doVetoProcess    = false;
  doVetoPartons    = false;
  retryPartonLevel = false;
  if (hasUserHooks) {
    userHooksPtr->initPtr( &pState.info, &pState.settings, &pState.particleData, &pState.rndm, &beamA,
      &beamB, &beamPomA, &beamPomB, pState.couplings, &partonSystems, &sigmaTot);
    if (!userHooksPtr->initAfterBeams()) {
      pState.info.errorMsg("Abort from Pythia::init: could not initialise UserHooks");
      return false;
    }
    doVetoProcess    = userHooksPtr->canVetoProcessLevel();
    doVetoPartons    = userHooksPtr->canVetoPartonLevel();
    retryPartonLevel = userHooksPtr->retryPartonLevel();
  }

  // Back to common initialization. Reset error counters.
  nErrEvent = 0;
  pState.info.setTooLowPTmin(false);
  pState.info.sigmaReset();

  // Initialize data members extracted from database.
  doPartonLevel    = pState.settings.get(Flag::PartonLevel_all);
  doHadronLevel    = pState.settings.get(Flag::HadronLevel_all);
  doDiffraction    = pState.settings.get(Flag::SoftQCD_all)
                  || pState.settings.get(Flag::SoftQCD_singleDiffractive)
                  || pState.settings.get(Flag::SoftQCD_doubleDiffractive)
                  || pState.settings.get(Flag::SoftQCD_centralDiffractive)
                  || pState.settings.get(Flag::SoftQCD_inelastic);
  doHardDiff       = pState.settings.get(Flag::Diffraction_doHard);
  doResDec         = pState.settings.get(Flag::ProcessLevel_resonanceDecays);
  doFSRinRes       = doPartonLevel && pState.settings.get(Flag::PartonLevel_FSR)
                  && pState.settings.get(Flag::PartonLevel_FSRinResonances);
  decayRHadrons    = pState.settings.get(Flag::RHadrons_allowDecay);
  doMomentumSpread = pState.settings.get(Flag::Beams_allowMomentumSpread);
  doVertexSpread   = pState.settings.get(Flag::Beams_allowVertexSpread);
  abortIfVeto      = pState.settings.get(Flag::Check_abortIfVeto);
  checkEvent       = pState.settings.get(Flag::Check_event);
  checkHistory     = pState.settings.get(Flag::Check_history);
  nErrList         = pState.settings.get(Mode::Check_nErrList);
  epTolErr         = pState.settings.get(Param::Check_epTolErr);
  epTolWarn        = pState.settings.get(Param::Check_epTolWarn);
  mTolErr          = pState.settings.get(Param::Check_mTolErr);
  mTolWarn         = pState.settings.get(Param::Check_mTolWarn);

  // Initialise merging hooks.
  if ( doMerging && (hasMergingHooks || hasOwnMergingHooks) )
    mergingHooksPtr->init( pState.settings, &pState.info, &pState.particleData, &partonSystems );

  // Initialize the random number generator.
  if ( pState.settings.get(Flag::Random_setSeed) )
    pState.rndm.init( pState.settings.get(Mode::Random_seed) ); // !!!
  pState.rndm.init( 123 );

  // Check that combinations of settings are allowed; change if not.
  checkSettings();

  Benchmark_stop(Pythia0init_readGenSettings);
  Benchmark_start(Pythia0init_initSMgaugeCoup); // trivial

  // Initialize the SM couplings (needed to initialize resonances).
  pState.couplings->init( pState.settings, &pState.rndm );

  Benchmark_stop(Pythia0init_initSMgaugeCoup);
  Benchmark_start(Pythia0init_initSUSYCoup); // trivial

  // Initialize SLHA interface (including SLHA/BSM couplings).
  bool useSLHAcouplings = false;
  pState.slhaInterface.setPtr( &pState.info );

  pState.slhaInterface.init( pState.settings, &pState.rndm, pState.couplings, &pState.particleData,
    useSLHAcouplings, particleDataBuffer );
  if (useSLHAcouplings) pState.couplings = pState.slhaInterface.couplingsPtr;
  
  Benchmark_stop(Pythia0init_initSUSYCoup);
  Benchmark_start(Pythia0init_initResonanceWidths);

  // Reset couplingsPtr to the correct memory address.
  // (in case we are using SUSY coup rather than SM coup)
  pState.particleData.initPtr( &pState.info, &pState.settings, &pState.rndm, pState.couplings);
  if (hasUserHooks) userHooksPtr->initPtr( &pState.info, &pState.settings, &pState.particleData,
    &pState.rndm, &beamA, &beamB, &beamPomA, &beamPomB, pState.couplings,
    &partonSystems, &sigmaTot);

  // Set headers to distinguish the two event listing kinds.
  // (nothing major here; just setting the Event headers)
  int startColTag = pState.settings.get(Mode::Event_startColTag);
  process.init("(hard process)", &pState.particleData, startColTag);
  event.init("(complete event)", &pState.particleData, startColTag);

  // Final setup stage of particle data, notably resonance widths.
  // (just creates a list of pointers to all the Resonance base classes)
  pState.particleData.initWidths( resonancePtrs);

  Benchmark_stop(Pythia0init_initResonanceWidths);
  Benchmark_start(Pythia0init_setupRhadrons); // trivial

  // Set up R-hadrons particle data, where relevant.
  rHadrons.init( &pState.info, pState.settings, &pState.particleData, &pState.rndm);

  Benchmark_stop(Pythia0init_setupRhadrons);
  Benchmark_start(Pythia0init_initShowers); // trivial

  // Set up objects for timelike and spacelike showers.
  // Note: second timeShower for decays only (in case we want a different algorithm for this)
  // Note: the user can set there own Time/Space shower prior to this point
  if (timesDecPtr == 0 || timesPtr == 0) {
    TimeShower* timesNow = new TimeShower();
    if (timesDecPtr == 0) {
      timesDecPtr    = timesNow;
      useNewTimesDec = true;
    }
    if (timesPtr == 0) {
      timesPtr    = timesNow;
      useNewTimes = true;
    }
  }
  if (spacePtr == 0) {
    spacePtr    = new SpaceShower();
    useNewSpace = true;
  }

  // Initialize pointers in showers.
  timesPtr->initPtr( &pState.info, &pState.settings, &pState.particleData, &pState.rndm, pState.couplings,
    &partonSystems, userHooksPtr, mergingHooksPtr);
  timesDecPtr->initPtr( &pState.info, &pState.settings, &pState.particleData, &pState.rndm, pState.couplings,
    &partonSystems, userHooksPtr, mergingHooksPtr);
  spacePtr->initPtr( &pState.info, &pState.settings, &pState.particleData, &pState.rndm, pState.couplings,
    &partonSystems, userHooksPtr, mergingHooksPtr);

  Benchmark_stop(Pythia0init_initShowers);
  Benchmark_start(Pythia0init_initBeamShape); // trivial

  // Set up values related to beam shape.
  // (nothing major here; just copying settings)
  if (beamShapePtr == 0) {
    beamShapePtr    = new BeamShape();
    useNewBeamShape = true;
  }
  beamShapePtr->init( pState.settings, &pState.rndm);

  // Check that beams and beam combination can be handled.
  if (!checkBeams()) {
    pState.info.errorMsg("Abort from Pythia::init: "
      "checkBeams initialization failed");
    return false;
  }

  Benchmark_stop(Pythia0init_initBeamShape);
  Benchmark_start(Pythia0init_initBeams); // trivial

  // (again this is just copying settings for each beam)

  // Do not set up beam kinematics when no process level.
  if (!doProcessLevel) boostType = 1;
  else {

    // Set up beam kinematics.
    // (calcs the beam CM energy etc.) !!
    if (!initKinematics()) {
      pState.info.errorMsg("Abort from Pythia::init: "
        "kinematics initialization failed");
      return false;
    }

    // Set up pointers to PDFs.
    // gets the PDF derived class for each beam
    if (!initPDFs()) {
      pState.info.errorMsg("Abort from Pythia::init: PDF initialization failed");
      return false;
    }

    // Set up the two beams and the common remnant system.
    // Note: StringFlav class is used to select parton/hadron flavours
    // it belongs to hadronLevel, so need to extract it here
    StringFlav* flavSelPtr = hadronLevel.getStringFlavPtr();
    beamA.init( idA, pzAcm, eA, mA, &pState.info, pState.settings, &pState.particleData, &pState.rndm,
      pdfAPtr, pdfHardAPtr, isUnresolvedA, flavSelPtr);
    beamB.init( idB, pzBcm, eB, mB, &pState.info, pState.settings, &pState.particleData, &pState.rndm,
      pdfBPtr, pdfHardBPtr, isUnresolvedB, flavSelPtr);

    // Optionally set up new alternative beams for these Pomerons.
    if ( doDiffraction || doHardDiff) {
      beamPomA.init( 990,  0.5 * eCM, 0.5 * eCM, 0., &pState.info, pState.settings,
        &pState.particleData, &pState.rndm, pdfPomAPtr, pdfPomAPtr, false, flavSelPtr);
      beamPomB.init( 990, -0.5 * eCM, 0.5 * eCM, 0., &pState.info, pState.settings,
        &pState.particleData, &pState.rndm, pdfPomBPtr, pdfPomBPtr, false, flavSelPtr);
    }
  }

  Benchmark_stop(Pythia0init_initBeams);
  Benchmark_start(Pythia0init_initProcessLevel); // done

  // !!!
  // Send info/pointers to process level for initialization.
  if ( doProcessLevel && !processLevel.init(&pState, &pState.info, pState.settings, &pState.particleData,
    &pState.rndm, &beamA, &beamB, pState.couplings, &sigmaTot, doLHA, &pState.slhaInterface,
    userHooksPtr, sigmaPtrs, phaseSpacePtrs) ) {
    pState.info.errorMsg("Abort from Pythia::init: "
      "processLevel initialization failed");
    return false;
  }

  // Initialize timelike showers already here, since needed in decays.
  // The pointers to the beams are needed by some external plugin showers.
  // (not much here; just copying settings)
  timesDecPtr->init( &beamA, &beamB);

  // (note: above processLevel.init will do this)
  // Alternatively only initialize resonance decays.
  if ( !doProcessLevel) processLevel.initDecays( &pState.info, &pState.particleData,
    &pState.rndm, pState.lhaUp);

  Benchmark_stop(Pythia0init_initProcessLevel);
  Benchmark_start(Pythia0init_initPartonLevel); // do next...

  // Send info/pointers to parton level for initialization.
  if ( doPartonLevel && doProcessLevel && !partonLevel.init(&pState, &pState.info, pState.settings,
    &pState.particleData, &pState.rndm, &beamA, &beamB, &beamPomA, &beamPomB, pState.couplings,
    &partonSystems, &sigmaTot, timesDecPtr, timesPtr, spacePtr, &rHadrons,
    userHooksPtr, mergingHooksPtr, false) ) {
    pState.info.errorMsg("Abort from Pythia::init: "
      "partonLevel initialization failed" );
    return false;
  }

  Benchmark_stop(Pythia0init_initPartonLevel);
  Benchmark_start(Pythia0init_initTrialPartonLevel);

  // Make pointer to shower available for merging machinery.
  if ( doMerging && (hasMergingHooks || hasOwnMergingHooks) )
    mergingHooksPtr->setShowerPointer(&partonLevel);

  // Alternatively only initialize final-state showers in resonance decays.
  if ( !doProcessLevel || !doPartonLevel) partonLevel.init(&pState, &pState.info, pState.settings,
    &pState.particleData, &pState.rndm, 0, 0, 0, 0, pState.couplings, &partonSystems, 0,
    timesDecPtr, 0, 0, &rHadrons, 0, 0, false);

  // Send info/pointers to parton level for trial shower initialization.
  if ( doMerging && !trialPartonLevel.init(&pState, &pState.info, pState.settings, &pState.particleData,
      &pState.rndm, &beamA, &beamB, &beamPomA, &beamPomB, pState.couplings,
      &partonSystems, &sigmaTot, timesDecPtr, timesPtr, spacePtr, &rHadrons,
      userHooksPtr, mergingHooksPtr, true) ) {
    pState.info.errorMsg("Abort from Pythia::init: "
      "trialPartonLevel initialization failed");
    return false;
  }

  Benchmark_stop(Pythia0init_initTrialPartonLevel);
  Benchmark_start(Pythia0init_initHadronLevel);

  // Initialise the merging wrapper class.
  if (doMerging ) merging.init( &pState.settings, &pState.info, &pState.particleData, &pState.rndm,
    &beamA, &beamB, mergingHooksPtr, &trialPartonLevel );

  // !!!
  // Send info/pointers to hadron level for initialization.
  // Note: forceHadronLevel() can come, so we must always initialize.
  if ( !hadronLevel.init( &pState.info, pState.settings, &pState.particleData, &pState.rndm,
    pState.couplings, timesDecPtr, &rHadrons, decayHandlePtr,
    handledParticles, userHooksPtr) ) {
    pState.info.errorMsg("Abort from Pythia::init: "
      "hadronLevel initialization failed");
    return false;
  }

  Benchmark_stop(Pythia0init_initHadronLevel);
  Benchmark_start(Pythia0init_getMoreOptions);

  // Optionally check particle data table for inconsistencies.
  if ( pState.settings.get(Flag::Check_particleData) )
    pState.particleData.checkTable( pState.settings.get(Mode::Check_levelParticleData) );

  // Optionally show settings and particle data, changed or all.
  bool showCS  = pState.settings.get(Flag::Init_showChangedSettings);
  bool showAS  = pState.settings.get(Flag::Init_showAllSettings);
  bool showCPD = pState.settings.get(Flag::Init_showChangedParticleData);
  bool showCRD = pState.settings.get(Flag::Init_showChangedResonanceData);
  bool showAPD = pState.settings.get(Flag::Init_showAllParticleData);
  int  show1PD = pState.settings.get(Mode::Init_showOneParticleData);
  bool showPro = pState.settings.get(Flag::Init_showProcesses);
  if (showCS)      pState.settings.list(true);
  if (showAS)      pState.settings.list(false);
  if (show1PD > 0) pState.particleData.list(show1PD);
  if (showCPD)     pState.particleData.listChanged(showCRD);
  if (showAPD)     pState.particleData.listAll();

  // Listing options for the next() routine.
  nCount       = pState.settings.get(Mode::Next_numberCount);
  nShowLHA     = pState.settings.get(Mode::Next_numberShowLHA);
  nShowInfo    = pState.settings.get(Mode::Next_numberShowInfo);
  nShowProc    = pState.settings.get(Mode::Next_numberShowProcess);
  nShowEvt     = pState.settings.get(Mode::Next_numberShowEvent);
  showSaV      = pState.settings.get(Flag::Next_showScaleAndVertex);
  showMaD      = pState.settings.get(Flag::Next_showMothersAndDaughters);

  Benchmark_stop(Pythia0init_getMoreOptions);
  Benchmark_start(Pythia0init_initCOlorReconnection);

  // Init colour reconnection and junction splitting.
  colourReconnection.init( &pState.info, pState.settings, &pState.rndm, &pState.particleData,
    &beamA, &beamB, &partonSystems);

  Benchmark_stop(Pythia0init_initCOlorReconnection);
  Benchmark_start(Pythia0init_initCOlorJunctionSplitting);
  
  junctionSplitting.init(&pState.info, pState.settings, &pState.rndm, &pState.particleData);

  // Flags for colour reconnection.
  doReconnect        = pState.settings.get(Flag::ColourReconnection_reconnect);
  reconnectMode      = pState.settings.get(Mode::ColourReconnection_mode);
  forceHadronLevelCR = pState.settings.get(Flag::ColourReconnection_forceHadronLevelCR);

  // Succeeded.
  isInit = true;
  pState.info.addCounter(2);
  if (useNewLHA && showPro) pState.lhaUp->listInit();
  return true;

}

//--------------------------------------------------------------------------

// Check that combinations of settings are allowed; change if not.

// done
void Pythia::checkSettings() {

  // Double rescattering not allowed if ISR or FSR.
  if ((pState.settings.get(Flag::PartonLevel_ISR) || pState.settings.get(Flag::PartonLevel_FSR))
    && pState.settings.get(Flag::MultipartonInteractions_allowDoubleRescatter)) {
    pState.info.errorMsg("Warning in Pythia::checkSettings: "
        "double rescattering switched off since showering is on");
    pState.settings.set(Flag::MultipartonInteractions_allowDoubleRescatter, false);
  }

}

//--------------------------------------------------------------------------

// Check that beams and beam combination can be handled. Set up unresolved.

// done
bool Pythia::checkBeams() {

  // Absolute flavours. If not to do process level then no check needed.
  int idAabs = abs(idA);
  int idBabs = abs(idB);
  if (!doProcessLevel) return true;

  // Neutrino beams always unresolved, charged lepton ones conditionally.
  bool isLeptonA  = (idAabs > 10 && idAabs < 17);
  bool isLeptonB  = (idBabs > 10 && idBabs < 17);
  bool isUnresLep = !pState.settings.get(Flag::PDF_lepton);
  isUnresolvedA   = isLeptonA && (idAabs%2 == 0 || isUnresLep);
  isUnresolvedB   = isLeptonB && (idBabs%2 == 0 || isUnresLep);

  // Equate Dark Matter "beams" with incoming neutrinos.
  if (idAabs > 50 && idAabs < 61) isLeptonA = isUnresolvedA = true;
  if (idBabs > 50 && idBabs < 61) isLeptonB = isUnresolvedB = true;

  // Lepton-lepton collisions OK (including neutrinos) if both (un)resolved.
  if (isLeptonA && isLeptonB && isUnresolvedA == isUnresolvedB) return true;

  // MBR model only implemented for pp/ppbar/pbarp collisions.
  int PomFlux     = pState.settings.get(Mode::Diffraction_PomFlux);
  if (PomFlux == 5) {
    bool ispp       = (idAabs == 2212 && idBabs == 2212);
    bool ispbarpbar = (idA == -2212 && idB == -2212);
    if (ispp && !ispbarpbar) return true;
    pState.info.errorMsg("Error in Pythia::init: cannot handle this beam combination"
      " with PomFlux == 5");
    return false;
  }

  // Hadron-hadron collisions OK, with Pomeron counted as hadron.
  bool isHadronA = (idAabs == 2212) || (idAabs == 2112) || (idA == 111)
                || (idAabs == 211)  || (idA == 990);
  bool isHadronB = (idBabs == 2212) || (idBabs == 2112) || (idB == 111)
                || (idBabs == 211)  || (idB == 990);
  if (isHadronA && isHadronB) return true;

  // Lepton-hadron collisions OK for DIS processes or LHEF input,
  // although still primitive.
  if ( (isLeptonA && isHadronB) || (isHadronA && isLeptonB) ) {
    bool doDIS = pState.settings.get(Flag::WeakBosonExchange_all)
              || pState.settings.get(Flag::WeakBosonExchange_ff2ff_t_gmZ_)
              || pState.settings.get(Flag::WeakBosonExchange_ff2ff_t_W_)
              || (frameType == 4);
    if (doDIS) return true;
  }

  // If no case above then failed.
  pState.info.errorMsg("Error in Pythia::init: cannot handle this beam combination");
  return false;

}

//--------------------------------------------------------------------------

// Calculate kinematics at initialization. Store beam four-momenta.

// done
bool Pythia::initKinematics() {

  // Find masses. Initial guess that we are in CM frame.
  mA       = pState.particleData.m0(idA);
  mB       = pState.particleData.m0(idB);
  betaZ    = 0.;
  gammaZ   = 1.;

  // Collinear beams not in CM frame: find CM energy.
  if (boostType == 2) {
    eA     = max(eA, mA);
    eB     = max(eB, mB);
    pzA    = sqrt(eA*eA - mA*mA);
    pzB    = -sqrt(eB*eB - mB*mB);
    pAinit = Vec4( 0., 0., pzA, eA);
    pBinit = Vec4( 0., 0., pzB, eB);
    eCM    = sqrt( pow2(eA + eB) - pow2(pzA + pzB) );

    // Find boost to rest frame.
    betaZ  = (pzA + pzB) / (eA + eB);
    gammaZ = (eA + eB) / eCM;
    if (abs(betaZ) < 1e-10) boostType = 1;
  }

  // Completely general beam directions: find CM energy.
  else if (boostType == 3) {
    eA     = sqrt( pxA*pxA + pyA*pyA + pzA*pzA + mA*mA);
    eB     = sqrt( pxB*pxB + pyB*pyB + pzB*pzB + mB*mB);
    pAinit = Vec4( pxA, pyA, pzA, eA);
    pBinit = Vec4( pxB, pyB, pzB, eB);
    eCM = (pAinit + pBinit).mCalc();

    // Find boost+rotation needed to move from/to CM frame.
    MfromCM.reset();
    MfromCM.fromCMframe( pAinit, pBinit);
    MtoCM = MfromCM;
    MtoCM.invert();
  }

  // Fail if CM energy below beam masses.
  if (eCM < mA + mB) {
    pState.info.errorMsg("Error in Pythia::initKinematics: too low energy");
    return false;
  }

  // Set up CM-frame kinematics with beams along +-z axis.
  pzAcm    = 0.5 * sqrtpos( (eCM + mA + mB) * (eCM - mA - mB)
           * (eCM - mA + mB) * (eCM + mA - mB) ) / eCM;
  pzBcm    = -pzAcm;
  eA       = sqrt(mA*mA + pzAcm*pzAcm);
  eB       = sqrt(mB*mB + pzBcm*pzBcm);

  // If in CM frame then store beam four-vectors (else already done above).
  if (boostType != 2 && boostType != 3) {
    pAinit = Vec4( 0., 0., pzAcm, eA);
    pBinit = Vec4( 0., 0., pzBcm, eB);
  }

  // Store main info for access in process generation.
  pState.info.setBeamA( idA, pzAcm, eA, mA);
  pState.info.setBeamB( idB, pzBcm, eB, mB);
  pState.info.setECM( eCM);

  // Must allow for generic boost+rotation when beam momentum spread.
  if (doMomentumSpread) boostType = 3;

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Set up pointers to PDFs.

// done
bool Pythia::initPDFs() {

  // Delete any PDF's created in a previous init call.
  if (useNewPdfHard) {
    if (pdfHardAPtr != pdfAPtr) {
      delete pdfHardAPtr;
      pdfHardAPtr = 0;
    }
    if (pdfHardBPtr != pdfBPtr) {
      delete pdfHardBPtr;
      pdfHardBPtr = 0;
    }
    useNewPdfHard = false;
  }
  if (useNewPdfA) {
    delete pdfAPtr;
    useNewPdfA    = false;
    pdfAPtr       = 0;
  }
  if (useNewPdfB) {
    delete pdfBPtr;
    useNewPdfB    = false;
    pdfBPtr       = 0;
  }
  if (useNewPdfPomA) {
    delete pdfPomAPtr;
    useNewPdfPomA = false;
    pdfPomAPtr    = 0;
  }
  if (useNewPdfPomB) {
    delete pdfPomBPtr;
    useNewPdfPomB = false;
    pdfPomBPtr    = 0;
  }

  // Note: there are PDFs for both hard + non-hard 

  // Set up the PDF's, if not already done.
  if (pdfAPtr == 0) {
    // (this function returns the appropriate PDF base class)
    pdfAPtr     = getPDFPtr(idA);
    if (pdfAPtr == 0 || !pdfAPtr->isSetup()) {
      pState.info.errorMsg("Error in Pythia::init: "
        "could not set up PDF for beam A");
      return false;
    }
    pdfHardAPtr = pdfAPtr;
    useNewPdfA  = true;
  }
  if (pdfBPtr == 0) {
    pdfBPtr     = getPDFPtr(idB, 1, "B");
    if (pdfBPtr == 0 || !pdfBPtr->isSetup()) {
      pState.info.errorMsg("Error in Pythia::init: "
        "could not set up PDF for beam B");
      return false;
    }
    pdfHardBPtr = pdfBPtr;
    useNewPdfB  = true;
  }

  // Optionally set up separate PDF's for hard process.
  if (pState.settings.get(Flag::PDF_useHard) && useNewPdfA && useNewPdfB) {
    pdfHardAPtr = getPDFPtr(idA, 2);
    if (!pdfHardAPtr->isSetup()) return false;
    pdfHardBPtr = getPDFPtr(idB, 2, "B");
    if (!pdfHardBPtr->isSetup()) return false;
    useNewPdfHard = true;
  }

  // Optionally set up Pomeron PDF's for diffractive physics.
  if ( doDiffraction || doHardDiff) {
    if (pdfPomAPtr == 0) {
      pdfPomAPtr    = getPDFPtr(990);
      useNewPdfPomA = true;
    }
    if (pdfPomBPtr == 0) {
      pdfPomBPtr    = getPDFPtr(990);
      useNewPdfPomB = true;
    }
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Main routine to generate the next event, using internal machinery.


bool Pythia::next() {

  static int nEvents = 0;
  nEvents++;

  {

  Benchmark_start(Pythia0next);
  Benchmark_start(Pythia0next_setup);

  // Check that constructor worked.
  if (!isConstructed) return false;

  // Regularly print how many events have been generated.
  int nPrevious = pState.info.getCounter(3);
  if (nCount > 0 && nPrevious > 0 && nPrevious%nCount == 0)
    cout << "\n Pythia::next(): " << nPrevious
         << " events have been generated " << endl;

  // Set/reset info counters specific to each event.
  pState.info.addCounter(3);
  for (int i = 10; i < 13; ++i) pState.info.setCounter(i);

  // Simpler option when no hard process, i.e. mainly hadron level.
  if (!doProcessLevel) {

    // Optionally fetch in resonance decays from LHA interface.
    if (doLHA && !processLevel.nextLHAdec( event)) {
      if (pState.info.atEndOfFile()) pState.info.errorMsg("Abort from Pythia::next:"
        " reached end of Les Houches Events File");
      return false;
    }

    // Reset info array (while event record contains data).
    pState.info.clear();

    // Set correct energy for system.
    Vec4 pSum = 0.;
    for (int i = 1; i < event.size(); ++i)
      if (event[i].isFinal()) pSum += event[i].p();
    event[0].p( pSum );
    event[0].m( pSum.mCalc() );

    // Generate hadronization and decays.
    bool status = forceHadronLevel();
    if (status) pState.info.addCounter(4);
    if (status && nPrevious < nShowEvt) event.list(showSaV, showMaD);
    return status;
  }

  // Reset arrays.
  pState.info.clear();
  process.clear();
  event.clear();
  partonSystems.clear();
  beamA.clear();
  beamB.clear();
  beamPomA.clear();
  beamPomB.clear();

  Benchmark_stop(Pythia0next_setup);
  Benchmark_start(Pythia0next_pickValenceFlavours);

  // Pick current beam valence flavours (for pi0, K0S, K0L, Pomeron).
  beamA.newValenceContent();
  beamB.newValenceContent();
  if ( doDiffraction || doHardDiff) {
    beamPomA.newValenceContent();
    beamPomB.newValenceContent();
  }

  // Can only generate event if initialization worked.
  if (!isInit) {
    pState.info.errorMsg("Abort from Pythia::next: "
      "not properly initialized so cannot generate events");
    return false;
  }

  // Pick beam momentum spread and beam vertex.
  if (doMomentumSpread || doVertexSpread) beamShapePtr->pick();

  // Recalculate kinematics when beam momentum spread.
  if (doMomentumSpread) nextKinematics();

  Benchmark_stop(Pythia0next_pickValenceFlavours);

  Benchmark_loopStart(Pythia0next_loopVetoHardProcess);
  bool quitLoop = false;

  // Outer loop over hard processes; only relevant for user-set vetoes.
  for ( ; ; ) {

    Benchmark_loopCount(Pythia0next_loopVetoHardProcess);
    Benchmark_start(Pythia0next_setLHEF3EventInfo);

    pState.info.addCounter(10);
    bool hasVetoed = false;
    bool hasVetoedDiff = false;

    // Provide the hard process that starts it off. Only one try.
    pState.info.clear();

    // Reset the event information. Necessary if the previous event was read
    // from LHEF, while the current event is not read from LHEF.
    pState.info.setLHEF3EventInfo();
    process.clear();

    Benchmark_stop(Pythia0next_setLHEF3EventInfo);
    Benchmark_start(Pythia0next_processLevel);

    if ( !processLevel.next( process) ) {
      if (doLHA && pState.info.atEndOfFile()) pState.info.errorMsg("Abort from "
        "Pythia::next: reached end of Les Houches Events File");
      else pState.info.errorMsg("Abort from Pythia::next: "
        "processLevel failed; giving up");
      {quitLoop = true; break;};
    }

    Benchmark_stop(Pythia0next_processLevel);
    Benchmark_start(Pythia0next_copyProcess);

    pState.info.addCounter(11);

    // Possibility for a user veto of the process-level event.
    if (doVetoProcess) {
      hasVetoed = userHooksPtr->doVetoProcessLevel( process);
      if (hasVetoed) {
        if (abortIfVeto) {quitLoop = true; break;};
        continue;
      }
    }

    // Possibility to perform matrix element merging for this event.
    if (doMerging) {
      int veto = merging.mergeProcess( process );
      // Apply possible merging scale cut.
      if (veto == -1) {
        hasVetoed = true;
        if (abortIfVeto) {quitLoop = true; break;};
        continue;
      // Exit because of vanishing no-emission probability.
      } else if (veto == 0) {
        event = process;
        break;
      }

      // Redo resonance decays after the merging, in case the resonance
      // structure has been changed because of reclusterings.
      if (veto == 2 && doResDec) processLevel.nextDecays( process);
    }

    // Possibility to stop the generation at this stage.
    if (!doPartonLevel) {
      boostAndVertex( true, true);
      processLevel.accumulate();
      pState.info.addCounter(4);
      if (doLHA && nPrevious < nShowLHA) pState.lhaUp->listEvent();
      if (nPrevious < nShowInfo) pState.info.list();
      if (nPrevious < nShowProc) process.list(showSaV, showMaD);
      return true;
    }

    // Save spare copy of process record in case of problems.
    Event processSave = process;
    int sizeMPI       = pState.info.sizeMPIarrays();
    pState.info.addCounter(12);
    for (int i = 14; i < 19; ++i) pState.info.setCounter(i);

    Benchmark_stop(Pythia0next_copyProcess);

    Benchmark_loopStart(Pythia0next_loopTrialPartonHadron);

    // Allow up to ten tries for parton- and hadron-level processing.
    bool physical   = true;
    for (int iTry = 0; iTry < NTRY; ++iTry) {
      
      Benchmark_loopCount(Pythia0next_loopTrialPartonHadron);
      Benchmark_start(Pythia0next_restoreProcess);

      pState.info.addCounter(14);
      physical  = true;
      hasVetoed = false;

      // Restore original process record if problems.
      if (iTry > 0) process = processSave;
      if (iTry > 0) pState.info.resizeMPIarrays( sizeMPI);

      // Reset event record and (extracted partons from) beam remnants.
      event.clear();
      beamA.clear();
      beamB.clear();
      beamPomA.clear();
      beamPomB.clear();
      partonSystems.clear();

      Benchmark_stop(Pythia0next_restoreProcess);

      Benchmark_start(Pythia0next_partonLevel);

      // Parton-level evolution: ISR, FSR, MPI.
      if ( !partonLevel.next( process, event) ) {

        // Abort event generation if parton level is set to abort.
        if (pState.info.getAbortPartonLevel()) {quitLoop = true; break;};

        // Skip to next hard process for failure owing to deliberate veto,
        // or alternatively retry for the same hard process.
        hasVetoed = partonLevel.hasVetoed();
        if (hasVetoed) {
          if (retryPartonLevel) {
            --iTry;
            continue;
          }
          if (abortIfVeto) {quitLoop = true; break;};
          break;
        }

        // If hard diffractive event has been discarded retry partonLevel.
        hasVetoedDiff = partonLevel.hasVetoedDiff();
        if (hasVetoedDiff) {
          pState.info.errorMsg("Warning in Pythia::next: "
            "discarding hard diffractive event from partonLevel; try again");
          break;
        }

        // Else make a new try for other failures.
        pState.info.errorMsg("Error in Pythia::next: "
          "partonLevel failed; try again");
        physical = false;
        continue;
      }
      pState.info.addCounter(15);

      Benchmark_stop(Pythia0next_partonLevel);

      Benchmark_start(Pythia0next_boostAndVertex);

      // Possibility for a user veto of the parton-level event.
      if (doVetoPartons) {
        hasVetoed = userHooksPtr->doVetoPartonLevel( event);
        if (hasVetoed) {
          if (abortIfVeto) return false;
          break;
        }
      }

      // Boost to lab frame (before decays, for vertices).
      boostAndVertex( true, true);

      // Possibility to stop the generation at this stage.
      if (!doHadronLevel) {
        processLevel.accumulate();
        partonLevel.accumulate();
        // Optionally check final event for problems.
        if (checkEvent && !check()) {
          pState.info.errorMsg("Abort from Pythia::next: "
            "check of event revealed problems");
          {quitLoop = true; break;}
        }
        pState.info.addCounter(4);
        if (doLHA && nPrevious < nShowLHA) pState.lhaUp->listEvent();
        if (nPrevious < nShowInfo) pState.info.list();
        if (nPrevious < nShowProc) process.list(showSaV, showMaD);
        if (nPrevious < nShowEvt)  event.list(showSaV, showMaD);
        return true;
      }

      Benchmark_stop(Pythia0next_boostAndVertex);

      Benchmark_start(Pythia0next_hadronLevel);

      // Hadron-level: hadronization, decays.
      pState.info.addCounter(16);
      if ( !hadronLevel.next( event) ) {
        pState.info.errorMsg("Error in Pythia::next: "
          "hadronLevel failed; try again");
        physical = false;
        continue;
      }

      Benchmark_stop(Pythia0next_hadronLevel);

      Benchmark_start(Pythia0next_doRHadronDecays);

      // If R-hadrons have been formed, then (optionally) let them decay.
      if (decayRHadrons && rHadrons.exist() && !doRHadronDecays()) {
        pState.info.errorMsg("Error in Pythia::next: "
          "decayRHadrons failed; try again");
        physical = false;
        continue;
      }
      pState.info.addCounter(17);

      Benchmark_stop(Pythia0next_doRHadronDecays);

      Benchmark_start(Pythia0next_checkForProblems);

      // Optionally check final event for problems.
      if (checkEvent && !check()) {
        pState.info.errorMsg("Error in Pythia::next: "
          "check of event revealed problems");
        physical = false;
        continue;
      }

      // Stop parton- and hadron-level looping if you got this far.
      pState.info.addCounter(18);

      Benchmark_stop(Pythia0next_checkForProblems);

      break;
    }

    if (quitLoop) break;

    Benchmark_loopStop(Pythia0next_loopTrialPartonHadron);

    // If event vetoed then to make a new try.
    if (hasVetoed || hasVetoedDiff)  {
      if (abortIfVeto) {quitLoop = true; break;}
      continue;
    }

    // If event failed any other way (after ten tries) then give up.
    if (!physical) {
      pState.info.errorMsg("Abort from Pythia::next: "
        "parton+hadronLevel failed; giving up");
      {quitLoop = true; break;}
    }


    // Process- and parton-level statistics. Event scale.
    Benchmark_start(Pythia0next_accumulateProcess);
    processLevel.accumulate();
    Benchmark_stop(Pythia0next_accumulateProcess);
    Benchmark_start(Pythia0next_accumulateParton);
    partonLevel.accumulate();
    Benchmark_stop(Pythia0next_accumulateParton);
    Benchmark_start(Pythia0next_scale);
    event.scale( process.scale() );
    Benchmark_stop(Pythia0next_scale);


    // End of outer loop over hard processes. Done with normal option.
    pState.info.addCounter(13);
    break;
  }

  Benchmark_loopStop(Pythia0next_loopVetoHardProcess);

  if (quitLoop) return false;

  Benchmark_start(Pythia0next_listEvents);

  // List events.
  if (doLHA && nPrevious < nShowLHA) pState.lhaUp->listEvent();
  if (nPrevious < nShowInfo) pState.info.list();
  if (nPrevious < nShowProc) process.list(showSaV,showMaD);
  if (nPrevious < nShowEvt)  event.list(showSaV, showMaD);

  // Done.
  pState.info.addCounter(4);

  Benchmark_stop(Pythia0next_listEvents);
  Benchmark_stop(Pythia0next);

  // ---------------------------------------

  }

  

  if (nEvents == 100)
  {
    
  }

  // --------------------------------------

  return true;

}

//--------------------------------------------------------------------------

// Generate only the hadronization/decay stage, using internal machinery.
// The "event" instance should already contain a parton-level configuration.

bool Pythia::forceHadronLevel(bool findJunctions) {

  // Can only generate event if initialization worked.
  if (!isInit) {
    pState.info.errorMsg("Abort from Pythia::forceHadronLevel: "
      "not properly initialized so cannot generate events");
    return false;
  }

  // Check whether any junctions in system. (Normally done in ProcessLevel.)
  // Avoid it if there are no final-state coloured partons.
  if (findJunctions) {
    event.clearJunctions();
    for (int i = 0; i < event.size(); ++i)
    if (event[i].isFinal()
    && (event[i].col() != 0 || event[i].acol() != 0)) {
      processLevel.findJunctions( event);
      break;
    }
  }

  // Allow for CR before the hadronization.
  if (forceHadronLevelCR) {

    // Setup parton system for SK-I and SK-II colour reconnection.
    // Require all final state particles to have the Ws as mothers.
    if (reconnectMode == 3 || reconnectMode == 4) {
      partonSystems.clear();
      partonSystems.addSys();
      partonSystems.addSys();
      for (int i = 5;i < event.size();++i) {
        if (event[i].mother1() - 3 < 0 || event[i].mother1() - 3 > 1) {
          pState.info.errorMsg("Error from Pythia::forceHadronLevel: "
            " Event is not setup correctly for SK-I or SK-II CR");
          return false;
        }
        partonSystems.addOut(event[i].mother1() - 3,i);
      }
    }

    // save spare copy of event in case of failure.
    Event spareEvent = event;
    bool colCorrect = false;

    // Allow up to ten tries for CR.
    for (int iTry = 0; iTry < NTRY; ++ iTry) {
      colourReconnection.next(event, 0);
      if (junctionSplitting.checkColours(event)) {
        colCorrect = true;
        break;
      }
      else event = spareEvent;
    }

    if (!colCorrect) {
      pState.info.errorMsg("Error in Pythia::forceHadronLevel: "
        "Colour reconnection failed.");
      return false;
    }
  }

  // Save spare copy of event in case of failure.
  Event spareEvent = event;

  // Allow up to ten tries for hadron-level processing.
  bool physical = true;
  for (int iTry = 0; iTry < NTRY; ++ iTry) {
    physical = true;

    // Check whether any resonances need to be handled at process level.
    if (doResDec) {
      process = event;
      processLevel.nextDecays( process);

      // Allow for showers if decays happened at process level.
      if (process.size() > event.size()) {
        if (doFSRinRes) {
          partonLevel.setupShowerSys( process, event);
          partonLevel.resonanceShowers( process, event, false);
        } else event = process;
      }
    }

    // Hadron-level: hadronization, decays.
    if (hadronLevel.next( event)) break;

    // If failure then warn, restore original configuration and try again.
    pState.info.errorMsg("Error in Pythia::forceHadronLevel: "
      "hadronLevel failed; try again");
    physical = false;
    event    = spareEvent;
  }

  // Done for simpler option.
  if (!physical)  {
    pState.info.errorMsg("Abort from Pythia::forceHadronLevel: "
      "hadronLevel failed; giving up");
    return false;
  }

  // Optionally check final event for problems.
  if (checkEvent && !check()) {
    pState.info.errorMsg("Abort from Pythia::forceHadronLevel: "
      "check of event revealed problems");
    return false;
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Recalculate kinematics for each event when beam momentum has a spread.

void Pythia::nextKinematics() {

  // Read out momentum shift to give current beam momenta.
  pAnow = pAinit + beamShapePtr->deltaPA();
  pAnow.e( sqrt(pAnow.pAbs2() + mA * mA) );
  pBnow = pBinit + beamShapePtr->deltaPB();
  pBnow.e( sqrt(pBnow.pAbs2() + mB * mB) );

  // Construct CM frame kinematics.
  eCM   = (pAnow + pBnow).mCalc();
  pzAcm = 0.5 * sqrtpos( (eCM + mA + mB) * (eCM - mA - mB)
        * (eCM - mA + mB) * (eCM + mA - mB) ) / eCM;
  pzBcm = -pzAcm;
  eA    = sqrt(mA*mA + pzAcm*pzAcm);
  eB    = sqrt(mB*mB + pzBcm*pzBcm);

  // Set relevant info for other classes to use.
  pState.info.setBeamA( idA, pzAcm, eA, mA);
  pState.info.setBeamB( idB, pzBcm, eB, mB);
  pState.info.setECM( eCM);
  beamA.newPzE( pzAcm, eA);
  beamB.newPzE( pzBcm, eB);

  // Set boost/rotation matrices from/to CM frame.
  MfromCM.reset();
  MfromCM.fromCMframe( pAnow, pBnow);
  MtoCM = MfromCM;
  MtoCM.invert();

}

//--------------------------------------------------------------------------

// Boost from CM frame to lab frame, or inverse. Set production vertex.

// done
void Pythia::boostAndVertex( bool toLab, bool setVertex) {

  // Boost process from CM frame to lab frame.
  if (toLab) {
    if      (boostType == 2) process.bst(0., 0., betaZ, gammaZ);
    else if (boostType == 3) process.rotbst(MfromCM);

    // Boost nonempty event from CM frame to lab frame.
    if (event.size() > 0) {
      if      (boostType == 2) event.bst(0., 0., betaZ, gammaZ);
      else if (boostType == 3) event.rotbst(MfromCM);
    }

  // Boost process from lab frame to CM frame.
  } else {
    if      (boostType == 2) process.bst(0., 0., -betaZ, gammaZ);
    else if (boostType == 3) process.rotbst(MtoCM);

    // Boost nonempty event from lab frame to CM frame.
    if (event.size() > 0) {
      if      (boostType == 2) event.bst(0., 0., -betaZ, gammaZ);
      else if (boostType == 3) event.rotbst(MtoCM);
    }
  }

  // Set production vertex; assumes particles are in lab frame and at origin.
  if (setVertex && doVertexSpread) {
    Vec4 vertex = beamShapePtr->vertex();
    for (int i = 0; i < process.size(); ++i) process[i].vProd( vertex);
    for (int i = 0; i < event.size(); ++i) event[i].vProd( vertex);
  }

}

//--------------------------------------------------------------------------

// Perform R-hadron decays, either as part of normal evolution or forced.

bool Pythia::doRHadronDecays( ) {

  // Check if R-hadrons exist to be processed.
  if ( !rHadrons.exist() ) return true;

  // Do the R-hadron decay itself.
  if ( !rHadrons.decay( event) ) return false;

  // Perform showers in resonance decay chains.
  if ( !partonLevel.resonanceShowers( process, event, false) ) return false;

  // Subsequent hadronization and decays.
  if ( !hadronLevel.next( event) ) return false;

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Print statistics on event generation.

void Pythia::stat() {

  // Read out settings for what to include.
  bool showPrL = pState.settings.get(Flag::Stat_showProcessLevel);
  bool showPaL = pState.settings.get(Flag::Stat_showPartonLevel);
  bool showErr = pState.settings.get(Flag::Stat_showErrors);
  bool reset   = pState.settings.get(Flag::Stat_reset);

  // Statistics on cross section and number of events.
  if (doProcessLevel) {
    if (showPrL) processLevel.statistics(false);
    if (reset)   processLevel.resetStatistics();
  }

  // Statistics from other classes, currently multiparton interactions.
  if (showPaL) partonLevel.statistics(false);
  if (reset)   partonLevel.resetStatistics();

  // Merging statistics.
  if (doMerging) merging.statistics();

  // Summary of which and how many warnings/errors encountered.
  if (showErr) pState.info.errorStatistics();
  if (reset)   pState.info.errorReset();

}

//--------------------------------------------------------------------------

// Write the Pythia banner, with symbol and version information.

// done
void Pythia::banner(ostream& os) {

  // Read in version number and last date of change.
  double versionNumber = pState.settings.get(Param::Pythia_versionNumber);
  int versionDate = pState.settings.get(Mode::Pythia_versionDate);
  string month[12] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};

  // Get date and time.
  time_t t = time(0);
  char dateNow[12];
  strftime(dateNow,12,"%d %b %Y",localtime(&t));
  char timeNow[9];
  strftime(timeNow,9,"%H:%M:%S",localtime(&t));

  os << "\n"
     << " *-------------------------------------------"
     << "-----------------------------------------* \n"
     << " |                                           "
     << "                                         | \n"
     << " |  *----------------------------------------"
     << "--------------------------------------*  | \n"
     << " |  |                                        "
     << "                                      |  | \n"
     << " |  |                                        "
     << "                                      |  | \n"
     << " |  |   PPP   Y   Y  TTTTT  H   H  III    A  "
     << "    Welcome to the Lund Monte Carlo!  |  | \n"
     << " |  |   P  P   Y Y     T    H   H   I    A A "
     << "    This is PYTHIA version " << fixed << setprecision(3)
     << setw(5) << versionNumber << "      |  | \n"
     << " |  |   PPP     Y      T    HHHHH   I   AAAAA"
     << "    Last date of change: " << setw(2) << versionDate%100
     << ' ' << month[ (versionDate/100)%100 - 1 ]
     << ' ' << setw(4) << versionDate/10000 <<  "  |  | \n"
     << " |  |   P       Y      T    H   H   I   A   A"
     << "                                      |  | \n"
     << " |  |   P       Y      T    H   H  III  A   A"
     << "    Now is " << dateNow << " at " << timeNow << "    |  | \n"
     << " |  |                                        "
     << "                                      |  | \n"
     << " |  |   Torbjorn Sjostrand;  Department of As"
     << "tronomy and Theoretical Physics,      |  | \n"
     << " |  |      Lund University, Solvegatan 14A, S"
     << "E-223 62 Lund, Sweden;                |  | \n"
     << " |  |      e-mail: torbjorn@thep.lu.se       "
     << "                                      |  | \n"
     << " |  |   Jesper Roy Christiansen;  Department "
     << "of Astronomy and Theoretical Physics, |  | \n"
     << " |  |      Lund University, Solvegatan 14A, S"
     << "E-223 62 Lund, Sweden;                |  | \n"
     << " |  |      e-mail: Jesper.Roy.Christiansen@th"
     << "ep.lu.se                              |  | \n"
     << " |  |   Nishita Desai;  Institut fuer Theoret"
     << "ische Physik,                         |  | \n"
     << " |  |     Universitaet Heidelberg, Philosophe"
     << "nweg 16, D-69120 Heidelberg, Germany; |  | \n"
     << " |  |      e-mail: n.desai@thphys.uni-heidelb"
     << "erg.de                                |  | \n"
     << " |  |   Philip Ilten;  Massachusetts Institut"
     << "e of Technology,                      |  | \n"
     << " |  |      stationed at CERN, CH-1211 Geneva "
     << "23, Switzerland;                      |  | \n"
     << " |  |      e-mail: philten@cern.ch           "
     << "                                      |  | \n"
     << " |  |   Stephen Mrenna;  Computing Division, "
     << "Simulations Group,                    |  | \n"
     << " |  |      Fermi National Accelerator Laborat"
     << "ory, MS 234, Batavia, IL 60510, USA;  |  | \n"
     << " |  |      e-mail: mrenna@fnal.gov           "
     << "                                      |  | \n"
     << " |  |   Stefan Prestel;  Theoretical Physics "
     << "Group,                                |  | \n"
     << " |  |      SLAC National Accelerator Laborato"
     << "ry, Menlo Park, CA 94025, USA;        |  | \n"
     << " |  |      e-mail: prestel@slac.stanford.edu "
     << "                                      |  | \n"
     << " |  |   Christine O. Rasmussen;  Department o"
     << "f Astronomy and Theoretical Physics,  |  | \n"
     << " |  |      Lund University, Solvegatan 14A, S"
     << "E-223 62 Lund, Sweden;                |  | \n"
     << " |  |      e-mail: christine.rasmussen@thep.l"
     << "u.se                                  |  | \n"
     << " |  |   Peter Skands;  School of Physics,    "
     << "                                      |  | \n"
     << " |  |      Monash University, PO Box 27, 3800"
     << " Melbourne, Australia;                |  | \n"
     << " |  |      e-mail: peter.skands@monash.edu   "
     << "                                      |  | \n"
     << " |  |                                        "
     << "                                      |  | \n"
     << " |  |   The main program reference is 'An Int"
     << "roduction to PYTHIA 8.2',             |  | \n"
     << " |  |   T. Sjostrand et al, Comput. Phys. Com"
     << "mun. 191 (2005) 159                   |  | \n"
     << " |  |   [arXiv:1410.3012 [hep-ph]]           "
     << "                                      |  | \n"
     << " |  |                                        "
     << "                                      |  | \n"
     << " |  |   The main physics reference is the 'PY"
     << "THIA 6.4 Physics and Manual',         |  | \n"
     << " |  |   T. Sjostrand, S. Mrenna and P. Skands"
     << ", JHEP05 (2006) 026 [hep-ph/0603175]  |  | \n"
     << " |  |                                        "
     << "                                      |  | \n"
     << " |  |   An archive of program versions and do"
     << "cumentation is found on the web:      |  | \n"
     << " |  |   http://www.thep.lu.se/Pythia         "
     << "                                      |  | \n"
     << " |  |                                        "
     << "                                      |  | \n"
     << " |  |   This program is released under the GN"
     << "U General Public Licence version 2.   |  | \n"
     << " |  |   Please respect the MCnet Guidelines f"
     << "or Event Generator Authors and Users. |  | \n"
     << " |  |                                        "
     << "                                      |  | \n"
     << " |  |   Disclaimer: this program comes withou"
     << "t any guarantees.                     |  | \n"
     << " |  |   Beware of errors and use common sense"
     << " when interpreting results.           |  | \n"
     << " |  |                                        "
     << "                                      |  | \n"
     << " |  |   Copyright (C) 2015 Torbjorn Sjostrand"
     << "                                      |  | \n"
     << " |  |                                        "
     << "                                      |  | \n"
     << " |  |                                        "
     << "                                      |  | \n"
     << " |  *----------------------------------------"
     << "--------------------------------------*  | \n"
     << " |                                           "
     << "                                         | \n"
     << " *-------------------------------------------"
     << "-----------------------------------------* \n" << endl;

}

//--------------------------------------------------------------------------

// Check for lines in file that mark the beginning of new subrun.

int Pythia::readSubrun(stringref line, bool warn, ostream& os) {

  // If empty line then done.
  int subrunLine = SUBRUNDEFAULT;
  if (line.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos)
    return subrunLine;

  // If first character is not a letter, then done.
  string lineNow = line;
  int firstChar = lineNow.find_first_not_of(" \n\t\v\b\r\f\a");
  if (!isalpha(lineNow[firstChar])) return subrunLine;

  // Replace an equal sign by a blank to make parsing simpler.
  while (lineNow.find("=") != string::npos) {
    int firstEqual = lineNow.find_first_of("=");
    lineNow.replace(firstEqual, 1, " ");
  }

  // Get first word of a line.
  istringstream splitLine(lineNow);
  string name;
  splitLine >> name;

  // Replace two colons by one (:: -> :) to allow for such mistakes.
  while (name.find("::") != string::npos) {
    int firstColonColon = name.find_first_of("::");
    name.replace(firstColonColon, 2, ":");
  }


  // Convert to lowercase.
  for (int i = 0; i < int(name.length()); ++i) name[i] = tolower(name[i]);

  // If no match then done.
  if (name != "main:subrun") return subrunLine;

  // Else find new subrun number and return it.
  splitLine >> subrunLine;
  if (!splitLine) {
    if (warn) os << "\n PYTHIA Warning: Main:subrun number not"
        << " recognized; skip:\n   " << line << endl;
    subrunLine = SUBRUNDEFAULT;
  }
  return subrunLine;

}

//--------------------------------------------------------------------------

// Check for lines in file that mark the beginning or end of commented section.
// Return +1 for beginning, -1 for end, 0 else.

// done
int Pythia::readCommented(stringref line) {

  // If less than two nontrivial characters on line then done.
  if (line.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos) return 0;
  int firstChar = line.find_first_not_of(" \n\t\v\b\r\f\a");
  if (int(line.size()) < firstChar + 2) return 0;

  // If first two nontrivial characters are /* or */ then done.
  if (line.substr(firstChar, 2) == "/*") return +1;
  if (line.substr(firstChar, 2) == "*/") return -1;

  // Else done.
  return 0;

}

//--------------------------------------------------------------------------

// Check that the final event makes sense: no unknown id codes;
// charge and energy-momentum conserved.

// done
bool Pythia::check(ostream& os) {

  Benchmark_start(Pythia0check);
  Benchmark_start(Pythia0check_reset);

  // Reset.
  bool physical     = true;
  bool listVertices = false;
  bool listHistory  = false;
  bool listSystems  = false;
  bool listBeams    = false;
  iErrId.resize(0);
  iErrCol.resize(0);
  iErrEpm.resize(0);
  iErrNan.resize(0);
  iErrNanVtx.resize(0);
  Vec4 pSum;
  double chargeSum  = 0.;

  Benchmark_stop(Pythia0check_reset);
  Benchmark_start(Pythia0check_charge_and_psum);

  // Incoming beams counted with negative momentum and charge.
  if (doProcessLevel) {
    pSum      = - (event[1].p() + event[2].p());
    chargeSum = - (event[1].charge() + event[2].charge());

  // If no ProcessLevel then sum final state of process record.
  } else if (process.size() > 0) {
    pSum = - process[0].p();
    for (int i = 0; i < process.size(); ++i)
      if (process[i].isFinal()) chargeSum -= process[i].charge();

  // If process not filled, then use outgoing primary in event.
  } else {
    pSum = - event[0].p();
    for (int i = 1; i < event.size(); ++i)
      if (event[i].statusAbs() < 10 || event[i].statusAbs() == 23)
        chargeSum -= event[i].charge();
  }
  double eLab = abs(pSum.e());

  Benchmark_stop(Pythia0check_charge_and_psum);
  Benchmark_loopStart(Pythia0check_loopAllParticles);

  // Loop over particles in the event.
  for (int i = 0; i < event.size(); ++i) {

    Benchmark_loopCount(Pythia0check_loopAllParticles);
    Benchmark_start(Pythia0check_find_unrecognized);

    // Look for any unrecognized particle codes.
    int id = event[i].id();
    if (id == 0 || !pState.particleData.isParticle(id)) {
      ostringstream errCode;
      errCode << ", i = " << i << ", id = " << id;
      pState.info.errorMsg("Error in Pythia::check: "
        "unknown particle code", errCode.str());
      physical = false;
      iErrId.push_back(i);

    // Check that colour assignments are the expected ones.
    } else {
      int colType = event[i].colType();
      int col     = event[i].col();
      int acol    = event[i].acol();
      if ( (colType ==  0 && (col  > 0 || acol  > 0))
        || (colType ==  1 && (col <= 0 || acol  > 0))
        || (colType == -1 && (col  > 0 || acol <= 0))
        || (colType ==  2 && (col <= 0 || acol <= 0)) ) {
        ostringstream errCode;
        errCode << ", i = " << i << ", id = " << id << " cols = " << col
                << ' ' << acol;
        pState.info.errorMsg("Error in Pythia::check: "
          "incorrect colours", errCode.str());
        physical = false;
        iErrCol.push_back(i);
      }
    }

    Benchmark_stop(Pythia0check_find_unrecognized);
    Benchmark_start(Pythia0check_check_momentum);

    // Look for particles with mismatched or not-a-number energy/momentum/mass.
    if (abs(event[i].px()) >= 0. && abs(event[i].py()) >= 0.
      && abs(event[i].pz()) >= 0.  && abs(event[i].e()) >= 0.
      && abs(event[i].m()) >= 0.) {
      double errMass = abs(event[i].mCalc() - event[i].m())
        / max( 1.0, event[i].e());
      if (errMass > mTolErr) {
        pState.info.errorMsg("Error in Pythia::check: "
          "unmatched particle energy/momentum/mass");
        physical = false;
        iErrEpm.push_back(i);
      } else if (errMass > mTolWarn) {
        pState.info.errorMsg("Warning in Pythia::check: "
          "not quite matched particle energy/momentum/mass");
      }
    } else {
      pState.info.errorMsg("Error in Pythia::check: "
        "not-a-number energy/momentum/mass");
      physical = false;
      iErrNan.push_back(i);
    }

    Benchmark_stop(Pythia0check_check_momentum);
    Benchmark_start(Pythia0check_check_nans);

    // Look for particles with not-a-number vertex/lifetime.
    if (abs(event[i].xProd()) >= 0. && abs(event[i].yProd()) >= 0.
      && abs(event[i].zProd()) >= 0.  && abs(event[i].tProd()) >= 0.
      && abs(event[i].tau()) >= 0.) ;
    else {
      pState.info.errorMsg("Error in Pythia::check: "
        "not-a-number vertex/lifetime");
      physical     = false;
      listVertices = true;
      iErrNanVtx.push_back(i);
    }

    Benchmark_stop(Pythia0check_check_nans);
    Benchmark_start(Pythia0check_check_charge_total);

    // Add final-state four-momentum and charge.
    if (event[i].isFinal()) {
      pSum      += event[i].p();
      chargeSum += event[i].charge();
    }

  // End of particle loop.
  }

  Benchmark_loopStop(Pythia0check_loopAllParticles);
  Benchmark_start(Pythia0check_check_momentum_charge_conservation);

  // Check energy-momentum/charge conservation.
  double epDev = abs(pSum.e()) + abs(pSum.px()) + abs(pSum.py())
    + abs(pSum.pz());
  if (epDev > epTolErr * eLab) {
    pState.info.errorMsg("Error in Pythia::check: energy-momentum not conserved");
    physical = false;
  } else if (epDev > epTolWarn * eLab) {
    pState.info.errorMsg("Warning in Pythia::check: "
      "energy-momentum not quite conserved");
  }
  if (abs(chargeSum) > 0.1) {
    pState.info.errorMsg("Error in Pythia::check: charge not conserved");
    physical = false;
  }

  Benchmark_stop(Pythia0check_check_momentum_charge_conservation);
  Benchmark_start(Pythia0check_check_records);

  // Check that beams and event records agree on incoming partons.
  // Only meaningful for resolved beams.
  if (pState.info.isResolved() && !pState.info.hasUnresolvedBeams())
  for (int iSys = 0; iSys < beamA.sizeInit(); ++iSys) {
    int eventANw  = partonSystems.getInA(iSys);
    int eventBNw  = partonSystems.getInB(iSys);
    int beamANw   = (*&beamA)[iSys].iPos();
    int beamBNw   = (*&beamB)[iSys].iPos();
    if (eventANw != beamANw || eventBNw != beamBNw) {
      pState.info.errorMsg("Error in Pythia::check: "
        "event and beams records disagree");
      physical    = false;
      listSystems = true;
      listBeams   = true;
    }
  }

  Benchmark_stop(Pythia0check_check_records);
  Benchmark_start(Pythia0check_check_mothers_and_daughters_match);

  // Check that mother and daughter information match for each particle.
  vector<int> noMot;
  vector<int> noDau;
  vector< pair<int,int> > noMotDau;
  if (checkHistory) {

    // Loop through the event and check that there are beam particles.
    bool hasBeams = false;
    for (int i = 0; i < event.size(); ++i) {
      int status = event[i].status();
      if (abs(status) == 12) hasBeams = true;

      // Check that mother and daughter lists not empty where not expected to.
      vector<int> mList = event[i].motherList();
      vector<int> dList = event[i].daughterList();
      if (mList.size() == 0 && abs(status) != 11 && abs(status) != 12)
        noMot.push_back(i);
      if (dList.size() == 0 && status < 0 && status != -11)
        noDau.push_back(i);

      // Check that the particle appears in the daughters list of each mother.
      for (int j = 0; j < int(mList.size()); ++j) {
        if ( event[mList[j]].daughter1() <= i
          && event[mList[j]].daughter2() >= i ) continue;
        vector<int> dmList = event[mList[j]].daughterList();
        bool foundMatch = false;
        for (int k = 0; k < int(dmList.size()); ++k)
        if (dmList[k] == i) {
          foundMatch = true;
          break;
        }
        if (!hasBeams && mList.size() == 1 && mList[0] == 0) foundMatch = true;
        if (!foundMatch) {
          bool oldPair = false;
          for (int k = 0; k < int(noMotDau.size()); ++k)
          if (noMotDau[k].first == mList[j] && noMotDau[k].second == i) {
            oldPair = true;
            break;
          }
          if (!oldPair) noMotDau.push_back( make_pair( mList[j], i) );
        }
      }

      // Check that the particle appears in the mothers list of each daughter.
      for (int j = 0; j < int(dList.size()); ++j) {
        if ( event[dList[j]].statusAbs() > 80
          && event[dList[j]].statusAbs() < 90
          && event[dList[j]].mother1() <= i
          && event[dList[j]].mother2() >= i) continue;
        vector<int> mdList = event[dList[j]].motherList();
        bool foundMatch = false;
        for (int k = 0; k < int(mdList.size()); ++k)
        if (mdList[k] == i) {
          foundMatch = true;
          break;
        }
        if (!foundMatch) {
          bool oldPair = false;
          for (int k = 0; k < int(noMotDau.size()); ++k)
          if (noMotDau[k].first == i && noMotDau[k].second == dList[j]) {
            oldPair = true;
            break;
          }
          if (!oldPair) noMotDau.push_back( make_pair( i, dList[j]) );
        }
      }
    }

    // Warn if any errors were found.
    if (noMot.size() > 0 || noDau.size() > 0 || noMotDau.size() > 0) {
      pState.info.errorMsg("Error in Pythia::check: "
        "mismatch in daughter and mother lists");
      physical    = false;
      listHistory = true;
    }
  }

  // Done for sensible events.
  if (physical) return true;

  // Print (the first few) flawed events: local info.
  if (nErrEvent < nErrList) {
    os << "\n PYTHIA erroneous event info: \n";
    if (iErrId.size() > 0) {
      os << " unknown particle codes in lines ";
      for (int i = 0; i < int(iErrId.size()); ++i)
        os << iErrId[i] << ' ';
      os << "\n";
    }
    if (iErrCol.size() > 0) {
      os << " incorrect colour assignments in lines ";
      for (int i = 0; i < int(iErrCol.size()); ++i)
        os << iErrCol[i] << ' ';
      os << "\n";
    }
    if (iErrEpm.size() > 0) {
      os << " mismatch between energy/momentum/mass in lines ";
      for (int i = 0; i < int(iErrEpm.size()); ++i)
        os << iErrEpm[i] << ' ';
      os << "\n";
    }
    if (iErrNan.size() > 0) {
      os << " not-a-number energy/momentum/mass in lines ";
      for (int i = 0; i < int(iErrNan.size()); ++i)
        os << iErrNan[i] << ' ';
      os << "\n";
    }
    if (iErrNanVtx.size() > 0) {
      os << " not-a-number vertex/lifetime in lines ";
      for (int i = 0; i < int(iErrNanVtx.size()); ++i)
        os << iErrNanVtx[i] << ' ';
      os << "\n";
    }
    if (epDev > epTolErr * eLab) os << scientific << setprecision(3)
      << " total energy-momentum non-conservation = " << epDev << "\n";
    if (abs(chargeSum) > 0.1) os << fixed << setprecision(2)
      << " total charge non-conservation = " << chargeSum << "\n";
    if (noMot.size() > 0) {
      os << " missing mothers for particles ";
      for (int i = 0; i < int(noMot.size()); ++i) os << noMot[i] << ' ';
      os << "\n";
    }
    if (noDau.size() > 0) {
      os << " missing daughters for particles ";
      for (int i = 0; i < int(noDau.size()); ++i) os << noDau[i] << ' ';
      os << "\n";
    }
    if (noMotDau.size() > 0) {
      os << " inconsistent history for (mother,daughter) pairs ";
      for (int i = 0; i < int(noMotDau.size()); ++i)
        os << '(' << noMotDau[i].first << ',' << noMotDau[i].second << ") ";
      os << "\n";
    }

    // Print (the first few) flawed events: standard listings.
    pState.info.list();
    event.list(listVertices, listHistory);
    if (listSystems) partonSystems.list();
    if (listBeams) beamA.list();
    if (listBeams) beamB.list();
  }

  // Update error counter. Done also for flawed event.
  ++nErrEvent;
  return false;

}

//--------------------------------------------------------------------------

// Routine to set up a PDF pointer.

// done
PDF* Pythia::getPDFPtr(int idIn, int sequence, stringref beam) {

  // Temporary pointer to be returned.
  PDF* tempPDFPtr = 0;

  // One option is to treat a Pomeron like a pi0.
  if (idIn == 990 && pState.settings.get(Mode::PDF_PomSet) == 2) idIn = 111;

  // Proton beam, normal or hard choice. Also used for neutron.
  if (abs(idIn) == 2212 || abs(idIn) == 2112) {
    string pSet;
    if (beam == "")
    {
      if (sequence == 1) pSet = pState.settings.get(Word::PDF_pSet);
      if (sequence != 1) pSet = pState.settings.get(Word::PDF_pHardSet);
    }
    else if (beam == "B")
    {
      if (sequence == 1) pSet = pState.settings.get(Word::PDF_pSetB);
      if (sequence != 1) pSet = pState.settings.get(Word::PDF_pHardSetB);
    }
    if (pSet == "void" && sequence != 1 && beam == "B")
      pSet = pState.settings.get(Word::PDF_pHardSet);
    if (pSet == "void") pSet = pState.settings.get(Word::PDF_pSet);
    istringstream pSetStream(pSet);
    int pSetInt(0);
    pSetStream >> pSetInt;

    // Use sets from LHAPDF.
    if (pSetInt == 0)
      tempPDFPtr = new LHAPDF(idIn, pSet, &pState.info);

    // Use internal sets.
    else if (pSetInt == 1) tempPDFPtr = new GRV94L(idIn);
    else if (pSetInt == 2) tempPDFPtr = new CTEQ5L(idIn);
    else if (pSetInt <= 6)
      tempPDFPtr = new MSTWpdf(idIn, pSetInt - 2, xmlPath, &pState.info);
    else if (pSetInt <= 12)
      tempPDFPtr = new CTEQ6pdf(idIn, pSetInt - 6, xmlPath, &pState.info);
    else if (pSetInt <= 16)
      tempPDFPtr = new NNPDF(idIn, pSetInt - 12, xmlPath, &pState.info);
    else tempPDFPtr = 0;
  }

  // Pion beam (or, in one option, Pomeron beam).
  else if (abs(idIn) == 211 || idIn == 111) {
    string pSet;
    if (beam == "") pSet = pState.settings.get(Word::PDF_piSet);
    if (beam == "B") pSet = pState.settings.get(Word::PDF_piSetB);
    istringstream pSetStream(pSet);
    int pSetInt(0);
    pSetStream >> pSetInt;

    // Use sets from LHAPDF.
    if (pSetInt == 0)
      tempPDFPtr = new LHAPDF(idIn, pSet, &pState.info);

    // Use internal set.
    else if (pSetInt == 1) tempPDFPtr = new GRVpiL(idIn);
    else tempPDFPtr = 0;
  }

  // Pomeron beam, if not treated like a pi0 beam.
  else if (idIn == 990) {
    int    pomSet  = pState.settings.get(Mode::PDF_PomSet);
    double rescale = pState.settings.get(Param::PDF_PomRescale);

    // A generic Q2-independent parametrization.
    if (pomSet == 1) {
      double gluonA      = pState.settings.get(Param::PDF_PomGluonA);
      double gluonB      = pState.settings.get(Param::PDF_PomGluonB);
      double quarkA      = pState.settings.get(Param::PDF_PomQuarkA);
      double quarkB      = pState.settings.get(Param::PDF_PomQuarkB);
      double quarkFrac   = pState.settings.get(Param::PDF_PomQuarkFrac);
      double strangeSupp = pState.settings.get(Param::PDF_PomStrangeSupp);
      tempPDFPtr = new PomFix( 990, gluonA, gluonB, quarkA, quarkB,
        quarkFrac, strangeSupp);
    }

    // The H1 Q2-dependent parametrizations. Initialization requires files.
    else if (pomSet == 3 || pomSet == 4)
      tempPDFPtr = new PomH1FitAB( 990, pomSet - 2, rescale, xmlPath, &pState.info);
    else if (pomSet == 5)
      tempPDFPtr = new PomH1Jets( 990, rescale, xmlPath, &pState.info);
    else if (pomSet == 6)
      tempPDFPtr = new PomH1FitAB( 990, 3, rescale, xmlPath, &pState.info);
  }

  // Lepton beam: neutrino, resolved charged lepton or unresolved ditto.
  else if (abs(idIn) > 10 && abs(idIn) < 17) {
    if (abs(idIn)%2 == 0) tempPDFPtr = new NeutrinoPoint(idIn);
    else if (pState.settings.get(Flag::PDF_lepton)) tempPDFPtr = new Lepton(idIn);
    else tempPDFPtr = new LeptonPoint(idIn);
  }

  // Optionally allow extrapolation beyond x and Q2 limits.
  if (tempPDFPtr)
    tempPDFPtr->setExtrapolate( pState.settings.get(Flag::PDF_extrapolate) );

  // Done.
  return tempPDFPtr;
}

//==========================================================================

} // end namespace Pythia8
