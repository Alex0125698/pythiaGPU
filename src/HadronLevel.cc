// HadronLevel.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the HadronLevel class.

#include "Pythia8/HadronLevel.h"

namespace Pythia8 {

//==========================================================================

// The HadronLevel class.

//--------------------------------------------------------------------------

// Find settings. Initialize HadronLevel classes as required.

bool HadronLevel::init(Info* infoPtrIn, Settings& settings,
  ParticleData* particleDataPtrIn, Rndm* rndmPtrIn,
  Couplings* couplingsPtrIn, TimeShower* timesDecPtr,
  RHadrons* rHadronsPtrIn, DecayHandler* decayHandlePtr,
  const vector<int>& handledParticles, UserHooks* userHooksPtrIn) {

  // Save pointers.
  infoPtr         = infoPtrIn;
  particleDataPtr = particleDataPtrIn;
  rndmPtr         = rndmPtrIn;
  couplingsPtr    = couplingsPtrIn;
  rHadronsPtr     = rHadronsPtrIn;
  userHooksPtr    = userHooksPtrIn;

  // Main flags.
  doHadronize     = settings.flag("HadronLevel:Hadronize");
  doDecay         = settings.flag("HadronLevel:Decay");
  doBoseEinstein  = settings.flag("HadronLevel:BoseEinstein");

  // Boundary mass between string and ministring handling.
  mStringMin      = settings.parm("HadronLevel:mStringMin");

  // For junction processing.
  eNormJunction   = settings.parm("StringFragmentation:eNormJunction");

  // Allow R-hadron formation.
  allowRH         = settings.flag("RHadrons:allow");

  // Particles that should decay or not before Bose-Einstein stage.
  widthSepBE      = settings.parm("BoseEinstein:widthSep");

  // Hadron scattering --rjc
  doHadronScatter = settings.flag("HadronScatter:scatter");
  hsAfterDecay    = settings.flag("HadronScatter:afterDecay");

  // Initialize auxiliary fragmentation classes.
  flavSel.init(settings, rndmPtr);
  pTSel.init(settings, *particleDataPtr, rndmPtr);
  zSel.init(settings, *particleDataPtr, rndmPtr);

  // Initialize auxiliary administrative class.
  colConfig.init(infoPtr, settings, &flavSel);

  // Initialize string and ministring fragmentation.
  stringFrag.init(infoPtr, settings, particleDataPtr, rndmPtr,
    &flavSel, &pTSel, &zSel, userHooksPtr);
  ministringFrag.init(infoPtr, settings, particleDataPtr, rndmPtr,
    &flavSel, &pTSel, &zSel);

  // Initialize particle decays.
  decays.init(infoPtr, settings, particleDataPtr, rndmPtr, couplingsPtr,
    timesDecPtr, &flavSel, decayHandlePtr, handledParticles);

  // Initialize BoseEinstein.
  boseEinstein.init(infoPtr, settings, *particleDataPtr);

  // Initialize HadronScatter --rjc
  if (doHadronScatter)
    hadronScatter.init(infoPtr, settings, rndmPtr, particleDataPtr);

  // Initialize Hidden-Valley fragmentation, if necessary.
  useHiddenValley = hiddenvalleyFrag.init(infoPtr, settings,
    particleDataPtr, rndmPtr);

  // Send flavour and z selection pointers to R-hadron machinery.
  rHadronsPtr->fragPtrs( &flavSel, &zSel);

  // Initialize the colour tracing class.
  colTrace.init(infoPtr);

  // Initialize the junction splitting class.
  junctionSplitting.init(infoPtr, settings, rndmPtr, particleDataPtr);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Hadronize and decay the next parton-level.

bool HadronLevel::next( Event& event) {

  Benchmark_start(HadronLevel0next);

  // Store current event size to mark Parton Level content.
  event.savePartonLevelSize();

  Benchmark_start(HadronLevel0next_hiddenValleyFrag);

  // Do Hidden-Valley fragmentation, if necessary.
  if (useHiddenValley) hiddenvalleyFrag.fragment(event);
  
  Benchmark_stop(HadronLevel0next_hiddenValleyFrag);
  Benchmark_start(HadronLevel0next_oniaDecay);

  // Colour-octet onia states must be decayed to singlet + gluon.
  if (!decayOctetOnia(event)) return false;
  
  Benchmark_stop(HadronLevel0next_oniaDecay);
  Benchmark_start(HadronLevel0next_cleanJunctions);

  // remove junction structures.
  if (!junctionSplitting.checkColours(event)) {
    infoPtr->errorMsg("Error in HadronLevel::next: "
        "failed colour/junction check");
    return false;
  }

  Benchmark_stop(HadronLevel0next_cleanJunctions);
  Benchmark_loopStart(HadronLevel0next_loopMultiHadr);
  bool failed = false;

  // Possibility of hadronization inside decay, but then no BE second time.
  // Hadron scattering, first pass only --rjc
  bool moreToDo, firstPass = true;
  bool doBoseEinsteinNow = doBoseEinstein;
  do 
  {
    moreToDo = false;

    Benchmark_loopCount(HadronLevel0next_loopMultiHadr);

    // First part: string fragmentation.
    if (doHadronize) 
    {
      Benchmark_start(HadronLevel0next_fragFindSinglets);

      // Find the complete colour singlet configuration of the event.
      if (!findSinglets( event)) { failed=true; break; };

      Benchmark_stop(HadronLevel0next_fragFindSinglets);
      Benchmark_start(HadronLevel0next_fragRHadrons);

      // Fragment off R-hadrons, if necessary.
      if (allowRH && !rHadronsPtr->produce( colConfig, event))
        { failed=true; break; };

      Benchmark_stop(HadronLevel0next_fragRHadrons);
      Benchmark_loopStart(HadronLevel0next_loopColorSinglets);

      // Process all colour singlet (sub)systems.
      for (int iSub = 0; iSub < colConfig.size(); ++iSub) {

        Benchmark_loopCount(HadronLevel0next_loopColorSinglets);
        Benchmark_start(HadronLevel0next_fragCollectPartons);

        // Collect sequentially all partons in a colour singlet subsystem.
        colConfig.collect(iSub, event);

        Benchmark_stop(HadronLevel0next_fragCollectPartons);
        Benchmark_placeholder(HadronLevel0next_fragStringFrag);
        Benchmark_placeholder(HadronLevel0next_fragMiniStringFrag);

        // String fragmentation of each colour singlet (sub)system.
        if ( colConfig[iSub].massExcess > mStringMin ) 
        {
          Benchmark_start(HadronLevel0next_fragStringFrag);
          if (!stringFrag.fragment( iSub, colConfig, event)) { failed=true; break; };
          Benchmark_stop(HadronLevel0next_fragStringFrag);

        }


        // Low-mass string treated separately. Tell if diffractive system.
        else 
        {
          Benchmark_placeholder(HadronLevel0next_fixBug);
          Benchmark_start(HadronLevel0next_fragMiniStringFrag);
          bool isDiff = infoPtr->isDiffractiveA() || infoPtr->isDiffractiveB();
          if (!ministringFrag.fragment( iSub, colConfig, event, isDiff))
            { failed=true; break; };
          Benchmark_stop(HadronLevel0next_fragMiniStringFrag);
        }
      }

      Benchmark_loopStop(HadronLevel0next_loopColorSinglets);
    }

    Benchmark_start(HadronLevel0next_fragHadronScatter);

    // Hadron scattering --rjc
    if (doHadronScatter && !hsAfterDecay && firstPass)
      hadronScatter.scatter(event);
    
    Benchmark_stop(HadronLevel0next_fragHadronScatter);
    Benchmark_start(HadronLevel0next_sequentialDecayShort);

    // Second part: sequential decays of short-lived particles (incl. K0).
    if (doDecay) {

      // Loop through all entries to find those that should decay.
      int iDec = 0;
      do 
      {
        Particle& decayer = event[iDec];
        if ( decayer.isFinal() && decayer.canDecay() && decayer.mayDecay()
             && (decayer.mWidth() > widthSepBE || decayer.idAbs() == 311) ) 
        {
          decays.decay( iDec, event);
          if (decays.moreToDo()) moreToDo = true;
        }
      } while (++iDec < event.size());
    }

    Benchmark_stop(HadronLevel0next_sequentialDecayShort);
    Benchmark_start(HadronLevel0next_hadronScatter2);

    // Hadron scattering --rjc
    if (doHadronScatter && hsAfterDecay && firstPass)
      hadronScatter.scatter(event);

    Benchmark_stop(HadronLevel0next_hadronScatter2);
    Benchmark_start(HadronLevel0next_boseEinsteinEffects);

    // Third part: include Bose-Einstein effects among current particles.
    if (doBoseEinsteinNow) {
      if (!boseEinstein.shiftEvent(event)) { failed=true; break; };
      doBoseEinsteinNow = false;
    }

    Benchmark_stop(HadronLevel0next_boseEinsteinEffects);
    Benchmark_start(HadronLevel0next_sequentialDecayAll);

    // Fourth part: sequential decays also of long-lived particles.
    if (doDecay) {

      // Loop through all entries to find those that should decay.
      int iDec = 0;
      do {
        Particle& decayer = event[iDec];
        if ( decayer.isFinal() && decayer.canDecay() && decayer.mayDecay() ) {
          decays.decay( iDec, event);
          if (decays.moreToDo()) moreToDo = true;
        }
      } while (++iDec < event.size());
    }

    Benchmark_stop(HadronLevel0next_sequentialDecayAll);

  // Normally done first time around, but sometimes not (e.g. Upsilon).
  } while (moreToDo);

  Benchmark_loopStop(HadronLevel0next_loopMultiHadr);

  // Done.
  return !failed;

}

//--------------------------------------------------------------------------

// Allow more decays if on/off switches changed.
// Note: does not do sequential hadronization, e.g. for Upsilon.

bool HadronLevel::moreDecays( Event& event) {

  // Colour-octet onia states must be decayed to singlet + gluon.
  if (!decayOctetOnia(event)) return false;

  // Loop through all entries to find those that should decay.
  int iDec = 0;
  do {
    if ( event[iDec].isFinal() && event[iDec].canDecay()
      && event[iDec].mayDecay() ) decays.decay( iDec, event);
  } while (++iDec < event.size());

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Decay colour-octet onium states.

bool HadronLevel::decayOctetOnia(Event& event) {

  Benchmark_start(HadronLevel0decayOctetOnia);
  Benchmark_loopStart(HadronLevel0decayOctetOnia);

  // Loop over particles and decay any onia encountered.
  for (int iDec = 0; iDec < event.size(); ++iDec)
  {
    Benchmark_loopCount(HadronLevel0decayOctetOnia);

    Benchmark_start(HadronLevel0decayOctetOnia_decay);
    
    if (event[iDec].isFinal()
      && particleDataPtr->isOctetHadron(event[iDec].id())) 
    {
      if (!decays.decay( iDec, event)) return false;

      // Set colour flow by hand: gluon inherits octet-onium state.
      int iGlu = event.size() - 1;
      event[iGlu].cols( event[iDec].col(), event[iDec].acol() );
    }
  }

  Benchmark_loopStop(HadronLevel0decayOctetOnia);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Trace colour flow in the event to form colour singlet subsystems.

bool HadronLevel::findSinglets(Event& event) {

  Benchmark_start(HadronLevel0findSinglets);

  // Clear up storage.
  colConfig.clear();

  Benchmark_start(HadronLevel0findSinglets_setupColList);
  Benchmark_placeholder(HadronLevel0findSinglets_arrange);
  Benchmark_placeholder(HadronLevel0findSinglets_openStrings);
  Benchmark_placeholder(HadronLevel0findSinglets_closedStrings);

  // Find a list of final partons and of all colour ends and gluons.
  if (colTrace.setupColList(event)) return true;

  Benchmark_stop(HadronLevel0findSinglets_setupColList);
  Benchmark_start(HadronLevel0findSinglets_arrange);

  // Begin arrange the partons into separate colour singlets.

  // Junctions: loop over them, and identify kind.
  for (int iJun = 0; iJun < event.sizeJunction(); ++iJun)
  if (event.remainsJunction(iJun)) {
    event.remainsJunction(iJun, false);
    int kindJun = event.kindJunction(iJun);
    iParton.resize(0);

    // Loop over junction legs.
    for (int iCol = 0; iCol < 3; ++iCol) {
      int indxCol = event.colJunction(iJun, iCol);
      iParton.push_back( -(10 + 10 * iJun + iCol) );
      // Junctions: find color ends.
      if (kindJun % 2 == 1 && !colTrace.traceFromAcol(indxCol, event, iJun,
        iCol, iParton)) return false;
      // Antijunctions: find anticolor ends.
      if (kindJun % 2 == 0 && !colTrace.traceFromCol(indxCol, event, iJun,
        iCol, iParton)) return false;
    }

    // A junction may be eliminated by insert if two quarks are nearby.
    int nJunOld = event.sizeJunction();
    if (!colConfig.insert(iParton, event)) return false;
    if (event.sizeJunction() < nJunOld) --iJun;
  }

  Benchmark_stop(HadronLevel0findSinglets_arrange);
  Benchmark_start(HadronLevel0findSinglets_openStrings);

  // Open strings: pick up each colour end and trace to its anticolor end.
  while (!colTrace.colFinished()) {
    iParton.resize(0);
    if (!colTrace.traceFromCol( -1, event, -1, -1, iParton)) return false;

    // Store found open string system. Analyze its properties.
    if (!colConfig.insert(iParton, event)) return false;
  }

  Benchmark_stop(HadronLevel0findSinglets_openStrings);
  Benchmark_start(HadronLevel0findSinglets_closedStrings);

  // Closed strings : begin at any gluon and trace until back at it.
  while (!colTrace.finished()) {
    iParton.resize(0);
    if (!colTrace.traceInLoop(event, iParton)) return false;

    // Store found closed string system. Analyze its properties.
    if (!colConfig.insert(iParton, event)) return false;
  }

  // Done.
  return true;

}

//==========================================================================

} // end namespace Pythia8
