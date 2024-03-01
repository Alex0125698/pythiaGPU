// ProcessContainer.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// ProcessContainer and SetupContainers classes.

#include "Pythia8/ProcessContainer.h"

// Internal headers for special processes.
#include "Pythia8/SigmaCompositeness.h"
#include "Pythia8/SigmaEW.h"
#include "Pythia8/SigmaExtraDim.h"
#include "Pythia8/SigmaGeneric.h"
#include "Pythia8/SigmaHiggs.h"
#include "Pythia8/SigmaLeftRightSym.h"
#include "Pythia8/SigmaLeptoquark.h"
#include "Pythia8/SigmaNewGaugeBosons.h"
#include "Pythia8/SigmaQCD.h"
#include "Pythia8/SigmaSUSY.h"

namespace Pythia8 {

//==========================================================================

// ProcessContainer class.
// Information allowing the generation of a specific process.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Number of event tries to check maximization finding reliability.
const int ProcessContainer::N12SAMPLE = 100;

// Ditto, but increased for 2 -> 3 processes.
const int ProcessContainer::N3SAMPLE  = 1000;

//--------------------------------------------------------------------------

// Initialize phase space and counters.
// Argument isFirst distinguishes two hard processes in same event.

bool ProcessContainer::init(bool isFirst, Info* infoPtrIn,
  Settings& settings, ParticleData* particleDataPtrIn, Rndm* rndmPtrIn,
  BeamParticle* beamAPtr, BeamParticle* beamBPtr, Couplings* couplingsPtr,
  SigmaTotal* sigmaTotPtr, ResonanceDecays* resDecaysPtrIn,
  SLHAinterface* slhaInterfacePtr, UserHooks* userHooksPtrIn, PythiaState* pState) {

  Benchmark_start(ProcessContainer0init);
  Benchmark_start(ProcessContainer0init_setup); // trivial

  // Extract info about current process from SigmaProcess object.
  isLHA       = sigmaProcessPtr->isLHA;
  isNonDiff   = sigmaProcessPtr->isNonDiff;
  isResolved  = sigmaProcessPtr->isResolved;
  isDiffA     = sigmaProcessPtr->isDiffA;
  isDiffB     = sigmaProcessPtr->isDiffB;
  isDiffC     = sigmaProcessPtr->isDiffC;
  isQCD3body  = sigmaProcessPtr->isQCD3body;
  int nFin    = sigmaProcessPtr->nFinal;
  lhaStrat    = (isLHA) ? lhaUpPtr->strategy() : 0;
  lhaStratAbs = abs(lhaStrat);
  allowNegSig = sigmaProcessPtr->allowNegativeSigma;

  useStrictLHEFscales = settings.get(Flag::Beams_strictLHEFscale);

  // Flag for maximum violation handling.
  increaseMaximum = settings.get(Flag::PhaseSpace_increaseMaximum);

  // Note: the phase space pointer setup here

  // Pick and create phase space generator. Send pointers where required.
  if (phaseSpacePtr != 0) ;
  else if (isLHA)       phaseSpacePtr = new PhaseSpaceLHA();
  else if (isNonDiff)   phaseSpacePtr = new PhaseSpace2to2nondiffractive();
  else if (!isResolved && !isDiffA  && !isDiffB  && !isDiffC )
                        phaseSpacePtr = new PhaseSpace2to2elastic();
  else if (!isResolved && !isDiffA  && !isDiffB && isDiffC)
                        phaseSpacePtr = new PhaseSpace2to3diffractive();
  else if (!isResolved) phaseSpacePtr = new PhaseSpace2to2diffractive(
                                        isDiffA, isDiffB);
  else if (nFin == 1)   phaseSpacePtr = new PhaseSpace2to1tauy();
  else if (nFin == 2)   phaseSpacePtr = new PhaseSpace2to2tauyz();
  else if (isQCD3body == true)  phaseSpacePtr = new PhaseSpace2to3yyycyl();
  else                  phaseSpacePtr = new PhaseSpace2to3tauycyl();

  // Store pointers and perform simple initialization.
  infoPtr         = infoPtrIn;
  particleDataPtr = particleDataPtrIn;
  rndmPtr         = rndmPtrIn;
  resDecaysPtr    = resDecaysPtrIn;
  userHooksPtr    = userHooksPtrIn;
  canVetoResDecay = (userHooksPtr != 0)
                  ? userHooksPtr->canVetoResonanceDecays() : false;
  if (isLHA) {
    // sigmaProcessPtr->setLHAPtr(lhaUpPtr);
    phaseSpacePtr->setLHAPtr(lhaUpPtr);
  }

  Benchmark_stop(ProcessContainer0init_setup);
  Benchmark_start(ProcessContainer0init_sigmaProcess); // todo

  // SigmaProcesses are initialized here
  // I think this is where the strings get created

  pState->beamA = beamAPtr;
  pState->beamB = beamBPtr;
  sigmaProcessPtr->init(pState);
  
  Benchmark_stop(ProcessContainer0init_sigmaProcess);
  Benchmark_start(ProcessContainer0init_phaseSpace); // trivial
  
  phaseSpacePtr->init( isFirst, sigmaProcessPtr, infoPtr, &settings,
    particleDataPtr, rndmPtr, beamAPtr,  beamBPtr, couplingsPtr, sigmaTotPtr,
    userHooksPtr);
  
  Benchmark_stop(ProcessContainer0init_phaseSpace);
  Benchmark_start(ProcessContainer0init_initProcess); // todo

  // Reset cross section statistics.
  nTry      = 0;
  nSel      = 0;
  nAcc      = 0;
  nTryStat  = 0;
  sigmaMx   = 0.;
  sigmaSum  = 0.;
  sigma2Sum = 0.;
  sigmaNeg  = 0.;
  sigmaAvg  = 0.;
  sigmaFin  = 0.;
  deltaFin  = 0.;
  wtAccSum  = 0.;

  // Initialize process and allowed incoming partons.
  sigmaProcessPtr->initProc(); // TODO: look at the actual function (this is virtual)

  Benchmark_stop(ProcessContainer0init_initProcess);
  Benchmark_start(ProcessContainer0init_setupIncomingPartonFlux); // trivial

  // TODO: what is it for?

  if (!sigmaProcessPtr->initFlux()) return false;
  
  Benchmark_stop(ProcessContainer0init_setupIncomingPartonFlux);
  Benchmark_start(ProcessContainer0init_setupPhaseSpaceSampling); // todo

  // Find maximum of differential cross section * phasespace.
  bool physical       = phaseSpacePtr->setupSampling(); // this is virtual
  
  
  sigmaMx             = phaseSpacePtr->sigmaMax();
  double sigmaHalfWay = sigmaMx;

  // Separate signed maximum needed for LHA with negative weight.
  sigmaSgn            = phaseSpacePtr->sigmaSumSigned();

  Benchmark_stop(ProcessContainer0init_setupPhaseSpaceSampling);

  // Check maximum by a few events, and extrapolate a further increase.
  if (physical & !isLHA) 
  {
    Benchmark_loopStart(ProcessContainer0init_loopOverSamples);

    int nSample = (nFin < 3) ? N12SAMPLE : N3SAMPLE;

    for (int iSample = 0; iSample < nSample; ++iSample) 
    {
      Benchmark_loopCount(ProcessContainer0init_loopOverSamples);
      Benchmark_loopStart(ProcessContainer0init_trialKin);

      bool test = false;
      while (!test)
      {
        Benchmark_loopCount(ProcessContainer0init_trialKin);
        Benchmark_start(ProcessContainer0init_phasetrialKin); // todo
        test = phaseSpacePtr->trialKin(false); // this is virtual
      }

      Benchmark_loopStop(ProcessContainer0init_trialKin);

      if (iSample == nSample/2) sigmaHalfWay = phaseSpacePtr->sigmaMax();
    }

    Benchmark_loopStop(ProcessContainer0init_loopOverSamples);

    double sigmaFullWay = phaseSpacePtr->sigmaMax();
    sigmaMx = (sigmaHalfWay > 0.) ? pow2(sigmaFullWay) / sigmaHalfWay
                                  : sigmaFullWay;

    phaseSpacePtr->setSigmaMax(sigmaMx);
  }

  Benchmark_start(ProcessContainer0init_moreSettings);

  // Allow Pythia to overwrite incoming beams or parts of Les Houches input.
  idRenameBeams = settings.get(Mode::LesHouches_idRenameBeams);
  setLifetime   = settings.get(Mode::LesHouches_setLifetime);
  setLeptonMass = settings.get(Mode::LesHouches_setLeptonMass);
  mRecalculate  = settings.get(Param::LesHouches_mRecalculate);
  matchInOut    = settings.get(Flag::LesHouches_matchInOut);
  idLep[0] = 11; mLep[0] = particleDataPtr->m0(11);
  idLep[1] = 13; mLep[1] = particleDataPtr->m0(13);
  idLep[2] = 15; mLep[2] = particleDataPtr->m0(15);

  // Done.
  return physical;
}

//--------------------------------------------------------------------------

// Generate a trial event; selected or not.

bool ProcessContainer::trialProcess() {

  Benchmark_start(trialProcess);
  Benchmark_placeholder(trialProcess_dummyA)
  Benchmark_loopStart(trialProcess_loopTrialProcess);

  bool quitLoop = false;

  // Loop over tries only occurs for Les Houches strategy = +-2.
  for (int iTry = 0;  ; ++iTry) {

    Benchmark_loopCount(trialProcess_loopTrialProcess);

    // Generate a trial phase space point, if meaningful.
    if (sigmaMx == 0.) { quitLoop = true; break;}
    infoPtr->setEndOfFile(false);
    bool repeatSame = (iTry > 0);
    
    Benchmark_start(trialProcess_trialKin);

    bool physical = phaseSpacePtr->trialKin(true, repeatSame);

    Benchmark_stop(trialProcess_trialKin);
    Benchmark_start(trialProcess_rest);

    // Note if at end of Les Houches file, else do statistics.
    if (isLHA && !physical) infoPtr->setEndOfFile(true);
    else {
      ++nTry;
      // Statistics for separate Les Houches process codes. Insert new codes.
      if (isLHA) {
        int codeLHANow = lhaUpPtr->idProcess();
        int iFill = -1;
        for (int i = 0; i < int(codeLHA.size()); ++i)
          if (codeLHANow == codeLHA[i]) iFill = i;
        if (iFill >= 0) {
          ++nTryLHA[iFill];
        } else {
          codeLHA.push_back(codeLHANow);
          nTryLHA.push_back(1);
          nSelLHA.push_back(0);
          nAccLHA.push_back(0);
          for (int i = int(codeLHA.size()) - 1; i > 0; --i) {
            if (codeLHA[i] < codeLHA[i - 1]) {
              swap(codeLHA[i], codeLHA[i - 1]);
              swap(nTryLHA[i], nTryLHA[i - 1]);
              swap(nSelLHA[i], nSelLHA[i - 1]);
              swap(nAccLHA[i], nAccLHA[i - 1]);
            }
            else break;
          }
        }
      }
    }

    // Possibly fail, else cross section.
    if (!physical) { quitLoop = true; break;}
    double sigmaNow = phaseSpacePtr->sigmaNow();

    // Tell if this event comes with weight from cross section.
    double sigmaWeight = 1.;
    if (!isLHA && !increaseMaximum && sigmaNow > sigmaMx)
      sigmaWeight = sigmaNow / sigmaMx;
    if ( lhaStrat < 0 && sigmaNow < 0.) sigmaWeight = -1.;
    if ( lhaStratAbs == 4) sigmaWeight = sigmaNow;

    // Also compensating weight from biased phase-space selection.
    double biasWeight = phaseSpacePtr->biasSelectionWeight();
    weightNow = sigmaWeight * biasWeight;
    infoPtr->setWeight( weightNow, lhaStrat);

    // Check that not negative cross section when not allowed.
    if (!allowNegSig) {
      if (sigmaNow < sigmaNeg) {
        infoPtr->errorMsg("Warning in ProcessContainer::trialProcess: neg"
          "ative cross section set 0", "for " +  sigmaProcessPtr->name );
        sigmaNeg = sigmaNow;
      }
      if (sigmaNow < 0.) sigmaNow = 0.;
    }

    // Cross section of process may come with a weight. Update sum.
    double sigmaAdd = sigmaNow * biasWeight;
    if (lhaStratAbs == 2 || lhaStratAbs == 3) sigmaAdd = sigmaSgn;
    sigmaSum  += sigmaAdd;
    sigma2Sum += pow2(sigmaAdd);

    // Check if maximum violated.
    newSigmaMx = phaseSpacePtr->newSigmaMax();
    if (newSigmaMx) sigmaMx = phaseSpacePtr->sigmaMax();

    // Select or reject trial point. Statistics.
    bool select = true;
    if (lhaStratAbs < 3) select
      = (newSigmaMx || rndmPtr->flat() * abs(sigmaMx) < abs(sigmaNow));
    if (select) {
      ++nSel;
      if (isLHA) {
        int codeLHANow = lhaUpPtr->idProcess();
        int iFill = -1;
        for (int i = 0; i < int(codeLHA.size()); ++i)
          if (codeLHANow == codeLHA[i]) iFill = i;
        if (iFill >= 0) ++nSelLHA[iFill];
      }
    }
    if (select || lhaStratAbs != 2) {quitLoop = !select; break;}

  }

  Benchmark_loopStop(trialProcess_loopTrialProcess);
  Benchmark_placeholder(trialProcess_dummyB)
  // seems to be weird un-indent here...
  if (quitLoop) return false;

  return true; // !!!!!
}

//--------------------------------------------------------------------------

// Accumulate statistics after user veto, including LHA code.

void ProcessContainer::accumulate() {

  ++nAcc;
  wtAccSum += weightNow;
  if (isLHA) {
    int codeLHANow = lhaUpPtr->idProcess();
    int iFill = -1;
    for (int i = 0; i < int(codeLHA.size()); ++i)
      if (codeLHANow == codeLHA[i]) iFill = i;
    if (iFill >= 0) ++nAccLHA[iFill];
  }

}

//--------------------------------------------------------------------------

// Pick flavours and colour flow of process.

bool ProcessContainer::constructState() {

  // Construct flavour and colours for selected event.
  if (isResolved && !isNonDiff) sigmaProcessPtr->pickInState();
  sigmaProcessPtr->setIdColAcol();

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Give the hard subprocess with kinematics.

bool ProcessContainer::constructProcess( Event& process, bool isHardest) {

  // Construct kinematics from selected phase space point.
  if (!phaseSpacePtr->finalKin()) return false;
  int nFin = sigmaProcessPtr->nFinal;

  // Basic info on process.
  if (isHardest) infoPtr->setType( name(), code(), nFin, isNonDiff,
    isResolved, isDiffA, isDiffB, isDiffC, isLHA);

  // Let hard process record begin with the event as a whole and
  // the two incoming beam particles.
  int idA = infoPtr->idA();
  if (abs(idA) == idRenameBeams) idA = 16;
  int idB = infoPtr->idB();
  if (abs(idB) == idRenameBeams) idB = -16;
  process.append( 90, -11, 0, 0, 0, 0, 0, 0,
    Vec4(0., 0., 0., infoPtr->eCM()), infoPtr->eCM(), 0. );
  process.append( idA, -12, 0, 0, 0, 0, 0, 0,
    Vec4(0., 0., infoPtr->pzA(), infoPtr->eA()), infoPtr->mA(), 0. );
  process.append( idB, -12, 0, 0, 0, 0, 0, 0,
    Vec4(0., 0., infoPtr->pzB(), infoPtr->eB()), infoPtr->mB(), 0. );

  // For nondiffractive process no interaction selected so far, so done.
  if (isNonDiff) return true;

  // Entries 3 and 4, now to be added, come from 1 and 2.
  process[1].daughter1(3);
  process[2].daughter1(4);
  double scale  = 0.;
  double scalup = 0.;

  // For DiffC entries 3 - 5 come jointly from 1 and 2 (to keep HepMC happy).
  if (isDiffC) {
    process[1].daughters(3, 5);
    process[2].daughters(3, 5);
  }

  // Insert the subprocess partons - resolved processes.
  int idRes = sigmaProcessPtr->idSChannel;
  if (isResolved && !isLHA) {

    // NOAM: Mothers and daughters without/with intermediate state.
    int m_M1 = 3;
    int m_M2 = 4;
    int m_D1 = 5;
    int m_D2 = (nFin == 1) ? 0 : 4 + nFin;
    if (idRes != 0) {
      m_M1   = 5;
      m_M2   = 0;
      m_D1   = 5;
      m_D2   = 0;
    }

    // Find scale from which to begin MPI/ISR/FSR evolution.
    scale = sqrt(Q2Fac());
    process.scale( scale );

    // Loop over incoming and outgoing partons.
    int colOffset = process.lastColTag();
    for (int i = 1; i <= 2 + nFin; ++i) {

      // Read out particle info from SigmaProcess object.
      int id        = sigmaProcessPtr->idSave[i];
      int status    = (i <= 2) ? -21 : 23;
      int mother1   = (i <= 2) ? i : m_M1 ;
      int mother2   = (i <= 2) ? 0 : m_M2 ;
      int daughter1 = (i <= 2) ? m_D1 : 0;
      int daughter2 = (i <= 2) ? m_D2 : 0;
      int col       = sigmaProcessPtr->colSave[i];
      if      (col > 0) col += colOffset;
      else if (col < 0) col -= colOffset;
      int acol      = sigmaProcessPtr->acolSave[i];
      if      (acol > 0) acol += colOffset;
      else if (acol < 0) acol -= colOffset;

      // Append to process record.
      int iNow = process.append( id, status, mother1, mother2,
        daughter1, daughter2, col, acol, phaseSpacePtr->p(i),
        phaseSpacePtr->m(i), scale);

      // NOAM: If there is an intermediate state, insert the it in
      // the process record after the two incoming particles.
      if (i == 2 && idRes != 0) {

        // Sign of intermediate state: go by charge.
        if (particleDataPtr->hasAnti(idRes)
          && process[3].chargeType() + process[4].chargeType() < 0)
          idRes *= -1;

        // The colour configuration of the intermediate state has to be
        // resolved separately.
        col         = 0;
        acol        = 0;
        int m_col1  = sigmaProcessPtr->colSave[1];
        int m_acol1 = sigmaProcessPtr->acolSave[1];
        int m_col2  = sigmaProcessPtr->colSave[2];
        int m_acol2 = sigmaProcessPtr->acolSave[2];
        if (m_col1 == m_acol2 && m_col2 != m_acol1) {
          col       = m_col2;
          acol       = m_acol1;
        } else if (m_col2 == m_acol1 && m_col1 != m_acol2) {
          col        = m_col1;
          acol       = m_acol2;
        }
        if      ( col > 0)  col += colOffset;
        else if ( col < 0)  col -= colOffset;
        if      (acol > 0) acol += colOffset;
        else if (acol < 0) acol -= colOffset;

        // Insert the intermediate state into the event record.
        Vec4 pIntMed = phaseSpacePtr->p(1) + phaseSpacePtr->p(2);
        process.append( idRes, -22, 3, 4,  6, 5 + nFin, col, acol,
          pIntMed, pIntMed.mCalc(), scale);
      }

      // Pick lifetime where relevant, else not.
      if (process[iNow].tau0() > 0.) process[iNow].tau(
        process[iNow].tau0() * rndmPtr->exp() );
    }
  }

  // Insert the outgoing particles - unresolved processes.
  else if (!isLHA) {
    int id3     = sigmaProcessPtr->idSave[3];
    int status3 = (id3 == process[1].id()) ? 14 : 15;
    process.append( id3, status3, 1, 0, 0, 0, 0, 0,
      phaseSpacePtr->p(3), phaseSpacePtr->m(3));
    int id4     = sigmaProcessPtr->idSave[4];
    int status4 = (id4 == process[2].id()) ? 14 : 15;
    process.append( id4, status4, 2, 0, 0, 0, 0, 0,
      phaseSpacePtr->p(4), phaseSpacePtr->m(4));

    // For central diffraction: two scattered protons inserted so far,
    // but modify mothers, add also centrally-produced hadronic system.
    if (isDiffC) {
      process[3].mothers( 1, 2);
      process[4].mothers( 1, 2);
      int id5     = sigmaProcessPtr->idSave[5];
      int status5 = 15;
      process.append( id5, status5, 1, 2, 0, 0, 0, 0,
        phaseSpacePtr->p(5), phaseSpacePtr->m(5));
    }
  }

  // Insert the outgoing particles - Les Houches Accord processes.
  else {

    // Since LHA partons may be out of order, determine correct one.
    // (Recall that zeroth particle is empty.)
    vector<int> newPos;
    newPos.reserve(lhaUpPtr->sizePart());
    newPos.push_back(0);
    for (int iNew = 0; iNew < lhaUpPtr->sizePart(); ++iNew) {
      // For iNew == 0 look for the two incoming partons, then for
      // partons having them as mothers, and so on layer by layer.
      for (int i = 1; i < lhaUpPtr->sizePart(); ++i)
        if (lhaUpPtr->mother1(i) == newPos[iNew]) newPos.push_back(i);
      if (int(newPos.size()) <= iNew) break;
    }

    // Find scale from which to begin MPI/ISR/FSR evolution.
    scalup = lhaUpPtr->scale();
    scale  = scalup;
    double scalePr = (scale < 0.) ? sqrt(Q2Fac()) : scale;
    process.scale( scalePr);

    // Copy over info from LHA event to process, in proper order.
    Vec4 pSumOut;
    for (int i = 1; i < lhaUpPtr->sizePart(); ++i) {
      int iOld = newPos[i];
      int id = lhaUpPtr->id(iOld);
      if (i == 1 && abs(id) == idRenameBeams) id = 16;
      if (i == 2 && abs(id) == idRenameBeams) id = -16;
      int idAbs = abs(id);

      // Translate from LHA status codes.
      int lhaStatus =  lhaUpPtr->status(iOld);
      int status = -21;
      if (lhaStatus == 2 || lhaStatus == 3) status = -22;
      if (lhaStatus == 1) status = 23;

      // Find where mothers have been moved by reordering.
      int mother1Old = lhaUpPtr->mother1(iOld);
      int mother2Old = lhaUpPtr->mother2(iOld);
      int mother1 = 0;
      int mother2 = 0;
      for (int im = 1; im < i; ++im) {
        if (mother1Old == newPos[im]) mother1 = im + 2;
        if (mother2Old == newPos[im]) mother2 = im + 2;
      }
      if (i <= 2) mother1 = i;

      // Ensure that second mother = 0 except for bona fide carbon copies.
      if (mother1 > 0 && mother2 == mother1) {
        int sister1 = process[mother1].daughter1();
        int sister2 = process[mother1].daughter2();
        if (sister2 != sister1 && sister2 != 0) mother2 = 0;
      }

      // Find daughters and where they have been moved by reordering.
      // (Values shifted two steps to account for inserted beams.)
      int daughter1 = 0;
      int daughter2 = 0;
      for (int im = i + 1; im < lhaUpPtr->sizePart(); ++im) {
        if (lhaUpPtr->mother1(newPos[im]) == iOld
          || lhaUpPtr->mother2(newPos[im]) == iOld) {
          if (daughter1 == 0 || im + 2 < daughter1) daughter1 = im + 2;
          if (daughter2 == 0 || im + 2 > daughter2) daughter2 = im + 2;
        }
      }
      // For 2 -> 1 hard scatterings reset second daughter to 0.
      if (daughter2 == daughter1) daughter2 = 0;

      // Colour trivial, except reset irrelevant colour indices.
      int colType = particleDataPtr->colType(id);
      int col1   = (colType == 1 || colType == 2 || abs(colType) == 3)
                 ? lhaUpPtr->col1(iOld) : 0;
      int col2   = (colType == -1 || colType == 2 || abs(colType) == 3)
                 ?  lhaUpPtr->col2(iOld) : 0;

      // Copy momentum, ensure lepton masses and consistent (E, p m) set.
      double px  = lhaUpPtr->px(iOld);
      double py  = lhaUpPtr->py(iOld);
      double pz  = lhaUpPtr->pz(iOld);
      double e   = lhaUpPtr->e(iOld);
      double m   = lhaUpPtr->m(iOld);
      for (int idL = 0; idL < 3; ++idL)
        if (idAbs == idLep[idL] && setLeptonMass > 0 && (setLeptonMass == 2
          || m < 0.9 * mLep[idL] || m > 1.1 * mLep[idL])) m = mLep[idL];
      if (mRecalculate > 0. && m > mRecalculate)
        m = sqrtpos( e*e - px*px - py*py - pz*pz);
      else e = sqrt( m*m + px*px + py*py + pz*pz);

      // Momentum sum for outgoing particles.
      if (matchInOut && i > 2 && lhaStatus == 1)
        pSumOut += Vec4( px, py, pz, e);

      // Polarization.
      double pol = lhaUpPtr->spin(iOld);

      // Allow scale setting for generic partons.
      double scaleShow = lhaUpPtr->scale(iOld);

      // For resonance decay products use resonance mass as scale.
      double scaleNow = scalePr;
      if (mother1 > 4 && !useStrictLHEFscales) scaleNow = process[mother1].m();
      if (scaleShow >= 0.0) scaleNow = scaleShow;

      // Store Les Houches Accord partons.
      int iNow = process.append( id, status, mother1, mother2, daughter1,
        daughter2, col1, col2, Vec4(px, py, pz, e), m, scaleNow, pol);

      // Check if need to store lifetime.
      double tau = lhaUpPtr->tau(iOld);
      if ( (setLifetime == 1 && idAbs == 15) || setLifetime == 2)
         tau = process[iNow].tau0() * rndmPtr->exp();
      if (tau > 0.) process[iNow].tau(tau);
    }

    // Reassign momenta and masses for incoming partons.
    if (matchInOut) {
      double e1 = 0.5 * (pSumOut.e() + pSumOut.pz());
      double e2 = 0.5 * (pSumOut.e() - pSumOut.pz());
      process[3].pz( e1);
      process[3].e(  e1);
      process[3].m(  0.);
      process[4].pz(-e2);
      process[4].e(  e2);
      process[4].m(  0.);
    }
  }

  // Loop through decay chains and set secondary vertices when needed.
  for (int i = 3; i < process.size(); ++i) {
    int iMother  = process[i].mother1();

    // If sister to already assigned vertex then assign same.
    if ( process[i - 1].mother1() == iMother && process[i - 1].hasVertex() )
      process[i].vProd( process[i - 1].vProd() );

    // Else if mother already has vertex and/or lifetime then assign.
    else if ( process[iMother].hasVertex() || process[iMother].tau() > 0.)
      process[i].vProd( process[iMother].vDec() );
  }

  // Further info on process. Reset quantities that may or may not be known.
  int    id1Now  = process[3].id();
  int    id2Now  = process[4].id();
  int    id1pdf  = 0;
  int    id2pdf  = 0;
  double x1pdf   = 0.;
  double x2pdf   = 0.;
  double pdf1    = 0.;
  double pdf2    = 0.;
  double tHat    = 0.;
  double uHat    = 0.;
  double pTHatL  = 0.;
  double m3      = 0.;
  double m4      = 0.;
  double theta   = 0.;
  double phi     = 0.;
  double x1Now, x2Now, Q2FacNow, alphaEM, alphaS, Q2Ren, sHat;

  // Internally generated and stored information.
  if (!isLHA) {
    id1pdf       = id1Now;
    id2pdf       = id2Now;
    x1Now        = phaseSpacePtr->x1();
    x2Now        = phaseSpacePtr->x2();
    x1pdf        = x1Now;
    x2pdf        = x2Now;
    pdf1         = sigmaProcessPtr->pdf1Save;
    pdf2         = sigmaProcessPtr->pdf2Save;
    Q2FacNow     = sigmaProcessPtr->Q2FacSave;
    alphaEM      = sigmaProcessPtr->alpEM;
    alphaS       = sigmaProcessPtr->alpS;
    Q2Ren        = sigmaProcessPtr->Q2RenSave;
    sHat         = phaseSpacePtr->sHat();
    tHat         = phaseSpacePtr->tHat();
    uHat         = phaseSpacePtr->uHat();
    pTHatL       = phaseSpacePtr->pTHat();
    m3           = phaseSpacePtr->m(3);
    m4           = phaseSpacePtr->m(4);
    theta        = phaseSpacePtr->thetaHat();
    phi          = phaseSpacePtr->phiHat();
  }

  // Les Houches Accord process partly available, partly to be constructed.
  else {
    x1Now        = 2. * process[3].e() / infoPtr->eCM();
    x2Now        = 2. * process[4].e() / infoPtr->eCM();
    Q2FacNow     = (scale < 0.) ? sigmaProcessPtr->Q2FacSave : pow2(scale);
    alphaEM      = lhaUpPtr->alphaQED();
    if (alphaEM < 0.001) alphaEM = sigmaProcessPtr->alpEM;
    alphaS       = lhaUpPtr->alphaQCD();
    if (alphaS  < 0.001) alphaS  = sigmaProcessPtr->alpS;
    Q2Ren        = (scale < 0.) ? sigmaProcessPtr->Q2RenSave : pow2(scale);
    Vec4 pSum    = process[3].p() + process[4].p();
    sHat         = pSum * pSum;
    int nFinLH   = 0;
    double pTLH  = 0.;
    for (int i = 5; i < process.size(); ++i)
    if (process[i].mother1() == 3 && process[i].mother2() == 4) {
      ++nFinLH;
      pTLH      += process[i].pT();
    }
    if (nFinLH > 0) pTHatL = pTLH / nFinLH;

    // Read info on parton densities if provided.
    id1pdf       = lhaUpPtr->id1pdf();
    id2pdf       = lhaUpPtr->id2pdf();
    x1pdf        = lhaUpPtr->x1pdf();
    x2pdf        = lhaUpPtr->x2pdf();
    if (lhaUpPtr->pdfIsSet()) {
      pdf1       = lhaUpPtr->pdf1();
      pdf2       = lhaUpPtr->pdf2();
      Q2FacNow   = pow2(lhaUpPtr->scalePDF());
    }

    // Reconstruct kinematics of 2 -> 2 processes from momenta.
    if (nFin == 2) {
      Vec4 pDifT = process[3].p() - process[5].p();
      tHat       = pDifT * pDifT;
      Vec4 pDifU = process[3].p() - process[6].p();
      uHat       = pDifU * pDifU;
      pTHatL     = process[5].pT();
      m3         = process[5].m();
      m4         = process[6].m();
      Vec4 p5    = process[5].p();
      p5.bstback(pSum);
      theta      = p5.theta();
      phi        = process[5].phi();
    }
  }

  // Store information.
  if (isHardest) {
    infoPtr->setPDFalpha( 0, id1pdf, id2pdf, x1pdf, x2pdf, pdf1, pdf2,
      Q2FacNow, alphaEM, alphaS, Q2Ren, scalup);
    infoPtr->setKin( 0, id1Now, id2Now, x1Now, x2Now, sHat, tHat, uHat,
      pTHatL, m3, m4, theta, phi);
  }
  infoPtr->setTypeMPI( code(), pTHatL);

  // For Les Houches event store subprocess classification.
  if (isLHA) {
    int codeSub  = lhaUpPtr->idProcess();
    ostringstream nameSub;
    nameSub << "user process " << codeSub;
    infoPtr->setSubType( 0, nameSub.str(), codeSub, nFin);
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Give the hard resonance decay chain from Les Houches input.
// Note: mildly modified copy of some constructProcess code; to streamline.

bool ProcessContainer::constructDecays( Event& process) {

  // Let hard process record begin with the event as a whole.
  process.clear();
  process.append( 90, -11, 0, 0, 0, 0, 0, 0, Vec4(0., 0., 0., 0.), 0., 0. );

  // Since LHA partons may be out of order, determine correct one.
  // (Recall that zeroth particle is empty.)
  vector<int> newPos;
  newPos.reserve(lhaUpPtr->sizePart());
  newPos.push_back(0);
  for (int iNew = 0; iNew < lhaUpPtr->sizePart(); ++iNew) {
    // For iNew == 0 look for the two incoming partons, then for
    // partons having them as mothers, and so on layer by layer.
    for (int i = 1; i < lhaUpPtr->sizePart(); ++i)
      if (lhaUpPtr->mother1(i) == newPos[iNew]) newPos.push_back(i);
    if (int(newPos.size()) <= iNew) break;
  }

  // Find scale from which to begin MPI/ISR/FSR evolution.
  double scale = lhaUpPtr->scale();
  process.scale( scale);
  Vec4 pSum;

  // Copy over info from LHA event to process, in proper order.
  for (int i = 1; i < lhaUpPtr->sizePart(); ++i) {
    int iOld = newPos[i];
    int id = lhaUpPtr->id(iOld);

    // Translate from LHA status codes.
    int lhaStatus =  lhaUpPtr->status(iOld);
    int status = -21;
    if (lhaStatus == 2 || lhaStatus == 3) status = -22;
    if (lhaStatus == 1) status = 23;

    // Find where mothers have been moved by reordering.
    int mother1Old = lhaUpPtr->mother1(iOld);
    int mother2Old = lhaUpPtr->mother2(iOld);
    int mother1 = 0;
    int mother2 = 0;
    for (int im = 1; im < i; ++im) {
      if (mother1Old == newPos[im]) mother1 = im;
      if (mother2Old == newPos[im]) mother2 = im;
    }

    // Ensure that second mother = 0 except for bona fide carbon copies.
    if (mother1 > 0 && mother2 == mother1) {
      int sister1 = process[mother1].daughter1();
      int sister2 = process[mother1].daughter2();
      if (sister2 != sister1 && sister2 != 0) mother2 = 0;
    }

    // Find daughters and where they have been moved by reordering.
    int daughter1 = 0;
    int daughter2 = 0;
    for (int im = i + 1; im < lhaUpPtr->sizePart(); ++im) {
      if (lhaUpPtr->mother1(newPos[im]) == iOld
        || lhaUpPtr->mother2(newPos[im]) == iOld) {
        if (daughter1 == 0 || im < daughter1) daughter1 = im;
        if (daughter2 == 0 || im > daughter2) daughter2 = im;
      }
    }
    // For 2 -> 1 hard scatterings reset second daughter to 0.
    if (daughter2 == daughter1) daughter2 = 0;

    // Colour trivial, except reset irrelevant colour indices.
    int colType = particleDataPtr->colType(id);
    int col1   = (colType == 1 || colType == 2 || abs(colType) == 3)
               ? lhaUpPtr->col1(iOld) : 0;
    int col2   = (colType == -1 || colType == 2 || abs(colType) == 3)
               ?  lhaUpPtr->col2(iOld) : 0;

    // Momentum trivial.
    double px  = lhaUpPtr->px(iOld);
    double py  = lhaUpPtr->py(iOld);
    double pz  = lhaUpPtr->pz(iOld);
    double e   = lhaUpPtr->e(iOld);
    double m   = lhaUpPtr->m(iOld);
    if (status > 0) pSum += Vec4( px, py, pz, e);

    // Polarization.
    double pol = lhaUpPtr->spin(iOld);

    // For resonance decay products use resonance mass as scale.
    double scaleNow = scale;
    if (mother1 > 0) scaleNow = process[mother1].m();

    // Store Les Houches Accord partons.
    int iNow = process.append( id, status, mother1, mother2, daughter1,
      daughter2, col1, col2, Vec4(px, py, pz, e), m, scaleNow, pol);

    // Check if need to store lifetime.
    double tau = lhaUpPtr->tau(iOld);
    if ( (setLifetime == 1 && abs(id) == 15) || setLifetime == 2)
       tau = process[iNow].tau0() * rndmPtr->exp();
    if (tau > 0.) process[iNow].tau(tau);
  }

  // Update four-momentum of system as a whole.
  process[0].p( pSum);
  process[0].m( pSum.mCalc());

  // Loop through decay chains and set secondary vertices when needed.
  for (int i = 1; i < process.size(); ++i) {
    int iMother  = process[i].mother1();

    // If sister to already assigned vertex then assign same.
    if ( process[i - 1].mother1() == iMother && process[i - 1].hasVertex()
      && i > 1) process[i].vProd( process[i - 1].vProd() );

    // Else if mother already has vertex and/or lifetime then assign.
    else if ( process[iMother].hasVertex() || process[iMother].tau() > 0.)
      process[i].vProd( process[iMother].vDec() );
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Handle resonance decays.

bool ProcessContainer::decayResonances( Event& process) {

  Benchmark_start(decayResonances);
  Benchmark_start(decayResonances_setup);

  // Save current event-record size and status codes.
  process.saveSize();
  vector<int> statusSave( process.size());
  for (int i = 0; i < process.size(); ++i)
    statusSave[i] = process[i].status();
  bool physical    = true;
  bool newChain    = false;
  bool newFlavours = false;

  Benchmark_stop(decayResonances_setup);

  // Do loop over user veto.
  do {

    // Do sequential chain of uncorrelated isotropic decays.
    do {
      Benchmark_start(decayResonances_next);
      physical = resDecaysPtr->next( process);
      Benchmark_stop(decayResonances_next);
      
      if (!physical) return false;

      Benchmark_start(decayResonances_weightDecayFlav);
      // Check whether flavours should be correlated.
      // (Currently only relevant for f fbar -> gamma*/Z0 gamma*/Z0.)
      newFlavours = ( sigmaProcessPtr->weightDecayFlav( process)
                    < rndmPtr->flat() );
      Benchmark_stop(decayResonances_weightDecayFlav);

      Benchmark_start(decayResonances_resetForVeto);
      // Reset the decay chains if have to redo.
      if (newFlavours) {
        process.restoreSize();
        for (int i = 0; i < process.size(); ++i)
          process[i].status( statusSave[i]);
      }
      Benchmark_stop(decayResonances_resetForVeto);

    // Loop back where required to generate new decays with new flavours.
    } while (newFlavours);

    Benchmark_start(decayResonances_decayKinematics);

    // Correct to nonisotropic decays.
    phaseSpacePtr->decayKinematics( process);

    Benchmark_stop(decayResonances_decayKinematics);

    Benchmark_start(decayResonances_rest);

    // Optionally user hooks check/veto on decay chain.
    if (canVetoResDecay)
      newChain = userHooksPtr->doVetoResonanceDecays( process);

    // Reset the decay chains if have to redo.
    if (newChain) {
      process.restoreSize();
      for (int i = 0; i < process.size(); ++i)
        process[i].status( statusSave[i]);
    }

  // Loop back where required to generate new decay chain.
  } while(newChain);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Reset event generation statistics; but NOT maximum of cross section.

void ProcessContainer::reset() {

  nTry      = 0;
  nSel      = 0;
  nAcc      = 0;
  nTryStat  = 0;
  sigmaSum  = 0.;
  sigma2Sum = 0.;
  sigmaNeg  = 0.;
  sigmaAvg  = 0.;
  sigmaFin  = 0.;
  deltaFin  = 0.;
  wtAccSum  = 0.;

}

//--------------------------------------------------------------------------

// Estimate integrated cross section and its uncertainty.

void ProcessContainer::sigmaDelta() {

  // Initial values. No analysis meaningful unless accepted events.
  nTryStat = nTry;
  sigmaAvg = 0.;
  sigmaFin = 0.;
  deltaFin = 0.;
  if (nAcc == 0) return;

  // Average value. No error analysis unless at least two events.
  double nTryInv  = 1. / nTry;
  double nSelInv  = 1. / nSel;
  double nAccInv  = 1. / nAcc;
  sigmaAvg        = sigmaSum * nTryInv;
  double fracAcc  = nAcc * nSelInv;
  sigmaFin        = sigmaAvg * fracAcc;
  deltaFin        = sigmaFin;
  if (nAcc == 1) return;

  // Estimated error. Quadratic sum of cross section term and
  // binomial from accept/reject step.
  double delta2Sig   = (lhaStratAbs != 3)
    ? (sigma2Sum * nTryInv - pow2(sigmaAvg)) * nTryInv / pow2(sigmaAvg)
    : pow2( lhaUpPtr->xErrSum() / lhaUpPtr->xSecSum());
  double delta2Veto  = (nSel - nAcc) * nAccInv * nSelInv;
  double delta2Sum   = delta2Sig + delta2Veto;
  deltaFin           = sqrtpos(delta2Sum) * sigmaFin;

}

//==========================================================================

// SetupContainer class.
// Turns list of user-desired processes into a vector of containers.

//--------------------------------------------------------------------------

// Main routine to initialize list of processes.

bool SetupContainers::init(vector<ProcessContainer*>& containerPtrs,
       Info *infoPtr, Settings& settings, ParticleData* particleDataPtr,
       Couplings* couplings) {

  // Reset process list, if filled in previous subrun.
  if (containerPtrs.size() > 0) {
    for (int i = 0; i < int(containerPtrs.size()); ++i)
      delete containerPtrs[i];
    containerPtrs.clear();
  }
  SigmaProcess* sigmaPtr;

  // Set up requested objects for soft QCD processes.
  bool softQCD = settings.get(Flag::SoftQCD_all);
  bool inelastic = settings.get(Flag::SoftQCD_inelastic);
  if (softQCD || inelastic || settings.get(Flag::SoftQCD_nonDiffractive)) {
    sigmaPtr = new Sigma0nonDiffractive;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (softQCD || settings.get(Flag::SoftQCD_elastic)) {
    sigmaPtr = new Sigma0AB2AB;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (softQCD || inelastic || settings.get(Flag::SoftQCD_singleDiffractive)) {
    sigmaPtr = new Sigma0AB2XB;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma0AB2AX;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (softQCD || inelastic || settings.get(Flag::SoftQCD_doubleDiffractive)) {
    sigmaPtr = new Sigma0AB2XX;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (softQCD || inelastic || settings.get(Flag::SoftQCD_centralDiffractive)) {
    sigmaPtr = new Sigma0AB2AXB;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Set up requested objects for hard QCD processes.
  bool hardQCD = settings.get(Flag::HardQCD_all);
  if (hardQCD || settings.get(Flag::HardQCD_gg2gg)) {
    sigmaPtr = new Sigma2gg2gg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD || settings.get(Flag::HardQCD_gg2qqbar)) {
    sigmaPtr = new Sigma2gg2qqbar;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD || settings.get(Flag::HardQCD_qg2qg)) {
    sigmaPtr = new Sigma2qg2qg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD || settings.get(Flag::HardQCD_qq2qq)) {
    sigmaPtr = new Sigma2qq2qq;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD || settings.get(Flag::HardQCD_qqbar2gg)) {
    sigmaPtr = new Sigma2qqbar2gg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD || settings.get(Flag::HardQCD_qqbar2qqbarNew)) {
    sigmaPtr = new Sigma2qqbar2qqbarNew;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Set up requested objects for c cbar and b bbar, also hard QCD.
  bool hardccbar = settings.get(Flag::HardQCD_hardccbar);
  if (hardQCD || hardccbar || settings.get(Flag::HardQCD_gg2ccbar)) {
    sigmaPtr = new Sigma2gg2QQbar(4, 121);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD || hardccbar || settings.get(Flag::HardQCD_qqbar2ccbar)) {
    sigmaPtr = new Sigma2qqbar2QQbar(4, 122);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  bool hardbbbar = settings.get(Flag::HardQCD_hardbbbar);
  if (hardQCD || hardbbbar || settings.get(Flag::HardQCD_gg2bbbar)) {
    sigmaPtr = new Sigma2gg2QQbar(5, 123);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD || hardbbbar || settings.get(Flag::HardQCD_qqbar2bbbar)) {
    sigmaPtr = new Sigma2qqbar2QQbar(5, 124);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Set up requested objects for hard QCD 2 -> 3 processes.
  bool hardQCD3parton = settings.get(Flag::HardQCD_3parton);
  if (hardQCD3parton || settings.get(Flag::HardQCD_gg2ggg)) {
    sigmaPtr = new Sigma3gg2ggg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD3parton || settings.get(Flag::HardQCD_qqbar2ggg)) {
    sigmaPtr = new Sigma3qqbar2ggg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD3parton || settings.get(Flag::HardQCD_qg2qgg)) {
    sigmaPtr = new Sigma3qg2qgg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD3parton || settings.get(Flag::HardQCD_qq2qqgDiff)) {
    sigmaPtr = new Sigma3qq2qqgDiff;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD3parton || settings.get(Flag::HardQCD_qq2qqgSame)) {
    sigmaPtr = new Sigma3qq2qqgSame;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD3parton || settings.get(Flag::HardQCD_qqbar2qqbargDiff)) {
    sigmaPtr = new Sigma3qqbar2qqbargDiff;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD3parton || settings.get(Flag::HardQCD_qqbar2qqbargSame)) {
    sigmaPtr = new Sigma3qqbar2qqbargSame;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD3parton || settings.get(Flag::HardQCD_gg2qqbarg)) {
    sigmaPtr = new Sigma3gg2qqbarg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD3parton || settings.get(Flag::HardQCD_qg2qqqbarDiff)) {
    sigmaPtr = new Sigma3qg2qqqbarDiff;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD3parton || settings.get(Flag::HardQCD_qg2qqqbarSame)) {
    sigmaPtr = new Sigma3qg2qqqbarSame;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Set up requested objects for prompt photon processes.
  bool promptPhotons = settings.get(Flag::PromptPhoton_all);
  if (promptPhotons
    || settings.get(Flag::PromptPhoton_qg2qgamma)) {
    sigmaPtr = new Sigma2qg2qgamma;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (promptPhotons
    || settings.get(Flag::PromptPhoton_qqbar2ggamma)) {
    sigmaPtr = new Sigma2qqbar2ggamma;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (promptPhotons
    || settings.get(Flag::PromptPhoton_gg2ggamma)) {
    sigmaPtr = new Sigma2gg2ggamma;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (promptPhotons
    || settings.get(Flag::PromptPhoton_ffbar2gammagamma)) {
    sigmaPtr = new Sigma2ffbar2gammagamma;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (promptPhotons
    || settings.get(Flag::PromptPhoton_gg2gammagamma)) {
    sigmaPtr = new Sigma2gg2gammagamma;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Set up requested objects for weak gauge boson t-channel exchange.
  bool weakBosonExchanges = settings.get(Flag::WeakBosonExchange_all);
  if (weakBosonExchanges
    || settings.get(Flag::WeakBosonExchange_ff2ff_t_gmZ_)) {
    sigmaPtr = new Sigma2ff2fftgmZ;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (weakBosonExchanges
    || settings.get(Flag::WeakBosonExchange_ff2ff_t_W_)) {
    sigmaPtr = new Sigma2ff2fftW;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Set up requested objects for weak gauge boson processes.
  bool weakSingleBosons = settings.get(Flag::WeakSingleBoson_all);
  if (weakSingleBosons
    || settings.get(Flag::WeakSingleBoson_ffbar2gmZ)) {
    sigmaPtr = new Sigma1ffbar2gmZ;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (weakSingleBosons
    || settings.get(Flag::WeakSingleBoson_ffbar2W)) {
    sigmaPtr = new Sigma1ffbar2W;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Set up requested object for s-channel gamma exchange.
  // Subset of gamma*/Z0 above, intended for multiparton interactions.
  if (settings.get(Flag::WeakSingleBoson_ffbar2ffbar_s_gm_)) {
    sigmaPtr = new Sigma2ffbar2ffbarsgm;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Set up requested object for s-channel gamma*/Z0 or W+- exchange.
  if (settings.get(Flag::WeakSingleBoson_ffbar2ffbar_s_gmZ_)) {
    sigmaPtr = new Sigma2ffbar2ffbarsgmZ;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::WeakSingleBoson_ffbar2ffbar_s_W_)) {
    sigmaPtr = new Sigma2ffbar2ffbarsW;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Set up requested objects for weak gauge boson pair processes.
  bool weakDoubleBosons = settings.get(Flag::WeakDoubleBoson_all);
  if (weakDoubleBosons
    || settings.get(Flag::WeakDoubleBoson_ffbar2gmZgmZ)) {
    sigmaPtr = new Sigma2ffbar2gmZgmZ;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (weakDoubleBosons
    || settings.get(Flag::WeakDoubleBoson_ffbar2ZW)) {
    sigmaPtr = new Sigma2ffbar2ZW;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (weakDoubleBosons
    || settings.get(Flag::WeakDoubleBoson_ffbar2WW)) {
    sigmaPtr = new Sigma2ffbar2WW;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Set up requested objects for weak gauge boson + parton processes.
  bool weakBosonAndPartons = settings.get(Flag::WeakBosonAndParton_all);
  if (weakBosonAndPartons
    || settings.get(Flag::WeakBosonAndParton_qqbar2gmZg)) {
    sigmaPtr = new Sigma2qqbar2gmZg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (weakBosonAndPartons
    || settings.get(Flag::WeakBosonAndParton_qg2gmZq)) {
    sigmaPtr = new Sigma2qg2gmZq;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (weakBosonAndPartons
    || settings.get(Flag::WeakBosonAndParton_ffbar2gmZgm)) {
    sigmaPtr = new Sigma2ffbar2gmZgm;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (weakBosonAndPartons
    || settings.get(Flag::WeakBosonAndParton_fgm2gmZf)) {
    sigmaPtr = new Sigma2fgm2gmZf;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (weakBosonAndPartons
    || settings.get(Flag::WeakBosonAndParton_qqbar2Wg)) {
    sigmaPtr = new Sigma2qqbar2Wg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (weakBosonAndPartons
    || settings.get(Flag::WeakBosonAndParton_qg2Wq)) {
    sigmaPtr = new Sigma2qg2Wq;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (weakBosonAndPartons
    || settings.get(Flag::WeakBosonAndParton_ffbar2Wgm)) {
    sigmaPtr = new Sigma2ffbar2Wgm;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (weakBosonAndPartons
    || settings.get(Flag::WeakBosonAndParton_fgm2Wf)) {
    sigmaPtr = new Sigma2fgm2Wf;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Set up requested objects for photon collision processes.
  bool photonCollisions = settings.get(Flag::PhotonCollision_all);
  if (photonCollisions || settings.get(Flag::PhotonCollision_gmgm2qqbar)) {
    sigmaPtr = new Sigma2gmgm2ffbar(1, 261);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (photonCollisions || settings.get(Flag::PhotonCollision_gmgm2ccbar)) {
    sigmaPtr = new Sigma2gmgm2ffbar(4, 262);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (photonCollisions || settings.get(Flag::PhotonCollision_gmgm2bbbar)) {
    sigmaPtr = new Sigma2gmgm2ffbar(5, 263);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (photonCollisions || settings.get(Flag::PhotonCollision_gmgm2ee)) {
    sigmaPtr = new Sigma2gmgm2ffbar(11, 264);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (photonCollisions || settings.get(Flag::PhotonCollision_gmgm2mumu)) {
    sigmaPtr = new Sigma2gmgm2ffbar(13, 265);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (photonCollisions || settings.get(Flag::PhotonCollision_gmgm2tautau)) {
    sigmaPtr = new Sigma2gmgm2ffbar(15, 266);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Set up requested objects for onia production.
  charmonium = SigmaOniaSetup(infoPtr, &settings, particleDataPtr, 4);
  bottomonium = SigmaOniaSetup(infoPtr, &settings, particleDataPtr, 5);
  vector<SigmaProcess*> charmoniumSigmaPtrs, bottomoniumSigmaPtrs;
  charmonium.setupSigma2gg(charmoniumSigmaPtrs);
  charmonium.setupSigma2qg(charmoniumSigmaPtrs);
  charmonium.setupSigma2qq(charmoniumSigmaPtrs);
  bottomonium.setupSigma2gg(bottomoniumSigmaPtrs);
  bottomonium.setupSigma2qg(bottomoniumSigmaPtrs);
  bottomonium.setupSigma2qq(bottomoniumSigmaPtrs);
  for (unsigned int i = 0; i < charmoniumSigmaPtrs.size(); ++i)
    containerPtrs.push_back( new ProcessContainer(charmoniumSigmaPtrs[i]) );
  for (unsigned int i = 0; i < bottomoniumSigmaPtrs.size(); ++i)
    containerPtrs.push_back( new ProcessContainer(bottomoniumSigmaPtrs[i]) );

  // Set up requested objects for top production.
  bool tops = settings.get(Flag::Top_all);
  if (tops || settings.get(Flag::Top_gg2ttbar)) {
    sigmaPtr = new Sigma2gg2QQbar(6, 601);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (tops || settings.get(Flag::Top_qqbar2ttbar)) {
    sigmaPtr = new Sigma2qqbar2QQbar(6, 602);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (tops || settings.get(Flag::Top_qq2tq_t_W_)) {
    sigmaPtr = new Sigma2qq2QqtW(6, 603);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (tops || settings.get(Flag::Top_ffbar2ttbar_s_gmZ_)) {
    sigmaPtr = new Sigma2ffbar2FFbarsgmZ(6, 604);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (tops || settings.get(Flag::Top_ffbar2tqbar_s_W_)) {
    sigmaPtr = new Sigma2ffbar2FfbarsW(6, 0, 605);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (tops || settings.get(Flag::Top_gmgm2ttbar)) {
    sigmaPtr = new Sigma2gmgm2ffbar(6, 606);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Set up requested objects for fourth-generation b' production
  bool bPrimes = settings.get(Flag::FourthBottom_all);
  if (bPrimes || settings.get(Flag::FourthBottom_gg2bPrimebPrimebar)) {
    sigmaPtr = new Sigma2gg2QQbar(7, 801);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (bPrimes || settings.get(Flag::FourthBottom_qqbar2bPrimebPrimebar)) {
    sigmaPtr = new Sigma2qqbar2QQbar(7, 802);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (bPrimes || settings.get(Flag::FourthBottom_qq2bPrimeq_t_W_)) {
    sigmaPtr = new Sigma2qq2QqtW(7, 803);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (bPrimes || settings.get(Flag::FourthBottom_ffbar2bPrimebPrimebar_s_gmZ_)) {
    sigmaPtr = new Sigma2ffbar2FFbarsgmZ(7, 804);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (bPrimes || settings.get(Flag::FourthBottom_ffbar2bPrimeqbar_s_W_)) {
    sigmaPtr = new Sigma2ffbar2FfbarsW(7, 0, 805);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (bPrimes || settings.get(Flag::FourthBottom_ffbar2bPrimetbar_s_W_)) {
    sigmaPtr = new Sigma2ffbar2FfbarsW(7, 6, 806);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Set up requested objects for fourth-generation t' production
  bool tPrimes = settings.get(Flag::FourthTop_all);
  if (tPrimes || settings.get(Flag::FourthTop_gg2tPrimetPrimebar)) {
    sigmaPtr = new Sigma2gg2QQbar(8, 821);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (tPrimes || settings.get(Flag::FourthTop_qqbar2tPrimetPrimebar)) {
    sigmaPtr = new Sigma2qqbar2QQbar(8, 822);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (tPrimes || settings.get(Flag::FourthTop_qq2tPrimeq_t_W_)) {
    sigmaPtr = new Sigma2qq2QqtW(8, 823);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (tPrimes || settings.get(Flag::FourthTop_ffbar2tPrimetPrimebar_s_gmZ_)) {
    sigmaPtr = new Sigma2ffbar2FFbarsgmZ(8, 824);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (tPrimes || settings.get(Flag::FourthTop_ffbar2tPrimeqbar_s_W_)) {
    sigmaPtr = new Sigma2ffbar2FfbarsW(8, 0, 825);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Set up requested objects for two different fourth-generation fermions.
  if (bPrimes || tPrimes
    || settings.get(Flag::FourthPair_ffbar2tPrimebPrimebar_s_W_)) {
    sigmaPtr = new Sigma2ffbar2FfbarsW(8, 7, 841);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::FourthPair_ffbar2tauPrimenuPrimebar_s_W_)) {
    sigmaPtr = new Sigma2ffbar2FfbarsW(17, 18, 842);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Flag for global choice between SM and BSM Higgses.
  bool useBSMHiggses = settings.get(Flag::Higgs_useBSM);

  // Set up requested objects for Standard-Model Higgs production.
  if (!useBSMHiggses) {
    bool HiggsesSM = settings.get(Flag::HiggsSM_all);
    if (HiggsesSM || settings.get(Flag::HiggsSM_ffbar2H)) {
      sigmaPtr = new Sigma1ffbar2H(0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesSM || settings.get(Flag::HiggsSM_gg2H)) {
     sigmaPtr = new Sigma1gg2H(0);
     containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesSM || settings.get(Flag::HiggsSM_gmgm2H)) {
      sigmaPtr = new Sigma1gmgm2H(0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesSM || settings.get(Flag::HiggsSM_ffbar2HZ)) {
      sigmaPtr = new Sigma2ffbar2HZ(0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesSM || settings.get(Flag::HiggsSM_ffbar2HW)) {
      sigmaPtr = new Sigma2ffbar2HW(0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesSM || settings.get(Flag::HiggsSM_ff2Hff_t_ZZ_)) {
      sigmaPtr = new Sigma3ff2HfftZZ(0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesSM || settings.get(Flag::HiggsSM_ff2Hff_t_WW_)) {
      sigmaPtr = new Sigma3ff2HfftWW(0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesSM || settings.get(Flag::HiggsSM_gg2Httbar)) {
      sigmaPtr = new Sigma3gg2HQQbar(6,0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesSM || settings.get(Flag::HiggsSM_qqbar2Httbar)) {
      sigmaPtr = new Sigma3qqbar2HQQbar(6,0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }

    // Further Standard-Model Higgs processes, not included in "all".
    if (settings.get(Flag::HiggsSM_qg2Hq)) {
      sigmaPtr = new Sigma2qg2Hq(4,0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
      sigmaPtr = new Sigma2qg2Hq(5,0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.get(Flag::HiggsSM_gg2Hbbbar)) {
      sigmaPtr = new Sigma3gg2HQQbar(5,0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.get(Flag::HiggsSM_qqbar2Hbbbar)) {
      sigmaPtr = new Sigma3qqbar2HQQbar(5,0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.get(Flag::HiggsSM_gg2Hg_l_t_)) {
      sigmaPtr = new Sigma2gg2Hglt(0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.get(Flag::HiggsSM_qg2Hq_l_t_)) {
      sigmaPtr = new Sigma2qg2Hqlt(0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.get(Flag::HiggsSM_qqbar2Hg_l_t_)) {
      sigmaPtr = new Sigma2qqbar2Hglt(0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
  }

  // Common switch for the group of Higgs production BSM.
  if (useBSMHiggses) {
    bool HiggsesBSM = settings.get(Flag::HiggsBSM_all);

    // Set up requested objects for BSM H1 production.
    bool HiggsesH1 = settings.get(Flag::HiggsBSM_allH1);
    if (HiggsesBSM || HiggsesH1 || settings.get(Flag::HiggsBSM_ffbar2H1)) {
      sigmaPtr = new Sigma1ffbar2H(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH1 || settings.get(Flag::HiggsBSM_gg2H1)) {
      sigmaPtr = new Sigma1gg2H(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH1 || settings.get(Flag::HiggsBSM_gmgm2H1)) {
      sigmaPtr = new Sigma1gmgm2H(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH1 || settings.get(Flag::HiggsBSM_ffbar2H1Z)) {
      sigmaPtr = new Sigma2ffbar2HZ(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH1 || settings.get(Flag::HiggsBSM_ffbar2H1W)) {
      sigmaPtr = new Sigma2ffbar2HW(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH1 || settings.get(Flag::HiggsBSM_ff2H1ff_t_ZZ_)) {
      sigmaPtr = new Sigma3ff2HfftZZ(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH1 || settings.get(Flag::HiggsBSM_ff2H1ff_t_WW_)) {
      sigmaPtr = new Sigma3ff2HfftWW(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH1 || settings.get(Flag::HiggsBSM_gg2H1ttbar)) {
      sigmaPtr = new Sigma3gg2HQQbar(6,1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH1 || settings.get(Flag::HiggsBSM_qqbar2H1ttbar)) {
      sigmaPtr = new Sigma3qqbar2HQQbar(6,1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }

    // Further BSM H1 processes, not included in "all".
    if (settings.get(Flag::HiggsBSM_qg2H1q)) {
      sigmaPtr = new Sigma2qg2Hq(4,1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
      sigmaPtr = new Sigma2qg2Hq(5,1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.get(Flag::HiggsBSM_gg2H1bbbar)) {
      sigmaPtr = new Sigma3gg2HQQbar(5,1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.get(Flag::HiggsBSM_qqbar2H1bbbar)) {
      sigmaPtr = new Sigma3qqbar2HQQbar(5,1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.get(Flag::HiggsBSM_gg2H1g_l_t_)) {
      sigmaPtr = new Sigma2gg2Hglt(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.get(Flag::HiggsBSM_qg2H1q_l_t_)) {
      sigmaPtr = new Sigma2qg2Hqlt(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.get(Flag::HiggsBSM_qqbar2H1g_l_t_)) {
      sigmaPtr = new Sigma2qqbar2Hglt(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }

    // Set up requested objects for BSM H2 production.
    bool HiggsesH2 = settings.get(Flag::HiggsBSM_allH2);
    if (HiggsesBSM || HiggsesH2 || settings.get(Flag::HiggsBSM_ffbar2H2)) {
      sigmaPtr = new Sigma1ffbar2H(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH2 || settings.get(Flag::HiggsBSM_gg2H2)) {
      sigmaPtr = new Sigma1gg2H(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH2 || settings.get(Flag::HiggsBSM_gmgm2H2)) {
      sigmaPtr = new Sigma1gmgm2H(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH2 || settings.get(Flag::HiggsBSM_ffbar2H2Z)) {
      sigmaPtr = new Sigma2ffbar2HZ(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH2 || settings.get(Flag::HiggsBSM_ffbar2H2W)) {
      sigmaPtr = new Sigma2ffbar2HW(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH2 || settings.get(Flag::HiggsBSM_ff2H2ff_t_ZZ_)) {
      sigmaPtr = new Sigma3ff2HfftZZ(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH2 || settings.get(Flag::HiggsBSM_ff2H2ff_t_WW_)) {
      sigmaPtr = new Sigma3ff2HfftWW(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH2 || settings.get(Flag::HiggsBSM_gg2H2ttbar)) {
      sigmaPtr = new Sigma3gg2HQQbar(6,2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH2 || settings.get(Flag::HiggsBSM_qqbar2H2ttbar)) {
      sigmaPtr = new Sigma3qqbar2HQQbar(6,2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }

    // Further BSM H2 processes, not included in "all".
   if (settings.get(Flag::HiggsBSM_qg2H2q)) {
      sigmaPtr = new Sigma2qg2Hq(4,2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
      sigmaPtr = new Sigma2qg2Hq(5,2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.get(Flag::HiggsBSM_gg2H2bbbar)) {
      sigmaPtr = new Sigma3gg2HQQbar(5,2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.get(Flag::HiggsBSM_qqbar2H2bbbar)) {
      sigmaPtr = new Sigma3qqbar2HQQbar(5,2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.get(Flag::HiggsBSM_gg2H2g_l_t_)) {
      sigmaPtr = new Sigma2gg2Hglt(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.get(Flag::HiggsBSM_qg2H2q_l_t_)) {
      sigmaPtr = new Sigma2qg2Hqlt(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.get(Flag::HiggsBSM_qqbar2H2g_l_t_)) {
      sigmaPtr = new Sigma2qqbar2Hglt(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }

    // Set up requested objects for BSM A3 production.
    bool HiggsesA3 = settings.get(Flag::HiggsBSM_allA3);
    if (HiggsesBSM || HiggsesA3 || settings.get(Flag::HiggsBSM_ffbar2A3)) {
      sigmaPtr = new Sigma1ffbar2H(3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesA3 || settings.get(Flag::HiggsBSM_gg2A3)) {
      sigmaPtr = new Sigma1gg2H(3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesA3 || settings.get(Flag::HiggsBSM_gmgm2A3)) {
      sigmaPtr = new Sigma1gmgm2H(3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesA3 || settings.get(Flag::HiggsBSM_ffbar2A3Z)) {
      sigmaPtr = new Sigma2ffbar2HZ(3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesA3 || settings.get(Flag::HiggsBSM_ffbar2A3W)) {
      sigmaPtr = new Sigma2ffbar2HW(3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesA3 || settings.get(Flag::HiggsBSM_ff2A3ff_t_ZZ_)) {
      sigmaPtr = new Sigma3ff2HfftZZ(3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesA3 || settings.get(Flag::HiggsBSM_ff2A3ff_t_WW_)) {
      sigmaPtr = new Sigma3ff2HfftWW(3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesA3 || settings.get(Flag::HiggsBSM_gg2A3ttbar)) {
      sigmaPtr = new Sigma3gg2HQQbar(6,3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesA3 || settings.get(Flag::HiggsBSM_qqbar2A3ttbar)) {
      sigmaPtr = new Sigma3qqbar2HQQbar(6,3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }

    // Further BSM A3 processes, not included in "all".
    if (settings.get(Flag::HiggsBSM_qg2A3q)) {
      sigmaPtr = new Sigma2qg2Hq(4,3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
      sigmaPtr = new Sigma2qg2Hq(5,3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.get(Flag::HiggsBSM_gg2A3bbbar)) {
      sigmaPtr = new Sigma3gg2HQQbar(5,3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.get(Flag::HiggsBSM_qqbar2A3bbbar)) {
      sigmaPtr = new Sigma3qqbar2HQQbar(5,3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.get(Flag::HiggsBSM_gg2A3g_l_t_)) {
      sigmaPtr = new Sigma2gg2Hglt(3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.get(Flag::HiggsBSM_qg2A3q_l_t_)) {
      sigmaPtr = new Sigma2qg2Hqlt(3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.get(Flag::HiggsBSM_qqbar2A3g_l_t_)) {
      sigmaPtr = new Sigma2qqbar2Hglt(3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }

    // Set up requested objects for Charged Higgs production
    bool HiggsesChg = settings.get(Flag::HiggsBSM_allHplusminus);
    if (HiggsesBSM || HiggsesChg || settings.get(Flag::HiggsBSM_ffbar2Hplusminus)) {
      sigmaPtr = new Sigma1ffbar2Hchg;
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesChg || settings.get(Flag::HiggsBSM_bg2Hplusminust)) {
      sigmaPtr = new Sigma2qg2Hchgq(6, 1062, "b g -> H+- t");
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }

    // Set up requested objects for Higgs pair-production
    bool HiggsesPairs = settings.get(Flag::HiggsBSM_allHpair);
    if (HiggsesBSM || HiggsesPairs || settings.get(Flag::HiggsBSM_ffbar2A3H1)) {
      sigmaPtr = new Sigma2ffbar2A3H12(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesPairs || settings.get(Flag::HiggsBSM_ffbar2A3H2)) {
      sigmaPtr = new Sigma2ffbar2A3H12(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesPairs || settings.get(Flag::HiggsBSM_ffbar2HplusminusH1)) {
      sigmaPtr = new Sigma2ffbar2HchgH12(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesPairs || settings.get(Flag::HiggsBSM_ffbar2HplusminusH2)) {
      sigmaPtr = new Sigma2ffbar2HchgH12(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesPairs || settings.get(Flag::HiggsBSM_ffbar2HplusHminus)) {
      sigmaPtr = new Sigma2ffbar2HposHneg();
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
  }

  // Set up requested objects for SUSY pair processes.
  if (couplings->isSUSY) {
    CoupSUSY* coupSUSY = (CoupSUSY *) couplings;

    bool SUSYs = settings.get(Flag::SUSY_all);
    bool nmssm = settings.get(Flag::SLHA_NMSSM);

    // Preselected SUSY codes.
    setupIdVecs( settings);

    // MSSM: 4 neutralinos
    int nNeut = 4;
    if (nmssm) nNeut = 5;

    // Gluino-gluino
    if (SUSYs || settings.get(Flag::SUSY_gg2gluinogluino)) {
      // Skip if outgoing codes not asked for
      if (allowIdVals( 1000021, 1000021)) {
        sigmaPtr = new Sigma2gg2gluinogluino();
        containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
      }
    }
    if (SUSYs || settings.get(Flag::SUSY_qqbar2gluinogluino)) {
      // Skip if outgoing codes not asked for
      if (allowIdVals( 1000021, 1000021)) {
        sigmaPtr = new Sigma2qqbar2gluinogluino();
        containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
      }
    }

    // Gluino-squark
    if (SUSYs || settings.get(Flag::SUSY_qg2squarkgluino)) {
      int iproc = 1202;
      for (int idx = 1; idx <= 6; ++idx) {
        for (int iso = 1; iso <= 2; ++iso) {
          iproc++;
          int id3 = iso + ((idx <= 3)
                  ? 1000000+2*(idx-1) : 2000000+2*(idx-4));
          int id4 = 1000021;
          // Skip if outgoing codes not asked for
          if (!allowIdVals( id3, id4)) continue;
          sigmaPtr = new Sigma2qg2squarkgluino(id3,iproc);
          containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
        }
      }
    }

    // Squark-antisquark (gg initiated)
    if (SUSYs || settings.get(Flag::SUSY_gg2squarkantisquark)) {
      int iproc = 1214;
      for (int idx = 1; idx <= 6; ++idx) {
        for (int iso = 1; iso <= 2; ++iso) {
          iproc++;
          int id = iso + ((idx <= 3)
                 ? 1000000+2*(idx-1) : 2000000+2*(idx-4));
          // Skip if outgoing codes not asked for
          if (!allowIdVals( id, id)) continue;
          sigmaPtr = new Sigma2gg2squarkantisquark(id,iproc);
          containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
        }
      }
    }

    // Squark-antisquark (qqbar initiated)
    if (SUSYs || settings.get(Flag::SUSY_qqbar2squarkantisquark)) {
      int iproc = 1230;
      for (int idx = 1; idx <= 6; ++idx) {
        for (int iso = 1; iso <= 2; ++iso) {
          for (int jso = iso; jso >= 1; --jso) {
            for (int jdx = 1; jdx <= 6; ++jdx) {
              if (iso == jso && jdx < idx) continue;
              int id1 = iso + ((idx <= 3) ? 1000000+2*(idx-1)
                               : 2000000+2*(idx-4));
              int id2 = jso + ((jdx <= 3) ? 1000000+2*(jdx-1)
                               : 2000000+2*(jdx-4));
              // Update process number counter (for ~q~q, +2 if not self-conj)
              //if (iproc == 1302) iproc=1310;
              iproc++;
              if (iso == jso && id1 != id2) iproc++;
              // Skip if outgoing codes not asked for
              if (!allowIdVals( id1, id2)) continue;
              if (iso == jso && id1 != id2) {
                sigmaPtr = new Sigma2qqbar2squarkantisquark(id1,-id2,iproc-1);
                containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
                sigmaPtr = new Sigma2qqbar2squarkantisquark(id2,-id1,iproc);
                containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
              } else {
                sigmaPtr = new Sigma2qqbar2squarkantisquark(id1,-id2,iproc);
                containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
              }
            }
          }
        }
      }
    }

    // Squark-squark
    if (SUSYs || settings.get(Flag::SUSY_qq2squarksquark)) {
      int iproc = 1350;
      for (int idx = 1; idx <= 6; ++idx) {
        for (int iso = 1; iso <= 2; ++iso) {
          for (int jso = iso; jso >= 1; jso--) {
            for (int jdx = 1; jdx <= 6; ++jdx) {
              if (iso == jso && jdx < idx) continue;
              iproc++;
              int id1 = iso + ((idx <= 3)
                      ? 1000000+2*(idx-1) : 2000000+2*(idx-4));
              int id2 = jso + ((jdx <= 3)
                      ? 1000000+2*(jdx-1) : 2000000+2*(jdx-4));
              // Skip if outgoing codes not asked for
              if (!allowIdVals( id1, id2)) continue;
              sigmaPtr = new Sigma2qq2squarksquark(id1,id2,iproc);
              containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
            }
          }
        }
      }
    }

    // Neutralino + squark
    if (SUSYs || settings.get(Flag::SUSY_qg2chi0squark)) {
      int iproc = 1430;
      for (int iNeut= 1; iNeut <= nNeut; iNeut++) {
        for (int idx = 1; idx <= 6; idx++) {
          bool isUp = false;
          for (int iso = 1; iso <= 2; iso++) {
            if (iso == 2) isUp = true;
            iproc++;
            int id3 = coupSUSY->idNeut(iNeut);
            int id4 = iso + ((idx <= 3)
                    ? 1000000+2*(idx-1) : 2000000+2*(idx-4));
            // Skip if outgoing codes not asked for
            if (!allowIdVals( id3, id4)) continue;
            sigmaPtr = new Sigma2qg2chi0squark(iNeut,idx,isUp,iproc);
            containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
          }
        }
      }
    }

    // Chargino + squark
    if (SUSYs || settings.get(Flag::SUSY_qg2chiplusminussquark)) {
      int iproc = 1490;
      for (int iChar = 1; iChar <= 2; iChar++) {
        for (int idx = 1; idx <= 6; idx++) {
          bool isUp = false;
          for (int iso = 1; iso <= 2; iso++) {
            if (iso == 2) isUp = true;
            iproc++;
            int id3 = coupSUSY->idChar(iChar);
            int id4 = iso + ((idx <= 3)
                    ? 1000000+2*(idx-1) : 2000000+2*(idx-4));
            // Skip if outgoing codes not asked for
            if (!allowIdVals( id3, id4)) continue;
            sigmaPtr = new Sigma2qg2charsquark(iChar,idx,isUp,iproc);
            containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
          }
        }
      }
    }

    // Neutralino pairs
    if (SUSYs || settings.get(Flag::SUSY_qqbar2chi0chi0)) {
      int iproc = 1550;
      for (int iNeut2 = 1; iNeut2 <= nNeut; iNeut2++) {
        for (int iNeut1 = 1; iNeut1 <= iNeut2; iNeut1++) {
          iproc++;
          // Skip if outgoing codes not asked for
          if (!allowIdVals( coupSUSY->idNeut(iNeut1),
            coupSUSY->idNeut(iNeut2) ) ) continue;
          sigmaPtr = new Sigma2qqbar2chi0chi0(iNeut1, iNeut2,iproc);
          containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
        }
      }
    }

    // Neutralino-Chargino
    if (SUSYs || settings.get(Flag::SUSY_qqbar2chiplusminuschi0)) {
      int iproc = 1570;
      for (int iNeut = 1; iNeut <= nNeut; iNeut++) {
        for (int iChar = 1; iChar <= 2; ++iChar) {
          iproc += 2;
          // Skip if outgoing codes not asked for
          if (!allowIdVals( coupSUSY->idNeut(iNeut),
            coupSUSY->idChar(iChar) ) ) continue;
          sigmaPtr = new Sigma2qqbar2charchi0( iChar, iNeut, iproc-1);
          containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
          sigmaPtr = new Sigma2qqbar2charchi0(-iChar, iNeut, iproc);
          containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
        }
      }
    }

    // Chargino-Chargino
    if (SUSYs || settings.get(Flag::SUSY_qqbar2chipluschiminus)) {
      int iproc = 1590;
      for (int i = 1; i <= 2; ++i) {
        for (int j = 1; j <= 2; ++j) {
          iproc++;
          // Skip if outgoing codes not asked for
          if (!allowIdVals( coupSUSY->idChar(i),
            coupSUSY->idChar(j) ) ) continue;
          sigmaPtr = new Sigma2qqbar2charchar( i,-j, iproc);
          containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
        }
      }
    }

    // RPV squark production
    if(SUSYs || settings.get(Flag::SUSY_qq2antisquark)) {
      for (int idx = 1; idx <= 6; ++idx) {
        for (int iso = 1; iso <= 2; ++iso) {
          int id1 = iso + ((idx <= 3) ? 1000000+2*(idx-1) : 2000000+2*(idx-4));
          // Skip if outgoing code not asked for
          if (!allowIdVals( id1, 0)) continue;
          sigmaPtr = new Sigma1qq2antisquark(id1);
          containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
        }
      }
    }

    // Neutralino-gluino
    if (SUSYs || settings.get(Flag::SUSY_qqbar2chi0gluino)) {
      int iproc = 1600;
      for (int iNeut = 1; iNeut <= nNeut; iNeut++) {
        iproc++;
        // Skip if outgoing codes not asked for
        if (!allowIdVals( coupSUSY->idNeut(iNeut), 1000021)) continue;
        sigmaPtr = new Sigma2qqbar2chi0gluino(iNeut, iproc);
        containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
      }
    }

    // Chargino-Gluino
    if (SUSYs || settings.get(Flag::SUSY_qqbar2chiplusminusgluino)) {
      int iproc = 1620;
      for (int iChar = 1; iChar <= 2; ++iChar) {
        iproc ++;
        // Skip if outgoing codes not asked for
        if (!allowIdVals( coupSUSY->idChar(iChar), 1000021)) continue;
        sigmaPtr = new Sigma2qqbar2chargluino( iChar, iproc);
        containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
      }
    }

    // Slepton-antislepton (qqbar initiated); Currently no RH sneutrinos
    if (SUSYs || settings.get(Flag::SUSY_qqbar2sleptonantislepton)) {
      int iproc = 1650;
      for (int idx = 1; idx <= 6; ++idx) {
        for (int iso = 1; iso <= 2; ++iso) {
          for (int jso = iso; jso >= 1; --jso) {
            for (int jdx = 1; jdx <= 6; ++jdx) {
              if (iso == jso && jdx < idx) continue;
              int id1 = iso + ((idx <= 3) ? 1000010+2*(idx-1)
                               : 2000010+2*(idx-4));
              int id2 = jso + ((jdx <= 3) ? 1000010+2*(jdx-1)
                               : 2000010+2*(jdx-4));
              // Update process number counter
              iproc++;
              if (iso == jso && id1 != id2) iproc++;
              // Exclude RH neutrinos from allowed final states
              if (abs(id1) >= 2000012 && id1 % 2 == 0) continue;
              if (abs(id2) >= 2000012 && id2 % 2 == 0) continue;
              // Skip if outgoing codes not asked for
              if (!allowIdVals( id1, id2)) continue;
              if (iso == jso && id1 != id2) {
                sigmaPtr
                  = new Sigma2qqbar2sleptonantislepton(id1,-id2,iproc-1);
                containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
                sigmaPtr = new Sigma2qqbar2sleptonantislepton(id2,-id1,iproc);
                containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
              } else {
                sigmaPtr = new Sigma2qqbar2sleptonantislepton(id1,-id2,iproc);
                containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
              }
            }
          }
        }
      }
    }

  } // End of SUSY processes.

  // // Set up requested objects for New-Gauge-Boson processes.
  // if (settings.get(Flag::NewGaugeBoson_ffbar2gmZZprime)) {
  //   sigmaPtr = new Sigma1ffbar2gmZZprime();
  //   containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  // }
  // if (settings.get(Flag::NewGaugeBoson_ffbar2Wprime)) {
  //   sigmaPtr = new Sigma1ffbar2Wprime();
  //   containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  // }
  // if (settings.get(Flag::NewGaugeBoson_ffbar2R0)) {
  //   sigmaPtr = new Sigma1ffbar2Rhorizontal();
  //   containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  // }

  // // Set up requested objects for Left-Right-Symmetry processes.
  // bool leftrights = settings.get(Flag::LeftRightSymmmetry_all);
  // if (leftrights || settings.get(Flag::LeftRightSymmmetry_ffbar2ZR)) {
  //   sigmaPtr = new Sigma1ffbar2ZRight();
  //   containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  // }
  // if (leftrights || settings.get(Flag::LeftRightSymmmetry_ffbar2WR)) {
  //   sigmaPtr = new Sigma1ffbar2WRight();
  //   containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  // }
  // if (leftrights || settings.get(Flag::LeftRightSymmmetry_ll2HL)) {
  //   sigmaPtr = new Sigma1ll2Hchgchg(1);
  //   containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  // }
  // if (leftrights || settings.get(Flag::LeftRightSymmmetry_lgm2HLe)) {
  //   sigmaPtr = new Sigma2lgm2Hchgchgl(1, 11);
  //   containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  // }
  // if (leftrights || settings.get(Flag::LeftRightSymmmetry_lgm2HLmu)) {
  //   sigmaPtr = new Sigma2lgm2Hchgchgl(1, 13);
  //   containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  // }
  // if (leftrights || settings.get(Flag::LeftRightSymmmetry_lgm2HLtau)) {
  //   sigmaPtr = new Sigma2lgm2Hchgchgl(1, 15);
  //   containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  // }
  // if (leftrights || settings.get(Flag::LeftRightSymmmetry_ff2HLff)) {
  //   sigmaPtr = new Sigma3ff2HchgchgfftWW(1);
  //   containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  // }
  // if (leftrights || settings.get(Flag::LeftRightSymmmetry_ffbar2HLHL)) {
  //   sigmaPtr = new Sigma2ffbar2HchgchgHchgchg(1);
  //   containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  // }
  // if (leftrights || settings.get(Flag::LeftRightSymmmetry_ll2HR)) {
  //   sigmaPtr = new Sigma1ll2Hchgchg(2);
  //   containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  // }
  // if (leftrights || settings.get(Flag::LeftRightSymmmetry_lgm2HRe)) {
  //   sigmaPtr = new Sigma2lgm2Hchgchgl(2, 11);
  //   containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  // }
  // if (leftrights || settings.get(Flag::LeftRightSymmmetry_lgm2HRmu)) {
  //   sigmaPtr = new Sigma2lgm2Hchgchgl(2, 13);
  //   containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  // }
  // if (leftrights || settings.get(Flag::LeftRightSymmmetry_lgm2HRtau)) {
  //   sigmaPtr = new Sigma2lgm2Hchgchgl(2, 15);
  //   containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  // }
  // if (leftrights || settings.get(Flag::LeftRightSymmmetry_ff2HRff)) {
  //   sigmaPtr = new Sigma3ff2HchgchgfftWW(2);
  //   containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  // }
  // if (leftrights || settings.get(Flag::LeftRightSymmmetry_ffbar2HRHR)) {
  //   sigmaPtr = new Sigma2ffbar2HchgchgHchgchg(2);
  //   containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  // }

  // // Set up requested objects for leptoquark LQ processes.
  // bool leptoquarks = settings.get(Flag::LeptoQuark_all);
  // if (leptoquarks || settings.get(Flag::LeptoQuark_ql2LQ)) {
  //   sigmaPtr = new Sigma1ql2LeptoQuark;
  //   containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  // }
  // if (leptoquarks || settings.get(Flag::LeptoQuark_qg2LQl)) {
  //   sigmaPtr = new Sigma2qg2LeptoQuarkl;
  //   containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  // }
  // if (leptoquarks || settings.get(Flag::LeptoQuark_gg2LQLQbar)) {
  //   sigmaPtr = new Sigma2gg2LQLQbar;
  //   containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  // }
  // if (leptoquarks || settings.get(Flag::LeptoQuark_qqbar2LQLQbar)) {
  //   sigmaPtr = new Sigma2qqbar2LQLQbar;
  //   containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  // }

  // Set up requested objects for excited-fermion processes.
  bool excitedfermions = settings.get(Flag::ExcitedFermion_all);
  if (excitedfermions || settings.get(Flag::ExcitedFermion_dg2dStar)) {
    sigmaPtr = new Sigma1qg2qStar(1);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions || settings.get(Flag::ExcitedFermion_ug2uStar)) {
    sigmaPtr = new Sigma1qg2qStar(2);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions || settings.get(Flag::ExcitedFermion_sg2sStar)) {
    sigmaPtr = new Sigma1qg2qStar(3);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions || settings.get(Flag::ExcitedFermion_cg2cStar)) {
    sigmaPtr = new Sigma1qg2qStar(4);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions || settings.get(Flag::ExcitedFermion_bg2bStar)) {
    sigmaPtr = new Sigma1qg2qStar(5);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions || settings.get(Flag::ExcitedFermion_egm2eStar)) {
    sigmaPtr = new Sigma1lgm2lStar(11);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions || settings.get(Flag::ExcitedFermion_mugm2muStar)) {
    sigmaPtr = new Sigma1lgm2lStar(13);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions || settings.get(Flag::ExcitedFermion_taugm2tauStar)) {
    sigmaPtr = new Sigma1lgm2lStar(15);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions || settings.get(Flag::ExcitedFermion_qq2dStarq)) {
    sigmaPtr = new Sigma2qq2qStarq(1);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions || settings.get(Flag::ExcitedFermion_qq2uStarq)) {
    sigmaPtr = new Sigma2qq2qStarq(2);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions || settings.get(Flag::ExcitedFermion_qq2sStarq)) {
    sigmaPtr = new Sigma2qq2qStarq(3);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions || settings.get(Flag::ExcitedFermion_qq2cStarq)) {
    sigmaPtr = new Sigma2qq2qStarq(4);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions || settings.get(Flag::ExcitedFermion_qq2bStarq)) {
    sigmaPtr = new Sigma2qq2qStarq(5);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions || settings.get(Flag::ExcitedFermion_qqbar2eStare)) {
    sigmaPtr = new Sigma2qqbar2lStarlbar(11);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions || settings.get(Flag::ExcitedFermion_qqbar2nueStarnue)) {
    sigmaPtr = new Sigma2qqbar2lStarlbar(12);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions || settings.get(Flag::ExcitedFermion_qqbar2muStarmu)) {
    sigmaPtr = new Sigma2qqbar2lStarlbar(13);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions || settings.get(Flag::ExcitedFermion_qqbar2numuStarnumu)) {
    sigmaPtr = new Sigma2qqbar2lStarlbar(14);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions || settings.get(Flag::ExcitedFermion_qqbar2tauStartau)) {
    sigmaPtr = new Sigma2qqbar2lStarlbar(15);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions
    || settings.get(Flag::ExcitedFermion_qqbar2nutauStarnutau)) {
    sigmaPtr = new Sigma2qqbar2lStarlbar(16);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions || settings.get(Flag::ExcitedFermion_qqbar2eStareStar)) {
    sigmaPtr = new Sigma2qqbar2lStarlStarBar(11);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions
    || settings.get(Flag::ExcitedFermion_qqbar2nueStarnueStar)) {
    sigmaPtr = new Sigma2qqbar2lStarlStarBar(12);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions || settings.get(Flag::ExcitedFermion_qqbar2muStarmuStar)) {
    sigmaPtr = new Sigma2qqbar2lStarlStarBar(13);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions
    || settings.get(Flag::ExcitedFermion_qqbar2numuStarnumuStar)) {
    sigmaPtr = new Sigma2qqbar2lStarlStarBar(14);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions
    || settings.get(Flag::ExcitedFermion_qqbar2tauStartauStar)) {
    sigmaPtr = new Sigma2qqbar2lStarlStarBar(15);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (excitedfermions
    || settings.get(Flag::ExcitedFermion_qqbar2nutauStarnutauStar)) {
    sigmaPtr = new Sigma2qqbar2lStarlStarBar(16);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Set up requested objects for contact interaction processes.
  if (settings.get(Flag::ContactInteractions_QCqq2qq)) {
    sigmaPtr = new Sigma2QCqq2qq();
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ContactInteractions_QCqqbar2qqbar)) {
    sigmaPtr = new Sigma2QCqqbar2qqbar();
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ContactInteractions_QCffbar2eebar)) {
    sigmaPtr = new Sigma2QCffbar2llbar(11, 4203);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ContactInteractions_QCffbar2mumubar)) {
    sigmaPtr = new Sigma2QCffbar2llbar(13, 4204);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ContactInteractions_QCffbar2tautaubar)) {
    sigmaPtr = new Sigma2QCffbar2llbar(15, 4205);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Set spin of particles in the Hidden Valley scenario.
  int spinFv = settings.get(Mode::HiddenValley_spinFv);
  for (int i = 1; i < 7; ++i) {
    if (particleDataPtr->spinType( 4900000 + i) != spinFv + 1)
        particleDataPtr->spinType( 4900000 + i,    spinFv + 1);
    if (particleDataPtr->spinType( 4900010 + i) != spinFv + 1)
        particleDataPtr->spinType( 4900010 + i,    spinFv + 1);
  }
  if (spinFv != 1) {
    if (particleDataPtr->spinType( 4900101) != 2)
       particleDataPtr->spinType( 4900101, 2);
  } else {
    int spinqv = settings.get(Mode::HiddenValley_spinqv);
    if (particleDataPtr->spinType( 4900101) != 2 * spinqv + 1)
        particleDataPtr->spinType( 4900101,    2 * spinqv + 1);
  }

  // Set up requested objects for HiddenValley processes.
  bool hiddenvalleys = settings.get(Flag::HiddenValley_all);
  if (hiddenvalleys || settings.get(Flag::HiddenValley_gg2DvDvbar)) {
    sigmaPtr = new Sigma2gg2qGqGbar( 4900001, 4901, spinFv,
      "g g -> Dv Dvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_gg2UvUvbar)) {
    sigmaPtr = new Sigma2gg2qGqGbar( 4900002, 4902, spinFv,
      "g g -> Uv Uvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_gg2SvSvbar)) {
    sigmaPtr = new Sigma2gg2qGqGbar( 4900003, 4903, spinFv,
      "g g -> Sv Svbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_gg2CvCvbar)) {
    sigmaPtr = new Sigma2gg2qGqGbar( 4900004, 4904, spinFv,
      "g g -> Cv Cvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_gg2BvBvbar)) {
    sigmaPtr = new Sigma2gg2qGqGbar( 4900005, 4905, spinFv,
      "g g -> Bv Bvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_gg2TvTvbar)) {
    sigmaPtr = new Sigma2gg2qGqGbar( 4900006, 4906, spinFv,
      "g g -> Tv Tvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_qqbar2DvDvbar)) {
    sigmaPtr = new Sigma2qqbar2qGqGbar( 4900001, 4911, spinFv,
      "q qbar -> Dv Dvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_qqbar2UvUvbar)) {
    sigmaPtr = new Sigma2qqbar2qGqGbar( 4900002, 4912, spinFv,
      "q qbar -> Uv Uvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_qqbar2SvSvbar)) {
    sigmaPtr = new Sigma2qqbar2qGqGbar( 4900003, 4913, spinFv,
      "q qbar -> Sv Svbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_qqbar2CvCvbar)) {
    sigmaPtr = new Sigma2qqbar2qGqGbar( 4900004, 4914, spinFv,
      "q qbar -> Cv Cvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_qqbar2BvBvbar)) {
    sigmaPtr = new Sigma2qqbar2qGqGbar( 4900005, 4915, spinFv,
      "q qbar -> Bv Bvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_qqbar2TvTvbar)) {
    sigmaPtr = new Sigma2qqbar2qGqGbar( 4900006, 4916, spinFv,
      "q qbar -> Tv Tvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_ffbar2DvDvbar)) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900001, 4921, spinFv,
      "f fbar -> Dv Dvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_ffbar2UvUvbar)) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900002, 4922, spinFv,
      "f fbar -> Uv Uvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_ffbar2SvSvbar)) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900003, 4923, spinFv,
      "f fbar -> Sv Svbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_ffbar2CvCvbar)) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900004, 4924, spinFv,
      "f fbar -> Cv Cvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_ffbar2BvBvbar)) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900005, 4925, spinFv,
      "f fbar -> Bv Bvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_ffbar2TvTvbar)) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900006, 4926, spinFv,
      "f fbar -> Tv Tvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_ffbar2EvEvbar)) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900011, 4931, spinFv,
      "f fbar -> Ev Evbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_ffbar2nuEvnuEvbar)) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900012, 4932, spinFv,
      "f fbar -> nuEv nuEvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_ffbar2MUvMUvbar)) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900013, 4933, spinFv,
      "f fbar -> MUv MUvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_ffbar2nuMUvnuMUvbar)) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900014, 4934, spinFv,
      "f fbar -> nuMUv nuMUvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_ffbar2TAUvTAUvbar)) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900015, 4935, spinFv,
      "f fbar -> TAUv TAUvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_ffbar2nuTAUvnuTAUvbar)) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900016, 4936, spinFv,
      "f fbar -> nuTAUv nuTAUvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hiddenvalleys || settings.get(Flag::HiddenValley_ffbar2Zv)) {
    sigmaPtr = new Sigma1ffbar2Zv();
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Set up requested objects for RS extra-dimensional G* processes.
  bool extraDimGstars = settings.get(Flag::ExtraDimensionsGstar_all);
  if (extraDimGstars || settings.get(Flag::ExtraDimensionsGstar_gg2Gstar)) {
    sigmaPtr = new Sigma1gg2GravitonStar;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (extraDimGstars || settings.get(Flag::ExtraDimensionsGstar_ffbar2Gstar)) {
    sigmaPtr = new Sigma1ffbar2GravitonStar;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsGstar_gg2Gstarg)) {
    sigmaPtr = new Sigma2gg2GravitonStarg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsGstar_qg2Gstarq)) {
    sigmaPtr = new Sigma2qg2GravitonStarq;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsGstar_qqbar2Gstarg)) {
    sigmaPtr = new Sigma2qqbar2GravitonStarg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  //  Set up requested objects for RS extra-dimensional KKgluon processes.
  if (settings.get(Flag::ExtraDimensionsGstar_qqbar2KKgluonstar)) {
    sigmaPtr = new Sigma1qqbar2KKgluonStar;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // NOAM: Set up requested objects for TEV extra-dimensional processes.
  if (settings.get(Flag::ExtraDimensionsTEV_ffbar2ddbar)) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(1, 5061);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsTEV_ffbar2uubar)) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(2, 5062);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsTEV_ffbar2ssbar)) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(3, 5063);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsTEV_ffbar2ccbar)) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(4, 5064);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsTEV_ffbar2bbbar)) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(5, 5065);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsTEV_ffbar2ttbar)) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(6, 5066);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsTEV_ffbar2epluseminus)) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(11, 5071);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsTEV_ffbar2nuenuebar)) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(12, 5072);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsTEV_ffbar2muplusmuminus)) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(13, 5073);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsTEV_ffbar2numunumubar)) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(14, 5074);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsTEV_ffbar2tauplustauminus)) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(15, 5075);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsTEV_ffbar2nutaunutaubar)) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(16, 5076);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Set up requested objects for large extra-dimensional G processes.
  bool extraDimLEDmono = settings.get(Flag::ExtraDimensionsLED_monojet);
  if (extraDimLEDmono || settings.get(Flag::ExtraDimensionsLED_gg2Gg)) {
    sigmaPtr = new Sigma2gg2LEDUnparticleg( true );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (extraDimLEDmono || settings.get(Flag::ExtraDimensionsLED_qg2Gq)) {
    sigmaPtr = new Sigma2qg2LEDUnparticleq( true );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (extraDimLEDmono || settings.get(Flag::ExtraDimensionsLED_qqbar2Gg)) {
    sigmaPtr = new Sigma2qqbar2LEDUnparticleg( true );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsLED_ffbar2GZ)) {
    sigmaPtr = new Sigma2ffbar2LEDUnparticleZ( true );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsLED_ffbar2Ggamma)) {
    sigmaPtr = new Sigma2ffbar2LEDUnparticlegamma( true );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsLED_ffbar2gammagamma)) {
    sigmaPtr = new Sigma2ffbar2LEDgammagamma( true );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsLED_gg2gammagamma)) {
    sigmaPtr = new Sigma2gg2LEDgammagamma( true );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsLED_ffbar2llbar)) {
    sigmaPtr = new Sigma2ffbar2LEDllbar( true );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsLED_gg2llbar)) {
    sigmaPtr = new Sigma2gg2LEDllbar( true );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  bool extraDimLEDdij = settings.get(Flag::ExtraDimensionsLED_dijets);
  if (extraDimLEDdij || settings.get(Flag::ExtraDimensionsLED_gg2DJgg)) {
    sigmaPtr = new Sigma2gg2LEDgg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (extraDimLEDdij || settings.get(Flag::ExtraDimensionsLED_gg2DJqqbar)) {
    sigmaPtr = new Sigma2gg2LEDqqbar;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (extraDimLEDdij || settings.get(Flag::ExtraDimensionsLED_qg2DJqg)) {
    sigmaPtr = new Sigma2qg2LEDqg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (extraDimLEDdij || settings.get(Flag::ExtraDimensionsLED_qq2DJqq)) {
    sigmaPtr = new Sigma2qq2LEDqq;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (extraDimLEDdij || settings.get(Flag::ExtraDimensionsLED_qqbar2DJgg)) {
    sigmaPtr = new Sigma2qqbar2LEDgg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (extraDimLEDdij || settings.get(Flag::ExtraDimensionsLED_qqbar2DJqqbarNew)) {
    sigmaPtr = new Sigma2qqbar2LEDqqbarNew;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Set up requested objects for unparticle processes.
  bool extraDimUnpartmono = settings.get(Flag::ExtraDimensionsUnpart_monojet);
  if (extraDimUnpartmono || settings.get(Flag::ExtraDimensionsUnpart_gg2Ug)) {
    sigmaPtr = new Sigma2gg2LEDUnparticleg( false );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (extraDimUnpartmono || settings.get(Flag::ExtraDimensionsUnpart_qg2Uq)) {
    sigmaPtr = new Sigma2qg2LEDUnparticleq( false );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (extraDimUnpartmono || settings.get(Flag::ExtraDimensionsUnpart_qqbar2Ug)) {
    sigmaPtr = new Sigma2qqbar2LEDUnparticleg( false );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsUnpart_ffbar2UZ)) {
    sigmaPtr = new Sigma2ffbar2LEDUnparticleZ( false );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsUnpart_ffbar2Ugamma)) {
    sigmaPtr = new Sigma2ffbar2LEDUnparticlegamma( false );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsUnpart_ffbar2gammagamma)) {
    sigmaPtr = new Sigma2ffbar2LEDgammagamma( false );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsUnpart_gg2gammagamma)) {
    sigmaPtr = new Sigma2gg2LEDgammagamma( false );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsUnpart_ffbar2llbar)) {
    sigmaPtr = new Sigma2ffbar2LEDllbar( false );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.get(Flag::ExtraDimensionsUnpart_gg2llbar)) {
    sigmaPtr = new Sigma2gg2LEDllbar( false );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Routine to initialize list of second hard processes.

bool SetupContainers::init2(vector<ProcessContainer*>& container2Ptrs,
  Settings& settings) {

  // Reset process list, if filled in previous subrun.
  if (container2Ptrs.size() > 0) {
    for (int i = 0; i < int(container2Ptrs.size()); ++i)
      delete container2Ptrs[i];
    container2Ptrs.clear();
  }
  SigmaProcess* sigmaPtr;

  // Two hard QCD jets.
  if (settings.get(Flag::SecondHard_TwoJets)) {
    sigmaPtr = new Sigma2gg2gg;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2qqbar;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qg2qg;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qq2qq;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2gg;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2qqbarNew;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2QQbar(4, 121);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbar(4, 122);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2QQbar(5, 123);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbar(5, 124);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // A prompt photon and a hard jet.
  if (settings.get(Flag::SecondHard_PhotonAndJet)) {
    sigmaPtr = new Sigma2qg2qgamma;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2ggamma;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2ggamma;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Two prompt photons.
  if (settings.get(Flag::SecondHard_TwoPhotons)) {
    sigmaPtr = new Sigma2ffbar2gammagamma;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2gammagamma;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Charmonium.
  if (settings.get(Flag::SecondHard_Charmonium)) {
    vector<SigmaProcess*> charmoniumSigmaPtrs;
    charmonium.setupSigma2gg(charmoniumSigmaPtrs, true);
    charmonium.setupSigma2qg(charmoniumSigmaPtrs, true);
    charmonium.setupSigma2qq(charmoniumSigmaPtrs, true);
    for (unsigned int i = 0; i < charmoniumSigmaPtrs.size(); ++i)
      container2Ptrs.push_back( new ProcessContainer(charmoniumSigmaPtrs[i]));
  }

  // Bottomonium.
  if (settings.get(Flag::SecondHard_Bottomonium)) {
    vector<SigmaProcess*> bottomoniumSigmaPtrs;
    bottomonium.setupSigma2gg(bottomoniumSigmaPtrs, true);
    bottomonium.setupSigma2qg(bottomoniumSigmaPtrs, true);
    bottomonium.setupSigma2qq(bottomoniumSigmaPtrs, true);
    for (unsigned int i = 0; i < bottomoniumSigmaPtrs.size(); ++i)
      container2Ptrs.push_back( new ProcessContainer(bottomoniumSigmaPtrs[i]));
  }

  // A single gamma*/Z0.
  if (settings.get(Flag::SecondHard_SingleGmZ)) {
    sigmaPtr = new Sigma1ffbar2gmZ;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // A single W+-.
  if (settings.get(Flag::SecondHard_SingleW)) {
    sigmaPtr = new Sigma1ffbar2W;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // A gamma*/Z0 and a hard jet.
  if (settings.get(Flag::SecondHard_GmZAndJet)) {
    sigmaPtr = new Sigma2qqbar2gmZg;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qg2gmZq;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // A W+- and a hard jet.
  if (settings.get(Flag::SecondHard_WAndJet)) {
    sigmaPtr = new Sigma2qqbar2Wg;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qg2Wq;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Top pair production.
  if (settings.get(Flag::SecondHard_TopPair)) {
    sigmaPtr = new Sigma2gg2QQbar(6, 601);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbar(6, 602);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2ffbar2FFbarsgmZ(6, 604);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Single top production.
  if (settings.get(Flag::SecondHard_SingleTop)) {
    sigmaPtr = new Sigma2qq2QqtW(6, 603);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2ffbar2FfbarsW(6, 0, 605);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Two b jets - already part of TwoJets sample above.
  if (settings.get(Flag::SecondHard_TwoBJets)) {
    sigmaPtr = new Sigma2gg2QQbar(5, 123);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbar(5, 124);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Set up arrays of allowed outgoing SUSY particles.

void SetupContainers::setupIdVecs( Settings& settings) {

  // First array either none, one or many particles.
  idVecA.resize(0);
  if (settings.get(Mode::SUSY_idA) != 0) {
    idVecA.push_back( abs(settings.get(Mode::SUSY_idA)) );
  } else {
    vector<int> idTmpA = settings.get(ModeList::SUSY_idVecA);
    for (int i = 0; i < int(idTmpA.size()); ++i)
      if (idTmpA[i] != 0) idVecA.push_back( abs(idTmpA[i]) );
  }
  nVecA = idVecA.size();

  // Second array either none, one or many particles.
  idVecB.resize(0);
  if (settings.get(Mode::SUSY_idB) != 0) {
    idVecB.push_back( abs(settings.get(Mode::SUSY_idB)) );
  } else {
    vector<int> idTmpB = settings.get(ModeList::SUSY_idVecB);
    for (int i = 0; i < int(idTmpB.size()); ++i)
    if (idTmpB[i] != 0) idVecB.push_back( abs(idTmpB[i]) );
  }
  nVecB = idVecB.size();

}

//--------------------------------------------------------------------------

// Check final state for allowed outgoing SUSY particles.
// Normally check two codes, but allow for only one.

bool SetupContainers::allowIdVals( int idCheck1, int idCheck2) {

  // If empty arrays or id's no need for checks. Else need absolute values.
  if (nVecA == 0 && nVecB == 0) return true;
  if (idCheck1 == 0 && idCheck2 == 0) return true;
  int idChk1 = abs(idCheck1);
  int idChk2 = abs(idCheck2);

  // If only one outgoing particle then check idVecA and idVecB.
  if (idChk1 == 0) swap(idChk1, idChk2);
  if (idChk2 == 0) {
    for (int i = 0; i < nVecA; ++i) if (idChk1 == idVecA[i]) return true;
    for (int i = 0; i < nVecB; ++i) if (idChk1 == idVecB[i]) return true;
    return false;
  }

  // If empty array idVecB then compare with idVecA.
  if (nVecB == 0) {
    for (int i = 0; i < nVecA; ++i)
      if (idChk1 == idVecA[i] || idChk2 == idVecA[i]) return true;
    return false;
  }

  // If empty array idVecA then compare with idVecB.
  if (nVecA == 0) {
    for (int i = 0; i < nVecB; ++i)
      if (idChk1 == idVecB[i] || idChk2 == idVecB[i]) return true;
    return false;
  }

  // Else check that pair matches allowed combinations.
  for (int i = 0; i < nVecA; ++i)
  for (int j = 0; j < nVecB; ++j)
    if ( (idChk1 == idVecA[i] && idChk2 == idVecB[j])
      || (idChk2 == idVecA[i] && idChk1 == idVecB[j]) ) return true;
  return false;

}

//==========================================================================

} // end namespace Pythia8
