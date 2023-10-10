
// re-writtten version of SigmaProcess for GPU parallelism
// original class structure
// > SigmaProcess
//   > Sigma0Process
//   > Sigma1Process
//   > Sigma2Process
//   > Sigma3Process
//   > SigmaLHAProcess



#ifndef Pythia8_SigmaProcessGPU_H
#define Pythia8_SigmaProcessGPU_H

#include "Pythia8/BasicsGPU.h"
// #include "Pythia8/BeamParticle.h"
#include "Pythia8/EventGPU.h"
#include "Pythia8/InfoGPU.h"
// #include "Pythia8/LesHouchesGPU.h"
#include "Pythia8/ParticleDataGPU.h"
// #include "Pythia8/PartonDistributions.h"
#include "Pythia8/PythiaComplexGPU.h"
#include "Pythia8/PythiaStdlibGPU.h"
#include "Pythia8/ResonanceWidthsGPU.h"
#include "Pythia8/SettingsGPU.h"
// #include "Pythia8/SigmaTotal.h"
#include "Pythia8/StandardModelGPU.h"
// #include "Pythia8/SLHAinterface.h"
#include "Pythia8/SusyLesHouchesGPU.h"

namespace Pythia8 
{

  ParticleData* particleDataPtr;

  struct InBeamDataGPU
  {
    int    id=0;
    double pdf=0;
  };

  struct InPairDataGPU
  {
    int    idA=0, idB=0;
    double pdfA=0, pdfB=0, pdfSigma=0;
  };

  struct SigmaProcessData 
  {

    // missing
    // Info* infoPtr;
    // Settings* settingsPtr;
    // Rndm* rndmPtr;
    // BeamParticle* beamAPtr;
    // BeamParticle* beamBPtr;
    // Couplings* couplingsPtr;
    // SigmaTotal* sigmaTotPtr;
    // SusyLesHouches* slhaPtr;
    // LHAup* lhaUpPtr;


    // index for GPU thread
    int index = -1;

    // Conversion of GeV^{-2} to mb for cross section.
    static constexpr double CONVERT2MB = 0.389380;
    
    // The sum of outgoing masses must not be too close to the cm energy.
    static constexpr double MASSMARGIN = 0.1;

    // Parameters of momentum rescaling procedure: maximally allowed
    // relative energy error and number of iterations.
    static constexpr double COMPRELERR = 1e-10;
    static constexpr int    NCOMPSTEP = 10;

    int_star     nQuarkIn, renormScale1, renormScale2, renormScale3, renormScale3VV, // Initialization data, normally only set once.
                 factorScale1, factorScale2, factorScale3, factorScale3VV;
    double_star  Kfactor, mcME, mbME, mmuME, mtauME, renormMultFac, renormFixScale,
                 factorMultFac, factorFixScale;
    int_star    higgsH1parity, higgsH2parity, higgsA3parity; // CP violation parameters for Higgs sector, normally only set once.
    double_star higgsH1eta, higgsH2eta, higgsA3eta, higgsH1phi, higgsH2phi,
                higgsA3phi;
    int_star    idA, idB; // Information on incoming beams.
    double_star mA, mB;
    bool_star   isLeptonA, isLeptonB, hasLeptonBeams;
    InPairDataGPU** inBeamA; // Partons in beams, with PDF's.
    InPairDataGPU** inBeamB;
    InPairDataGPU** inPair; // Allowed colliding parton pairs, with pdf's.
    double_star   mH, sH, sH2; // Store common subprocess kinematics quantities.
    double_star   Q2RenSave, alpEM, alpS, Q2FacSave, x1Save, x2Save, pdf1Save, // Store Q2 renormalization and factorization scales, and related values.
                  pdf2Save, sigmaSumSave;
    int_star      id1, id2, id3, id4, id5; // Store flavour, colour, anticolour, mass, angles and the whole particle.
    int_star      idSave[12], colSave[12], acolSave[12];
    double_star   mSave[12], cosTheta, sinTheta, phi, sHMass, sHBeta, pT2Mass, pTFin;
    Particle*     parton[12];
    Particle*     partonT[12]; // Minimal set of saved kinematics for trial interactions when using the x-dependent matter profile of multiparton interactions.
    double_star   mSaveT[12], pTFinT, cosThetaT, sinThetaT, phiT;
    double_star   mME[12];
    Vec4*         pME[12];
    bool_star     swapTU; // Store whether tHat and uHat are swapped (= same as swap 3 and 4).

    // Sigma0
    // (none)

    // Sigma1
    // (none)

    // Sigma2
    double_star   tH, uH, tH2, uH2, m3, s3, m4, s4, pT2, runBW3, runBW4; // Store subprocess kinematics quantities.

    // Sigma3
    double_star   m3, s3, m4, s4, m5, s5, runBW3, runBW4, runBW5; // Store subprocess kinematics quantities.
    Vec4*         p3cm, p4cm, p5cm;  

    // SigmaLHA
    // (none)


    // susy stuff

    CoupSUSY* coupSUSYPtr;
  
  };



  struct ProcessContainerData 
  {

    // --- inputs ---
    // (using constructor)
    // (generatet in ProcessLevel::init and passed in)
    // unique for each processContainer
    SigmaProcessData*    sigmaProcessPtr; // Pointer to the subprocess matrix element. Mark if external.
    
    // set to 0
    // also unique; there are several types of them
    PhaseSpace*      phaseSpacePtr; // Pointer to the phase space generator.
    // (using init function)
    ResonanceDecays* resDecaysPtr; // Pointer to ResonanceDecays object for sequential resonance decays.

    // --- private ---
    static const int N12SAMPLE, N3SAMPLE; // Constants: could only be changed in the code itself.
    bool             externalPtr;
    bool   matchInOut; // Possibility to modify Les Houches input.
    int    idRenameBeams, setLifetime, setLeptonMass, idLep[3];
    double mRecalculate, mLep[3];
    bool   isLHA, isNonDiff, isResolved, isDiffA, isDiffB, isDiffC, isQCD3body, // Info on process.
          allowNegSig, isSameSave, increaseMaximum, canVetoResDecay;
    int    lhaStrat, lhaStratAbs;
    bool   useStrictLHEFscales;
    bool   newSigmaMx; // Statistics on generation process. (Long integers just in case.)
    long   nTry, nSel, nAcc, nTryStat;
    double sigmaMx, sigmaSgn, sigmaSum, sigma2Sum, sigmaNeg, sigmaAvg,
          sigmaFin, deltaFin, weightNow, wtAccSum;
    vector<int> codeLHA; // Statistics for Les Houches event classification.
    vector<long> nTryLHA, nSelLHA, nAccLHA;
  };

}

#endif