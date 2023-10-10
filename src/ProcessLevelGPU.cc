#include "Pythia8/ProcessLevelGPU.h"


  namespace SigmaProcess
  {

    // --- public ---

    // Destructor.
    virtual ~SigmaProcess() {}

    // Perform simple initialization and store pointers.
    void init(Info* infoPtrIn, Settings* settingsPtrIn,
      ParticleData* particleDataPtrIn, Rndm* rndmPtrIn,
      BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn, Couplings* couplings,
      SigmaTotal* sigmaTotPtrIn = 0, SLHAinterface* slhaInterfacePtrIn = 0);

    // Store or replace Les Houches pointer.
    void setLHAPtr( LHAup* lhaUpPtrIn) {lhaUpPtr = lhaUpPtrIn;}

    // Initialize process. Only used for some processes.
    virtual void initProc() {}

    // Set up allowed flux of incoming partons. Default is no flux.
    virtual bool initFlux();

    // Input and complement kinematics for resolved 2 -> 1 process.
    // Usage: set1Kin( x1in, x2in, sHin).
    virtual void set1Kin( double , double , double ) {}

    // Input and complement kinematics for resolved 2 -> 2 process.
    // Usage: set2Kin( x1in, x2in, sHin, tHin, m3in, m4in, runBW3in, runBW4in).
    virtual void set2Kin( double , double , double , double , double ,
      double, double, double ) {}

    // Ditto, but for Multiparton Interactions applications, so different input.
    // Usage: set2KinMPI( x1in, x2in, sHin, tHin, uHin,
    //                   alpSin, alpEMin, needMasses, m3in, m4in)
    virtual void set2KinMPI( double , double , double , double ,
      double , double , double , bool , double , double ) {}

    // Input and complement kinematics for resolved 2 -> 3 process.
    // Usage: set3Kin( x1in, x2in, sHin, p3prel, p4prel, p5prel,
    //                 m3in, m4in, m5in, runBW3in, runBW4in, runBW5in);
    virtual void set3Kin( double , double , double , Vec4 , Vec4 , Vec4 ,
      double , double , double , double , double , double ) {}

    // Calculate flavour-independent parts of cross section.
    virtual void sigmaKin() {}

    // Evaluate sigma for unresolved, sigmaHat(sHat) for 2 -> 1 processes,
    // d(sigmaHat)/d(tHat) for (resolved) 2 -> 2 processes, and |M|^2 for
    // 2 -> 3 processes. Answer in "native" units, either mb or GeV^-2.
    virtual double sigmaHat() {return 0.;}

    // Wrapper to sigmaHat, to (a) store current incoming flavours and
    // (b) convert from GeV^-2 to mb where required.
    // For 2 -> 1/2 also (c) convert from from |M|^2 to d(sigmaHat)/d(tHat).
    virtual double sigmaHatWrap(int id1in = 0, int id2in = 0) {
      id1 = id1in; id2 = id2in;
      return ( convert2mb() ? CONVERT2MB * sigmaHat() : sigmaHat() ); }

    // Convolute above with parton flux and K factor. Sum over open channels.
    virtual double sigmaPDF();

    // Select incoming parton channel and extract parton densities (resolved).
    void pickInState(int id1in = 0, int id2in = 0);

    // Select flavour, colour and anticolour.
    virtual void setIdColAcol() {}

    // Perform kinematics for a Multiparton Interaction, in its rest frame.
    virtual bool final2KinMPI( int = 0, int = 0, Vec4 = 0., Vec4 = 0.,
      double = 0., double = 0.) {return true;}

    // Evaluate weight for simultaneous flavours (only gamma*/Z0 gamma*/Z0).
    // Usage: weightDecayFlav( process).
    virtual double weightDecayFlav( Event&) {return 1.;}

    // Evaluate weight for decay angular configuration.
    // Usage: weightDecay( process, iResBeg, iResEnd), where
    // iResBeg <= i < iResEnd is range of sister partons to test decays of.
    virtual double weightDecay( Event&, int, int) {return 1.;}

    // Set scale, when that is missing for an external LHA process.
    virtual void setScale() {}

    // Process name and code, and the number of final-state particles.
    virtual string name()            const {return "unnamed process";}
    virtual int    code()            const {return 0;}
    virtual int    nFinal()          const {return 2;}

    // Need to know which incoming partons to set up interaction for.
    virtual string inFlux()          const {return "unknown";}

    // Need to know whether to convert cross section answer from GeV^-2 to mb.
    virtual bool   convert2mb()      const {return true;}

    // For 2 -> 2 process optional conversion from |M|^2 to d(sigmaHat)/d(tHat).
    virtual bool   convertM2()       const {return false;}

    // Special treatment needed for Les Houches processes.
    virtual bool   isLHA()           const {return false;}

    // Special treatment needed for elastic and diffractive processes.
    virtual bool   isNonDiff()       const {return false;}
    virtual bool   isResolved()      const {return true;}
    virtual bool   isDiffA()         const {return false;}
    virtual bool   isDiffB()         const {return false;}
    virtual bool   isDiffC()         const {return false;}

    // Special treatment needed for SUSY processes.
    virtual bool   isSUSY()          const {return false;}

    // Special treatment needed if negative cross sections allowed.
    virtual bool   allowNegativeSigma() const {return false;}

    // Flavours in 2 -> 2/3 processes where masses needed from beginning.
    // (For a light quark masses will be used in the final kinematics,
    // but not at the matrix-element level. For a gluon no masses at all.)
    virtual int    id3Mass()         const {return 0;}
    virtual int    id4Mass()         const {return 0;}
    virtual int    id5Mass()         const {return 0;}

    // Special treatment needed if process contains an s-channel resonance.
    virtual int    resonanceA()      const {return 0;}
    virtual int    resonanceB()      const {return 0;}

    // 2 -> 2 and 2 -> 3 processes only through s-channel exchange.
    virtual bool   isSChannel()      const {return false;}

    // NOAM: Insert an intermediate resonance in 2 -> 1 -> 2 (or 3) listings.
    virtual int    idSChannel()      const {return 0;}

    // QCD 2 -> 3 processes need special phase space selection machinery.
    virtual bool   isQCD3body()      const {return false;}

    // Special treatment in 2 -> 3 with two massive propagators.
    virtual int    idTchan1()        const {return 0;}
    virtual int    idTchan2()        const {return 0;}
    virtual double tChanFracPow1()   const {return 0.3;}
    virtual double tChanFracPow2()   const {return 0.3;}
    virtual bool   useMirrorWeight() const {return false;}

    // Special process-specific gamma*/Z0 choice if >=0 (e.g. f fbar -> H0 Z0).
    virtual int    gmZmode()         const {return -1;}

    // Tell whether tHat and uHat are swapped (= same as swap 3 and 4).
    bool swappedTU()          const {return swapTU;}

    // Give back particle properties: flavours, colours, masses, or all.
    int    id(int i)          const {return idSave[i];}
    int    col(int i)         const {return colSave[i];}
    int    acol(int i)        const {return acolSave[i];}
    double m(int i)           const {return mSave[i];}
    Particle getParton(int i) const {return parton[i];}

    // Give back couplings and parton densities.
    // Not all known for nondiffractive.
    double Q2Ren()            const {return Q2RenSave;}
    double alphaEMRen()       const {return alpEM;}
    double alphaSRen()        const {return alpS;}
    double Q2Fac()            const {return Q2FacSave;}
    double pdf1()             const {return pdf1Save;}
    double pdf2()             const {return pdf2Save;}

    // Give back angles; relevant only for multipe-interactions processes.
    double thetaMPI()         const {return atan2( sinTheta, cosTheta);}
    double phiMPI()           const {return phi;}
    double sHBetaMPI()        const {return sHBeta;}
    double pT2MPI()           const {return pT2Mass;}
    double pTMPIFin()         const {return pTFin;}

    // Save and load kinematics for trial interactions
    void saveKin() {
      for (int i = 0; i < 12; i++) { partonT[i] = parton[i];
        mSaveT[i] = mSave[i]; }
      pTFinT = pTFin; phiT = phi; cosThetaT = cosTheta; sinThetaT = sinTheta; }
    void loadKin() {
      for (int i = 0; i < 12; i++) { parton[i] = partonT[i];
      mSave[i] = mSaveT[i]; }
      pTFin = pTFinT; cosTheta = cosThetaT; sinTheta = sinThetaT; phi = phiT;
    }
    void swapKin() {
      for (int i = 0; i < 12; i++) { swap(parton[i], partonT[i]);
                                    swap(mSave[i], mSaveT[i]); }
      swap(pTFin, pTFinT); swap(cosTheta, cosThetaT);
      swap(sinTheta, sinThetaT); swap(phi, phiT); }

    // --- private ---

    // Constructor.
    SigmaProcess() : infoPtr(0), settingsPtr(0), particleDataPtr(0),
      rndmPtr(0), beamAPtr(0), beamBPtr(0), couplingsPtr(0), sigmaTotPtr(0),
      slhaPtr(0), lhaUpPtr(0) {for (int i = 0; i < 12; ++i) mSave[i] = 0.;
      Q2RenSave = alpEM = alpS = Q2FacSave = pdf1Save = pdf2Save = 0.; }

    // Set flavour, colour and anticolour.
    void setId( int id1in = 0, int id2in = 0, int id3in = 0, int id4in = 0,
      int id5in = 0) {idSave[1] = id1in; idSave[2] = id2in; idSave[3] = id3in;
      idSave[4] = id4in; idSave[5] = id5in;}
    void setColAcol( int col1 = 0, int acol1 = 0,
      int col2 = 0, int acol2 = 0, int col3 = 0, int acol3 = 0,
      int col4 = 0, int acol4 = 0, int col5 = 0, int acol5 = 0) {
      colSave[1] = col1; acolSave[1] = acol1; colSave[2] = col2;
      acolSave[2] = acol2; colSave[3] = col3; acolSave[3] = acol3;
      colSave[4] = col4; acolSave[4] = acol4; colSave[5] = col5;
      acolSave[5] = acol5; }
    void swapColAcol() { swap(colSave[1], acolSave[1]);
      swap(colSave[2], acolSave[2]); swap(colSave[3], acolSave[3]);
      swap(colSave[4], acolSave[4]); swap(colSave[5], acolSave[5]);}
    void swapCol1234() { swap(colSave[1], colSave[2]);
      swap(colSave[3], colSave[4]); swap(acolSave[1], acolSave[2]);
      swap(acolSave[3], acolSave[4]);}
    void swapCol12() { swap(colSave[1], colSave[2]);
      swap(acolSave[1], acolSave[2]);}
    void swapCol34() { swap(colSave[3], colSave[4]);
      swap(acolSave[3], acolSave[4]);}
    double weightTopDecay( Event& process, int iResBeg, int iResEnd); // Common code for top and Higgs secondary decay angular weights.
    double weightHiggsDecay( Event& process, int iResBeg, int iResEnd);

    void addBeamA(int idIn) {inBeamA.push_back(InBeam(idIn));}
    void addBeamB(int idIn) {inBeamB.push_back(InBeam(idIn));}
    int sizeBeamA() const {return inBeamA.size();}
    int sizeBeamB() const {return inBeamB.size();}

    void addPair(int idAIn, int idBIn) {
      inPair.push_back(InPair(idAIn, idBIn));}
    int sizePair() const {return inPair.size();}


    virtual bool setupForME() {return true;} // Calculate and store all modified masses and four-vectors intended for matrix elements. Return false if failed.
    bool     setupForMEin();

    // -- Sigma0Process

    // -- public

    // Destructor.
    virtual ~Sigma0Process() {}

    // Number of final-state particles.
    virtual int    nFinal() const {return 2;};

    // No partonic flux to be set up.
    virtual bool   initFlux() {return true;}

    // Evaluate sigma for unresolved processes.
    virtual double sigmaHat() {return 0.;}

    // Since no PDF's there is no difference from above.
    virtual double sigmaPDF() {return sigmaHat();}

    // Answer for these processes already in mb, so do not convert.
    virtual bool convert2mb() const {return false;}

    // -- private

    // Constructor.
    Sigma0Process() {}



    // -- SigmaLHAProcess

    // Constructor.
    SigmaLHAProcess() {}

    // Destructor.
    virtual ~SigmaLHAProcess() {}

    // No partonic flux to be set up.
    virtual bool   initFlux() {return true;}

    // Dummy function: action is put in PhaseSpaceLHA.
    virtual double sigmaPDF() {return 1.;}

    // Evaluate weight for decay angular configuration, where relevant.
    virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

    // Set scale, when that is missing for an external LHA process.
    virtual void   setScale();

    // Info on the subprocess.
    virtual string name()     const {return "Les Houches User Process(es)";}
    virtual int    code()     const {return 9999;}

    // Number of final-state particles depends on current process choice.
    virtual int    nFinal()   const;

    // Answer for these processes not in GeV^-2, so do not do this conversion.
    virtual bool   convert2mb() const {return false;}

    // Ensure special treatment of Les Houches processes.
    virtual bool   isLHA()    const {return true;}

    // Special treatment needed if negative cross sections allowed.
    virtual bool   allowNegativeSigma() const {
      return (lhaUpPtr->strategy() < 0);}


    // -- Sigma1Process --


    // -- public

    // Destructor.
    virtual ~Sigma1Process() {}

    // Number of final-state particles.
    virtual int    nFinal() const {return 1;};

    // Input and complement kinematics for resolved 2 -> 1 process.
    virtual void   set1Kin( double x1in, double x2in, double sHin) {
      store1Kin( x1in, x2in, sHin); sigmaKin();}

    // Evaluate sigmaHat(sHat) for resolved 2 -> 1 processes.
    virtual double sigmaHat() {return 0.;}

    // Wrapper to sigmaHat, to (a) store current incoming flavours,
    // (b) convert from GeV^-2 to mb where required, and
    // (c) convert from |M|^2 to d(sigmaHat)/d(tHat) where required.
    virtual double sigmaHatWrap(int id1in = 0, int id2in = 0);

    // -- private

    // Constructor.
    Sigma1Process() {}

    // Store kinematics and set scales for resolved 2 -> 1 process.
    virtual void   store1Kin( double x1in, double x2in, double sHin);

    // Calculate modified masses and four-vectors for matrix elements.
    virtual bool   setupForME();

    // -- Sigma2Process


     // --- public ---

    // Destructor.
    virtual ~Sigma2Process() {}

    // Number of final-state particles.
    virtual int    nFinal() const {return 2;};

    // Input and complement kinematics for resolved 2 -> 2 process.
    virtual void   set2Kin( double x1in, double x2in, double sHin,
      double tHin, double m3in, double m4in, double runBW3in,
      double runBW4in) { store2Kin( x1in, x2in, sHin, tHin, m3in, m4in,
      runBW3in, runBW4in); sigmaKin();}

    // Ditto, but for Multiparton Interactions applications, so different input.
    virtual void   set2KinMPI( double x1in, double x2in, double sHin,
      double tHin, double uHin, double alpSin, double alpEMin,
      bool needMasses, double m3in, double m4in) {
      store2KinMPI( x1in, x2in, sHin, tHin, uHin, alpSin, alpEMin,
      needMasses, m3in, m4in); sigmaKin();}

    // Evaluate d(sigmaHat)/d(tHat) for resolved 2 -> 2 processes.
    virtual double sigmaHat() {return 0.;}

    // Wrapper to sigmaHat, to (a) store current incoming flavours,
    // (b) convert from GeV^-2 to mb where required, and
    // (c) convert from |M|^2 to d(sigmaHat)/d(tHat) where required.
    virtual double sigmaHatWrap(int id1in = 0, int id2in = 0) {
      id1 = id1in; id2 = id2in; double sigmaTmp = sigmaHat();
      if (convertM2())  sigmaTmp /= 16. * M_PI * sH2;
      if (convert2mb()) sigmaTmp *= CONVERT2MB; return sigmaTmp;}

    // Perform kinematics for a Multiparton Interaction, in its rest frame.
    virtual bool   final2KinMPI( int i1Res = 0, int i2Res = 0, Vec4 p1Res = 0.,
      Vec4 p2Res = 0., double m1Res = 0., double m2Res = 0.);

    // --- private --- 

    // Constructor.
    Sigma2Process() : tH(0.), uH(0.), tH2(0.), uH2(0.), m3(0.), s3(0.),
      m4(0.), s4(0.), pT2(0.), runBW3(0.), runBW4(0.) {}

    // Store kinematics and set scales for resolved 2 -> 2 process.
    virtual void   store2Kin( double x1in, double x2in, double sHin,
      double tHin, double m3in, double m4in, double runBW3in,
      double runBW4in);
    virtual void   store2KinMPI( double x1in, double x2in, double sHin,
      double tHin, double uHin, double alpSin, double alpEMin,
      bool needMasses, double m3in, double m4in);

    // Calculate modified masses and four-vectors for matrix elements.
    virtual bool   setupForME();



    // -- Sigma3Process

    // --- public ---

    // Destructor.
    virtual ~Sigma3Process() {}

    // Number of final-state particles.
    virtual int    nFinal() const {return 3;};

    // Input and complement kinematics for resolved 2 -> 3 process.
    virtual void   set3Kin( double x1in, double x2in, double sHin,
      Vec4 p3cmIn, Vec4 p4cmIn, Vec4 p5cmIn, double m3in, double m4in,
      double m5in, double runBW3in, double runBW4in, double runBW5in) {
      store3Kin( x1in, x2in, sHin, p3cmIn, p4cmIn, p5cmIn, m3in, m4in, m5in,
      runBW3in, runBW4in, runBW5in); sigmaKin();}

    // Evaluate d(sigmaHat)/d(tHat) for resolved 2 -> 3 processes.
    virtual double sigmaHat() {return 0.;}

    // --- private ---

    // Constructor.
    Sigma3Process() {}

    // Store kinematics and set scales for resolved 2 -> 3 process.
    virtual void   store3Kin( double x1in, double x2in, double sHin,
      Vec4 p3cmIn, Vec4 p4cmIn, Vec4 p5cmIn, double m3in, double m4in,
      double m5in, double runBW3in, double runBW4in, double runBW5in);

    // Calculate modified masses and four-vectors for matrix elements.
    virtual bool   setupForME();

  }

  namespace SetupContainers
  {
    // --- public ---

    // Constructor.
    SetupContainers() {}

    // Initialization assuming all necessary data already read.
    bool init(vector<ProcessContainer*>& containerPtrs, Info* infoPtr,
      Settings& settings, ParticleData* particleDataPtr, Couplings* couplings);

    // Initialization of a second hard process.
    bool init2(vector<ProcessContainer*>& container2Ptrs, Settings& settings);

    // --- private ---

    // Methods to check that outgoing SUSY particles are allowed ones.
    void setupIdVecs( Settings& settings);
    bool allowIdVals( int idCheck1, int idCheck2);
  }

  namespace ProcessContainer
  {
    // --- public ---

    // Constructor.
    ProcessContainer(SigmaProcess* sigmaProcessPtrIn = 0,
      bool externalPtrIn = false, PhaseSpace* phaseSpacePtrIn = 0) :
        sigmaProcessPtr(sigmaProcessPtrIn),
        externalPtr(externalPtrIn), phaseSpacePtr(phaseSpacePtrIn) {}

    // Destructor. Do not destroy external sigmaProcessPtr.
    ~ProcessContainer() {delete phaseSpacePtr;
      if (!externalPtr) delete sigmaProcessPtr;}

    // Initialize phase space and counters.
    bool init(bool isFirst, Info* infoPtrIn, Settings& settings,
      ParticleData* particleDataPtrIn, Rndm* rndmPtrIn, BeamParticle* beamAPtr,
      BeamParticle* beamBPtr, Couplings* couplings, SigmaTotal* sigmaTotPtr,
      ResonanceDecays* resDecaysPtrIn, SLHAinterface* slhaInterfacePtr,
      UserHooks* userHooksPtr);

    // Store or replace Les Houches pointer.
    void setLHAPtr( LHAup* lhaUpPtrIn,  ParticleData* particleDataPtrIn = 0)
      {lhaUpPtr = lhaUpPtrIn;
      if (particleDataPtrIn != 0) particleDataPtr = particleDataPtrIn;
      if (sigmaProcessPtr != 0) sigmaProcessPtr->setLHAPtr(lhaUpPtr);
      if (phaseSpacePtr != 0) phaseSpacePtr->setLHAPtr(lhaUpPtr);}

    // Update the CM energy of the event.
    void newECM(double eCM) {phaseSpacePtr->newECM(eCM);}

    // Generate a trial event; accepted or not.
    bool trialProcess();

    // Pick flavours and colour flow of process.
    bool constructState();

    // Give the hard subprocess (with option for a second hard subprocess).
    bool constructProcess( Event& process, bool isHardest = true);

    // Give the Les Houches decay chain for external resonance.
    bool constructDecays( Event& process);

    // Do resonance decays.
    bool decayResonances( Event& process);

    // Accumulate statistics after user veto.
    void accumulate();

    // Reset statistics on events generated so far.
    void reset();

    // Process name and code, and the number of final-state particles.
    string name()        const {return sigmaProcessPtr->name();}
    int    code()        const {return sigmaProcessPtr->code();}
    int    nFinal()      const {return sigmaProcessPtr->nFinal();}
    bool   isSUSY()      const {return sigmaProcessPtr->isSUSY();}

    // Member functions for info on generation process.
    bool   newSigmaMax() const {return newSigmaMx;}
    double sigmaMax()    const {return sigmaMx;}
    long   nTried()      const {return nTry;}
    long   nSelected()   const {return nSel;}
    long   nAccepted()   const {return nAcc;}
    double weightSum()   const {return wtAccSum;}
    double sigmaSelMC()  {if (nTry > nTryStat) sigmaDelta(); return sigmaAvg;}
    double sigmaMC()     {if (nTry > nTryStat) sigmaDelta(); return sigmaFin;}
    double deltaMC()     {if (nTry > nTryStat) sigmaDelta(); return deltaFin;}

    // Some kinematics quantities.
    int    id1()         const {return sigmaProcessPtr->id(1);}
    int    id2()         const {return sigmaProcessPtr->id(2);}
    double x1()          const {return phaseSpacePtr->x1();}
    double x2()          const {return phaseSpacePtr->x2();}
    double Q2Fac()       const {return sigmaProcessPtr->Q2Fac();}
    double mHat()        const {return sqrtpos(phaseSpacePtr->sHat());}
    double pTHat()       const {return phaseSpacePtr->pTHat();}

    // Tell whether container is for Les Houches events.
    bool   isLHAContainer() const {return isLHA;}
    int    lhaStrategy()    const {return lhaStrat;}

    // Info on Les Houches events.
    int    codeLHASize()       const {return codeLHA.size();}
    int    subCodeLHA(int i)   const {return codeLHA[i];}
    long   nTriedLHA(int i)    const {return nTryLHA[i];}
    long   nSelectedLHA(int i) const {return nSelLHA[i];}
    long   nAcceptedLHA(int i) const {return nAccLHA[i];}

    // When two hard processes set or get info whether process is matched.
    void   isSame( bool isSameIn) { isSameSave = isSameIn;}
    bool   isSame()      const {return isSameSave;}

    // -- private ---

    // Estimate integrated cross section and its uncertainty.
    void sigmaDelta();

  }

  namespace ProcessLevel
  {

    // --- public ---

    // Constructor.
    void ProcessLevel_constructor(ProcessLevelData& d) 
    {
      d.iLHACont = -1;
    }

    // Destructor to delete processes in containers.
    void ProcessLevel_destructor(ProcessLevelData& d)
    {
      // Run through list of first hard processes and delete them.
      for (int i = 0; i < int(d.containerPtrs.size()); ++i)
        delete d.containerPtrs[i];

      // Run through list of second hard processes and delete them.
      for (int i = 0; i < int(d.container2Ptrs.size()); ++i)
        delete d.container2Ptrs[i];
    }

// Main routine to initialize generation process.

bool ProcessLevel::init( Info* infoPtrIn, Settings& settings,
  ParticleData* particleDataPtrIn, Rndm* rndmPtrIn,
  BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
  Couplings* couplingsPtrIn, SigmaTotal* sigmaTotPtrIn, bool doLHA,
  SLHAinterface* slhaInterfacePtrIn, UserHooks* userHooksPtrIn,
  vector<SigmaProcess*>& sigmaPtrs, vector<PhaseSpace*>& phaseSpacePtrs,
  ostream& os) {

  // Store input pointers for future use.
  infoPtr          = infoPtrIn;
  particleDataPtr  = particleDataPtrIn;
  rndmPtr          = rndmPtrIn;
  beamAPtr         = beamAPtrIn;
  beamBPtr         = beamBPtrIn;
  couplingsPtr     = couplingsPtrIn;
  sigmaTotPtr      = sigmaTotPtrIn;
  userHooksPtr     = userHooksPtrIn;
  slhaInterfacePtr = slhaInterfacePtrIn;

  // Send on some input pointers.
  resonanceDecays.init( infoPtr, particleDataPtr, rndmPtr);

  // Set up SigmaTotal. Store sigma_nondiffractive for future use.
  sigmaTotPtr->init( infoPtr, settings, particleDataPtr);
  int    idA = infoPtr->idA();
  int    idB = infoPtr->idB();
  double eCM = infoPtr->eCM();
  sigmaTotPtr->calc( idA, idB, eCM);
  sigmaND = sigmaTotPtr->sigmaND();

  // Options to allow second hard interaction and resonance decays.
  doSecondHard  = settings.flag("SecondHard:generate");
  doSameCuts    = settings.flag("PhaseSpace:sameForSecond");
  doResDecays   = settings.flag("ProcessLevel:resonanceDecays");
  startColTag   = settings.mode("Event:startColTag");

  // Second interaction not to be combined with biased phase space.
  if (doSecondHard && userHooksPtr != 0
  && userHooksPtr->canBiasSelection()) {
    infoPtr->errorMsg("Error in ProcessLevel::init: "
      "cannot combine second interaction with biased phase space");
    return false;
  }

  // Mass and pT cuts for two hard processes.
  mHatMin1      = settings.parm("PhaseSpace:mHatMin");
  mHatMax1      = settings.parm("PhaseSpace:mHatMax");
  if (mHatMax1 < mHatMin1) mHatMax1 = eCM;
  pTHatMin1     = settings.parm("PhaseSpace:pTHatMin");
  pTHatMax1     = settings.parm("PhaseSpace:pTHatMax");
  if (pTHatMax1 < pTHatMin1) pTHatMax1 = eCM;
  mHatMin2      = settings.parm("PhaseSpace:mHatMinSecond");
  mHatMax2      = settings.parm("PhaseSpace:mHatMaxSecond");
  if (mHatMax2 < mHatMin2) mHatMax2 = eCM;
  pTHatMin2     = settings.parm("PhaseSpace:pTHatMinSecond");
  pTHatMax2     = settings.parm("PhaseSpace:pTHatMaxSecond");
  if (pTHatMax2 < pTHatMin2) pTHatMax2 = eCM;

  // Check whether mass and pT ranges agree or overlap.
  cutsAgree     = doSameCuts;
  if (mHatMin2 == mHatMin1 && mHatMax2 == mHatMax1 && pTHatMin2 == pTHatMin1
      && pTHatMax2 == pTHatMax1) cutsAgree = true;
  cutsOverlap   = cutsAgree;
  if (!cutsAgree) {
    bool mHatOverlap = (max( mHatMin1, mHatMin2)
                      < min( mHatMax1, mHatMax2));
    bool pTHatOverlap = (max( pTHatMin1, pTHatMin2)
                       < min( pTHatMax1, pTHatMax2));
    if (mHatOverlap && pTHatOverlap) cutsOverlap = true;
  }

  // Set up containers for all the internal hard processes.
  SetupContainers setupContainers;
  setupContainers.init(containerPtrs, infoPtr, settings, particleDataPtr,
                       couplingsPtr);

  // Append containers for external hard processes, if any.
  if (sigmaPtrs.size() > 0) {
    for (int iSig = 0; iSig < int(sigmaPtrs.size()); ++iSig)
      containerPtrs.push_back( new ProcessContainer(sigmaPtrs[iSig],
        true, phaseSpacePtrs[iSig]) );
  }

  // Append single container for Les Houches processes, if any.
  if (doLHA) {
    SigmaProcess* sigmaPtr = new SigmaLHAProcess();
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );

    // Store location of this container, and send in LHA pointer.
    iLHACont = containerPtrs.size() - 1;
    containerPtrs[iLHACont]->setLHAPtr(lhaUpPtr);
  }

  // If no processes found then refuse to do anything.
  if ( int(containerPtrs.size()) == 0) {
    infoPtr->errorMsg("Error in ProcessLevel::init: "
      "no process switched on");
    return false;
  }

  // Check whether pT-based weighting in 2 -> 2 is requested.
  if (settings.flag("PhaseSpace:bias2Selection")) {
    bool bias2Sel = false;
    if (sigmaPtrs.size() == 0 && !doLHA && !doSecondHard) {
      bias2Sel = true;
      for (int i = 0; i < int(containerPtrs.size()); ++i) {
        if (containerPtrs[i]->nFinal() != 2) bias2Sel = false;
        int code = containerPtrs[i]->code();
        if (code > 100 && code < 110) bias2Sel = false;
      }
    }
    if (!bias2Sel) {
      infoPtr->errorMsg("Error in ProcessLevel::init: "
        "requested event weighting not possible");
      return false;
    }
  }

  // Check that SUSY couplings were indeed initialized where necessary.
  bool hasSUSY = false;
  for (int i = 0; i < int(containerPtrs.size()); ++i)
    if (containerPtrs[i]->isSUSY()) hasSUSY = true;

  // If SUSY processes requested but no SUSY couplings present
  if(hasSUSY && !couplingsPtr->isSUSY) {
    infoPtr->errorMsg("Error in ProcessLevel::init: "
      "SUSY process switched on but no SUSY couplings found");
    return false;
  }

  // Fill SLHA blocks SMINPUTS and MASS from PYTHIA SM parameter values.
  slhaInterfacePtr->pythia2slha(particleDataPtr);

  // Initialize each process.
  int numberOn = 0;
  for (int i = 0; i < int(containerPtrs.size()); ++i)
    if (containerPtrs[i]->init(true, infoPtr, settings, particleDataPtr,
      rndmPtr, beamAPtr, beamBPtr, couplingsPtr, sigmaTotPtr,
      &resonanceDecays, slhaInterfacePtr, userHooksPtr)) ++numberOn;

  // Sum maxima for Monte Carlo choice.
  // NOTE: Gambit hack: Catch nans and infinities while summing up sigmas.
  sigmaMaxSum = 0.;
  bool valid = true;
  for (int i = 0; i < int(containerPtrs.size()); ++i) {
    if(std::isfinite(containerPtrs[i]->sigmaMax()))
      sigmaMaxSum += containerPtrs[i]->sigmaMax();
    else {
      std::cerr<<"\n\n\n ERROR: in Pythia8::ProcessLevel::init:\n";
      std::cerr<<"   Non-finite xsec: "<<containerPtrs[i]->sigmaMax()<<"\n";
      std::cerr<<"   Process code: "<<containerPtrs[i]->code();
      std::cerr<<",  Process: "<<containerPtrs[i]->name()<<"\n";
      std::cerr<<"This model is invalid.\n\n\n";
      delete containerPtrs[i];
      containerPtrs.erase(containerPtrs.begin() + i);
      i--;
      valid = false;
    }
  }
  if (!valid) {
    infoPtr->errorMsg("Error in ProcessLevel::init: Non-finite xsecs");
    return false;
  }

  // Option to pick a second hard interaction: repeat as above.
  int number2On = 0;
  if (doSecondHard) {
    setupContainers.init2(container2Ptrs, settings);
    if ( int(container2Ptrs.size()) == 0) {
      infoPtr->errorMsg("Error in ProcessLevel::init: "
        "no second hard process switched on");
      return false;
    }
    for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2)
      if (container2Ptrs[i2]->init(false, infoPtr, settings, particleDataPtr,
        rndmPtr, beamAPtr, beamBPtr, couplingsPtr, sigmaTotPtr,
        &resonanceDecays, slhaInterfacePtr, userHooksPtr)) ++number2On;
    sigma2MaxSum = 0.;
    for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2)
      sigma2MaxSum += container2Ptrs[i2]->sigmaMax();
  }

  // Printout during initialization is optional.
  if (settings.flag("Init:showProcesses")) {

    // Construct string with incoming beams and for cm energy.
    string collision = "We collide " + particleDataPtr->name(idA)
      + " with " + particleDataPtr->name(idB) + " at a CM energy of ";
    string pad( 51 - collision.length(), ' ');

    // Print initialization information: header.
    os << "\n *-------  PYTHIA Process Initialization  ---------"
       << "-----------------*\n"
       << " |                                                   "
       << "               |\n"
       << " | " << collision << scientific << setprecision(3)
       << setw(9) << eCM << " GeV" << pad << " |\n"
       << " |                                                   "
       << "               |\n"
       << " |---------------------------------------------------"
       << "---------------|\n"
       << " |                                                   "
       << " |             |\n"
       << " | Subprocess                                    Code"
       << " |   Estimated |\n"
       << " |                                                   "
       << " |    max (mb) |\n"
       << " |                                                   "
       << " |             |\n"
       << " |---------------------------------------------------"
       << "---------------|\n"
       << " |                                                   "
       << " |             |\n";

    // Loop over existing processes: print individual process info.
    map<int, double> sigmaMaxM;
    map<int, string> nameM;
    for (int i = 0; i < int(containerPtrs.size()); ++i) {
      int code = containerPtrs[i]->code();
      nameM[code] = containerPtrs[i]->name();
      sigmaMaxM[code] = containerPtrs[i]->sigmaMax() > sigmaMaxM[code] ?
        containerPtrs[i]->sigmaMax() : sigmaMaxM[code];
    }
    for (map<int, string>::iterator i = nameM.begin(); i != nameM.end(); ++i)
      os << " | " << left << setw(45) << i->second
         << right << setw(5) << i->first << " | "
         << scientific << setprecision(3) << setw(11)
         << sigmaMaxM[i->first] << " |\n";

    // Loop over second hard processes, if any, and repeat as above.
    if (doSecondHard) {
      os << " |                                                   "
         << " |             |\n"
         << " |---------------------------------------------------"
         <<"---------------|\n"
         << " |                                                   "
         <<" |             |\n";
      sigmaMaxM.clear();
      nameM.clear();
      for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2) {
        int code = container2Ptrs[i2]->code();
        nameM[code] = container2Ptrs[i2]->name();
        sigmaMaxM[code] = container2Ptrs[i2]->sigmaMax() > sigmaMaxM[code] ?
          container2Ptrs[i2]->sigmaMax() : sigmaMaxM[code];
      }
      for (map<int, string>::iterator i2 = nameM.begin(); i2 != nameM.end();
           ++i2)
        os << " | " << left << setw(45) << i2->second
           << right << setw(5) << i2->first << " | "
           << scientific << setprecision(3) << setw(11)
           << sigmaMaxM[i2->first] << " |\n";
    }

    // Listing finished.
    os << " |                                                     "
       << "             |\n"
       << " *-------  End PYTHIA Process Initialization ----------"
       <<"-------------*" << endl;
  }

  /* NOTE: Gambit hack: ColliderBit has its own xsec veto... remove this:
  // If sum of maxima vanishes then refuse to do anything.
  if ( numberOn == 0  || sigmaMaxSum <= 0.) {
    infoPtr->errorMsg("Error in ProcessLevel::init: "
      "all processes have vanishing cross sections");
    return false;
  }
  *** Gambit hack end */
  if ( doSecondHard && (number2On == 0  || sigma2MaxSum <= 0.) ) {
    infoPtr->errorMsg("Error in ProcessLevel::init: "
      "all second hard processes have vanishing cross sections");
    return false;
  }

  // If two hard processes then check whether some (but not all) agree.
  allHardSame  = true;
  noneHardSame = true;
  if (doSecondHard) {
    bool foundMatch = false;

    // Check for each first process if matched in second.
    for (int i = 0; i < int(containerPtrs.size()); ++i) {
      foundMatch = false;
      if (cutsOverlap)
      for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2)
        if (container2Ptrs[i2]->code() == containerPtrs[i]->code())
          foundMatch = true;
      containerPtrs[i]->isSame( foundMatch );
      if (!foundMatch)  allHardSame = false;
      if ( foundMatch) noneHardSame = false;
    }

    // Check for each second process if matched in first.
    for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2) {
      foundMatch = false;
      if (cutsOverlap)
      for (int i = 0; i < int(containerPtrs.size()); ++i)
        if (containerPtrs[i]->code() == container2Ptrs[i2]->code())
          foundMatch = true;
      container2Ptrs[i2]->isSame( foundMatch );
      if (!foundMatch)  allHardSame = false;
      if ( foundMatch) noneHardSame = false;
    }
  }

  // Concluding classification, including cuts.
  if (!cutsAgree) allHardSame = false;
  someHardSame = !allHardSame && !noneHardSame;

  // Reset counters for average impact-parameter enhancement.
  nImpact       = 0;
  sumImpactFac  = 0.;
  sum2ImpactFac = 0.;

  // Done.
  return true;
}

    // Store or replace Les Houches pointer.
    void setLHAPtr(ProcessLevelData& d, LHAup* lhaUpPtrIn) {lhaUpPtr = lhaUpPtrIn;
      if (iLHACont >= 0) containerPtrs[iLHACont]->setLHAPtr(lhaUpPtr);}

    // Generate the next "hard" process.
    bool next(ProcessLevelData& d, Event& process);

    // Special case: LHA input of resonance decay only.
    bool nextLHAdec(ProcessLevelData& d, Event& process);

    // Accumulate and update statistics (after possible user veto).
    void accumulate();

    // Print statistics on cross sections and number of events.
    void statistics(ProcessLevelData& d, bool reset = false, ostream& os = cout);

    // Reset statistics.
    void resetStatistics(ProcessLevelData& d);

    // Add any junctions to the process event record list.
    void findJunctions(ProcessLevelData& d, Event& junEvent);

    // Initialize and call resonance decays separately.
    void initDecays(ProcessLevelData& d, Info* infoPtrIn, ParticleData* particleDataPtrIn,
      Rndm* rndmPtrIn, LHAup* lhaUpPtrIn) { infoPtr = infoPtrIn;
      resonanceDecays.init( infoPtrIn, particleDataPtrIn, rndmPtrIn);
      containerLHAdec.setLHAPtr(lhaUpPtrIn, particleDataPtrIn); }

    bool nextDecays(ProcessLevelData& d, Event& process) { return resonanceDecays.next( process);}

    // --- private ---

    // Generate the next event with one interaction.
    bool nextOne(ProcessLevelData& d, Event& process);

    // Generate the next event with two hard interactions.
    bool nextTwo(ProcessLevelData& d, Event& process);

    // Append the second to the first process list.
    void combineProcessRecords(ProcessLevelData& d, Event& process, Event& process2);

    // Check that colours match up.
    bool checkColours(ProcessLevelData& d, Event& process);

    // Print statistics when two hard processes allowed.
    void statistics2(ProcessLevelData& d, bool reset, ostream& os = cout);


  }
