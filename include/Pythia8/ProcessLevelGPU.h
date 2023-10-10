



namespace Pythia8 {

  // all functions that access sigmaProcessPtr or phaseSpacePtr are template functions
  // that 

  struct ProcessLevelData 
  {
    // --- inputs ---

    // these are probably constant
    Info*           infoPtr; // Pointer to various information on the generation.
    Settings*       settingsPtr; // Pointer to the settings database.
    ParticleData*   particleDataPtr; // Pointer to the particle data table.
    BeamParticle*   beamAPtr; // Pointers to the two incoming beams.
    BeamParticle*   beamBPtr;
    Couplings*      couplingsPtr; // Pointer to Standard Model couplings, including alphaS and alphaEM.
    SLHAinterface*  slhaInterfacePtr; // Pointer to SusyLesHouches object for interface to SUSY spectra.
    UserHooks*      userHooksPtr; // Pointer to userHooks object for user interaction with program.
    
    SusyLesHouches* slhaPtr; // Pointer to an SLHA object.
    LHAup*          lhaUpPtr; // Pointer to LHAup for generating external events.

    // Todo: does this modify anything ??
    Rndm*           rndmPtr; // Pointer to the random number generator.
    SigmaTotal*     sigmaTotPtr;  // Pointer to SigmaTotal object needed to handle soft QCD processes.
    ResonanceDecays resonanceDecays; // ResonanceDecay object does sequential resonance decays.
    
    // --- private ---

    static const int MAXLOOP;
    bool   doSecondHard, doSameCuts, allHardSame, noneHardSame, // Generic info for process generation.
           someHardSame, cutsAgree, cutsOverlap, doResDecays;
    int    nImpact, startColTag;
    double mHatMin1, mHatMax1, pTHatMin1, pTHatMax1, mHatMin2, mHatMax2,
           pTHatMin2, pTHatMax2, sigmaND, sumImpactFac, sum2ImpactFac;
    int    iContainer, i2Container, iLHACont;
    double sigmaMaxSum, sigma2MaxSum;
    // (I think there is only 1)
    vector<ProcessContainerData*> containerPtrs; // Vector of containers of internally-generated processes.
    // (used for next two only)
    vector<ProcessContainerData*> container2Ptrs; // Ditto for optional choice of a second hard process.
    // (not sure if it is used)
    ProcessContainerData containerLHAdec; // Single half-dummy container for LHA input of resonance decay only.
  };

  struct SetupContainersData 
  {
    // --- private ---
    vector<int> idVecA, idVecB; // Arrays of allowed outgoing SUSY particles and their lengths.
    int nVecA, nVecB;
    SigmaOniaSetup charmonium, bottomonium; // Helper class to setup onia production.
  };

}