#ifndef PYTHIA_STATE_H
#define PYTHIA_STATE_H

#include "Pythia8/forwardDeclarations.h"
#include "Pythia8/Settings.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/Info.h"
#include "Pythia8/Basics.h"
#include "Pythia8/SigmaTotal.h"
#include "Pythia8/SLHAinterface.h"

namespace Pythia8 {

struct PythiaState
{
  Info info;
  Settings settings; // true singleton
  ParticleData particleData; // true singleton
  Rndm rndm;
  Couplings* couplings = nullptr;
  SLHAinterface slhaInterface;
  BeamParticle* beamA = nullptr;
  BeamParticle* beamB = nullptr;
  SigmaTotal sigmaTotal;
  LHAup* lhaUp = nullptr;
};

}

#endif // PYTHIA_STATE_H
