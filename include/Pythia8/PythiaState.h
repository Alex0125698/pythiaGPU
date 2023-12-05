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
  // true singletons
  Settings settings;
  ParticleData particleData;
  // duplicated per pythia event generator
  Info info;
  Rndm rndm;
  BeamParticle* beamA = nullptr;
  BeamParticle* beamB = nullptr;
  Couplings* couplings = nullptr;
  SigmaTotal sigmaTotal;
  SLHAinterface slhaInterface;
  LHAup* lhaUp = nullptr;
};

}
