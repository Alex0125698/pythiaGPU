// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the charged multiplicity distribution at the LHC.

#include "Pythia8/Pythia.h"
#include "Pythia8/timings.h"
using namespace Pythia8;

int main() 
{
  Benchmark::init();
  {

    Benchmark_start(Everything);
    Benchmark_start(Everything_constructor);

    // Settings (Beam CM energy, turn on all QCD processes, set phaseSpace pT min)
    Pythia pythia;
    pythia.readString("Beams:eCM = 8000.");
    pythia.readString("HardQCD:all = on");
    pythia.readString("PhaseSpace:pTHatMin = 20.");
    
    Benchmark_stop(Everything_constructor);
    Benchmark_start(Everything_init);
    
    // initialization
    pythia.init();
    
    Benchmark_stop(Everything_init);
    Benchmark_start(Everything_eventLoop);

    Hist mult("Number of charged particles", 100, -0.5, 799.5);
    
    // Begin event loop
    for (int iEvent = 0; iEvent < 100; ++iEvent) 
    {
      // ignore failed events
      if (!pythia.next()) continue;
      
      // Count number of final charged particles
      int nCharged = 0;
      for (int i = 0; i < pythia.event.size(); ++i)
      {
        if (pythia.event[i].isFinal() && pythia.event[i].isCharged())
        {
          ++nCharged;
        }
      }
      
      // add to histogram
      mult.fill( nCharged );

      Benchmark::recordLoopMinMaxAverage();
    }
    
    // Statistics. What is it?
    pythia.stat();

    // print histogram to terminal
    cout << mult;

    Benchmark_stop(Everything_eventLoop);
    Benchmark_stop(Everything);
    
  }
  
  Benchmark::printTimings();

  return 0;
}
