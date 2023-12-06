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
  // Benchmark::init();


  int N = 10;
  std::vector<double> t_cons(N,0.0);
  std::vector<double> t_init(N,0.0);
  std::vector<double> t_loop(N,0.0);

  for (int x=0; x<N; ++x)
  {
    Benchmark_start(Everything);
    Benchmark_start(Everything_constructor);

    Benchmark::Stopwatch t1;
    t1.start();

    // Settings (Beam CM energy, turn on all QCD processes, set phaseSpace pT min)
    Pythia pythia;
    pythia.readString("Beams:eCM = 8000.");
    pythia.readString("HardQCD:all = on");
    pythia.readString("PhaseSpace:pTHatMin = 20.");

    t1.stop();
    t_cons[x] = t1.totalDuration;
    
    Benchmark_stop(Everything_constructor);
    Benchmark_start(Everything_init);
    
    Benchmark::Stopwatch t2;
    t2.start();

    // initialization
    pythia.init();
    
    t2.stop();
    t_init[x] = t2.totalDuration;
    
    
    Benchmark_stop(Everything_init);
    Benchmark_start(Everything_eventLoop);

    // Hist mult("Number of charged particles", 100, -0.5, 799.5);
    
    Benchmark::Stopwatch t3;
    t3.start();

    // Begin event loop
    for (int iEvent = 0; iEvent < 200; ++iEvent) 
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
      
      // // add to histogram
      // mult.fill( nCharged );

      Benchmark::recordLoopMinMaxAverage();
    }
    
    t3.stop();
    t_loop[x] = t3.totalDuration;
    

    // Statistics. What is it?
    pythia.stat();

    // // print histogram to terminal
    // cout << mult;

    Benchmark_stop(Everything_eventLoop);
    Benchmark_stop(Everything);
    
  }
  
  // Benchmark::printTimings();

  double min_cons=1e99, max_cons=0., ave_cons=0.;
  double min_init=1e99, max_init=0., ave_init=0.;
  double min_loop=1e99, max_loop=0., ave_loop=0.;

  for (int x=0; x<N; ++x)
  {
    min_cons = std::min(min_cons, t_cons[x]);
    max_cons = std::max(max_cons, t_cons[x]);
    ave_cons += t_cons[x];
    min_init = std::min(min_init, t_init[x]);
    max_init = std::max(max_init, t_init[x]);
    ave_init += t_init[x];
    min_loop = std::min(min_loop, t_loop[x]);
    max_loop = std::max(max_loop, t_loop[x]);
    ave_loop += t_loop[x];
  }

  std::cerr << "cons " << std::fixed << std::setprecision(3) << min_cons << " " << ave_cons/N << " " << max_cons << std::endl;
  std::cerr << "init " << std::fixed << std::setprecision(3) << min_init << " " << ave_init/N << " " << max_init << std::endl;
  std::cerr << "loop " << std::fixed << std::setprecision(3) << min_loop << " " << ave_loop/N << " " << max_loop << std::endl;

  return 0;
}
