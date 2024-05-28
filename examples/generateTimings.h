#ifndef GENERATE_TIMINGS_H
#define GENERATE_TIMINGS_H

#include <vector>
#include <iostream>
#include <limits>
#include <cstring>
#include "Pythia8/timings.h"

using namespace Pythia8;
using Pythia8::Benchmark::dec;
using Pythia8::Benchmark::rjust;
using Pythia8::Benchmark::ljust;

constexpr double perScopeLimit = 0.5;
constexpr double perProgramLimit = 0.05;

inline std::string format(double val, int nDecimals=0)
{
  return rjust(dec(val,nDecimals),5);
}

inline void printTimings(std::vector<Benchmark::Data>& data)
{
  int nRuns = data.size();
  if (nRuns == 0)
  {
    std::cerr << "no runs" << std::endl;
    exit(0);
  }
  int nBenchmarks = data[0].nextIndex;
  if (nBenchmarks == 0)
  {
    std::cerr << "no benchmarks" << std::endl;
    exit(0);
  }

  std::cerr << "%All %Scope durAve durMin durMax " << std::endl;
  std::cerr << "-------------------------------- \n" << std::endl;

  std::vector<const char*> scopes;
  std::vector<double> durTotal;
  std::vector<double> durSummed;
  std::string prefix = "";
  double durAllPythia = -1;
  
  for (int i=0; i<nBenchmarks; ++i)
  {
    // parse begin-scope
    if (data[0].durations[i] == -0.2)
    {
      // skip begin followed by end
      if (i+1<nBenchmarks)
      {
        if(data[0].durations[i+1] == -0.3)
        {
          ++i; continue;
        }
      }

      scopes.push_back(data[0].names[i]);
      durTotal.push_back(-1.);
      durSummed.push_back(0.0);
      std::cerr << prefix << "> " << scopes.back() << std::endl;
      prefix += "  ";
      continue;
    }
    
    // parse end-scope
    if (data[0].durations[i] == -0.3)
    {
      if (scopes.empty())
      {
        std::cerr << "too many end-scopes" << std::endl;
        exit(0);
      }

      if (durTotal.back() != -1)
      {
        double durMiss = durTotal.back() - durSummed.back();
        double per = (100.*durMiss) / durTotal.back();
        double per2 = (100.*durMiss) / durAllPythia;
        if (durAllPythia == -1.) per2 = 0.;
        if (per > perScopeLimit && per2 > perProgramLimit)
        {
          std::cerr << prefix << ljust("missing",25) << " " << format(per2,1) << "% " << format(per,1) << "% " << format(durMiss) << ' ' << std::endl;
        }
      }

      prefix.resize(std::max((size_t)0,prefix.size()-2));
      std::cerr << prefix << "< " << scopes.back() << std::endl;
      scopes.pop_back();
      durTotal.pop_back();
      durSummed.pop_back();
      continue;
    }

    // parse loop
    if (data[0].durations[i] == -0.1)
    {
      double ave = data[0].loopCounts_average[i];
      int min = data[0].loopCounts_min[i];
      int max = data[0].loopCounts_max[i];
      std::cerr << prefix << ljust(data[0].names[i],25) << " " << format(ave,1) << " " << min << " " << max << ' ' << std::endl;
      continue;
    }

    // parse stopwatch
    if (data[0].durations[i] >= 0.0)
    {
      double durSum = 0.0;
      double durMin = +std::numeric_limits<double>::max();
      double durMax = -std::numeric_limits<double>::max();

      for (int j=0; j<nRuns; ++j)
      {
        if (j > 0 && data[j].nextIndex != 0)
        {
          std::cerr << "inconsistant number of benchmarks" << std::endl;
          exit(0);
        }

        double dur = 1000.*data[j].durations[i];
        durSum += dur;
        durMin = std::min(durMin,dur);
        durMax = std::max(durMax,dur);

      }

      durSum /= nRuns;

      if (std::strcmp(data[0].names[i], "all") == 0)
      {
        durTotal.back() = durSum;
        if (durAllPythia == -1) durAllPythia = durSum;
        // continue;
      }
      durSummed.back() += durSum;
      double per = (100.*durSum) / durTotal.back();
      if (durTotal.back() == -1.) per = 0.;
      double per2 = (100.*durSum) / durAllPythia;
      if (durAllPythia == -1.) per2 = 0.;

      int calls = data[0].callCounts[i];

      if (per > perScopeLimit && per2 > perProgramLimit)
      {
        std::cerr << prefix << ljust(data[0].names[i],25) << " " << format(per2,1) << "% " << format(per,1) << "% " << format(durSum) << " " << format(durMin) << " " << format(durMax) << " " << calls << ' ' << std::endl;
      }

      continue;
    }

    std::cerr << "ERROR " << data[0].names[i] << std::endl;
  }

}


typedef int (*mainfcn_t)();


inline int setupTimings(mainfcn_t fcn, int nIter)
{
  Benchmark::init();
  std::vector<Benchmark::Data> tdata;

  for (int i=0; i<nIter; ++i)
  {
    if (i>0) std::cout.setstate(std::ios_base::failbit);
    Benchmark::clearData();
    {
      Benchmark_function(everything)
      fcn();
    }
    tdata.push_back(Benchmark::getTimings());
  }
  std::cout.clear();

  printTimings(tdata);
  return 0;
}

#endif // GENERATE_TIMINGS_H
