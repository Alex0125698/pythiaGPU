// ----------------------------------------------

// timings source

#include "Pythia8/timings.h"

#include <sstream>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <algorithm>

namespace Pythia8 
{
  namespace Benchmark
  {
    // ---- utils

    bool startsWith(const std::string& str, const std::string& prefix, bool case_sensitive = true)
    {
      if (&prefix == &str) return true; // str and prefix are the same string
      if (prefix.length() > str.length()) return false;
      if (prefix.empty()) return false;
      for (size_t i = 0; i < prefix.length(); ++i) {
          if(case_sensitive)
          {
             if (prefix[i] != str[i]) return false;
          }
          else
          {
             if (tolower(prefix[i]) != tolower(str[i])) return false;
          }
      }
      return true;
    }

    std::string replace(const std::string mainstring, const std::string substring, const std::string replacement, const int num_occurences = -1)
    {
      if (mainstring.empty() || substring.empty() || substring == replacement || num_occurences == 0) 
        return mainstring;
      
      std::string result; 
      result.reserve(mainstring.size()*2);

      int remaining_replacements = num_occurences;
      size_t index = 0;

      while(remaining_replacements != 0)
      {
        --remaining_replacements;

        // get index of next substring
        size_t next_occurence = mainstring.find(substring,index);
        if (next_occurence == std::string::npos) break;

        // copy in everything before the substring
        std::copy(mainstring.begin()+index,mainstring.begin()+next_occurence,std::back_inserter(result));

        // copy in the replacement string
        std::copy(replacement.begin(),replacement.end(),std::back_inserter(result));

        // move the index past occurence
        index = next_occurence + substring.size();
      }

      // copy in the rest of the string
      std::copy(mainstring.begin()+index,mainstring.end(),std::back_inserter(result));

      return result;
    }

    std::string ljust(const std::string to_be_filled, const size_t width, const char fill_character = ' ')
    {
      std::string tmp(std::max(to_be_filled.size(),width),fill_character);
      std::copy(to_be_filled.begin(),to_be_filled.end(),tmp.begin());
      return tmp;
    }

    std::string rjust(const std::string to_be_filled, const size_t width, const char fill_character = ' ')
    {
      std::string tmp(std::max(to_be_filled.size(),width),fill_character);
      std::copy(to_be_filled.end()-to_be_filled.size(),to_be_filled.end(),tmp.begin());
      return tmp;
    }

    std::vector<size_t> argSort(const std::vector<std::string>& vec)
    {
      std::vector<std::pair<size_t,std::string>> tmp(vec.size());

      for (size_t i=0; i<vec.size(); ++i)
      {
        tmp[i].first = i;
        tmp[i].second = vec[i];
      }

      std::sort(tmp.begin(), tmp.end(), [&](const std::pair<size_t,std::string>& a, const std::pair<size_t,std::string>& b) { return a.second < b.second;});


      std::vector<size_t> result(vec.size());

      for (size_t i=0; i<vec.size(); ++i)
      {
        result[i] = tmp[i].first;
      }

      return result;
    };

    void permutate(std::vector<std::string>& arr, const std::vector<size_t>& indices)
    {
      auto arr2 = arr;
      for (size_t i=0; i<indices.size(); ++i)
        arr[i] = arr2[indices[i]];
    }

    void permutate(std::vector<int>& arr, const std::vector<size_t>& indices)
    {
      auto arr2 = arr;
      for (size_t i=0; i<indices.size(); ++i)
        arr[i] = arr2[indices[i]];
    }

    bool in(const int x, const std::vector<int>& v)
    {
      for (auto& xx : v) if (x == xx) return true;
      return false;
    };

    std::string dec(const double val, const int nDecimals)
    {
      std::stringstream s;
      s << std::fixed << std::setprecision(nDecimals) << val;
      return s.str();
    }

    // ---- actual benchmark code

    int nextIndex = 0;
    bool benchmark_started = false;
    constexpr int nPlaceholders = 10000;
    std::vector<const char*> names;
    std::vector<double> durations;
    std::vector<int> callCounts;

    std::vector<double> loopCounts_average;
    std::vector<double> loopCounts_min;
    std::vector<double> loopCounts_max;
    int nEvents = 0;

    void recordLoopMinMaxAverage()
    {
      int actualSize = nextIndex;
      for (int i=0; i<actualSize; i++)
      {
          if (durations[i] >= 0) continue;
          if (durations[i] == -0.1) continue;

          double ave = callCounts[i] / -durations[i];
          callCounts[i] = 0;
          durations[i] = 0.;
          loopCounts_average[i] = (ave + loopCounts_average[i] * nEvents) / (nEvents+1);
          loopCounts_min[i] = std::min(ave, loopCounts_min[i]);
          loopCounts_max[i] = std::max(ave, loopCounts_max[i]);
      }
      nEvents += 1;
    }

    bool started()
    {
      return benchmark_started;
    }
    
    int getUniqueIndex()
    {
        if (!benchmark_started)
        {
          std::cerr << "invalid call to getUniqueIndex before init" << std::endl;
          exit(0);
          return -1;
        }
        if (nextIndex >= nPlaceholders-2)
        {
          std::cerr << "too many calls to getUniqueIndex" << std::endl;
          exit(0);
          return -1;
        }
        return nextIndex++;
    }

    void init()
    {
      if (!benchmark_started)
      {
        benchmark_started = true;
        names.resize(nPlaceholders);
        durations.resize(nPlaceholders, 0.0);
        callCounts.resize(nPlaceholders, 0);
        loopCounts_average.resize(nPlaceholders, 0);
        loopCounts_min.resize(nPlaceholders, +std::numeric_limits<double>::max());
        loopCounts_max.resize(nPlaceholders, -std::numeric_limits<double>::max());
      }
    }

    void printTimings()
    {
      // merge duplicates

      int actualSize = nextIndex;
      std::vector<int> eraseIndices;

      for (int i=0; i<actualSize; ++i)
      {
        for (int j=i+1; j<actualSize; ++j)
        {
          if (std::strcmp(names[j],names[i]) == 0)
          {
            if (durations[i] == -0.1 || durations[j] == -0.1) continue;
            eraseIndices.push_back(j);
            durations[i] += durations[j];
            callCounts[i] += callCounts[j];
          }
        }
      }

      // erase duplicates

      int i = 0, j = 0;
      while (true)
      {
        while (in(j,eraseIndices)) j++;

        if (j >= actualSize) break;
        
        if (i != j)
        {
          names[i] = names[j];
          durations[i] = durations[j];
          callCounts[i] = callCounts[j];
          loopCounts_average[i] = loopCounts_average[j];
          loopCounts_min[i] = loopCounts_min[j];
          loopCounts_max[i] = loopCounts_max[j];
        }

        i++; j++;
      }
      actualSize = i;

      std::vector<std::string> pNames;
      std::vector<double> pTimes;
      std::vector<double> sumTimes;
      std::string indent = "";

      double durTotal = durations[0];

      std::cerr << "------ task hierarchy ------" << std::endl;
      for (int i=0; i<actualSize; i++)
      {
        auto dur = durations[i];
        auto name = names[i];
        auto calls = callCounts[i];
        
        if (dur == -0.1)
        {
          indent.resize(std::max((size_t)0,indent.size()-2));
          continue;
        }

        std::string pName = "";
        double pTime = 0.0;

        if (pNames.size() >= 2)
        {
          std::string ppName = *(pNames.end()-2);
          if (startsWith(name,ppName))
          {
            // before exiting scope, print the missing time
            double dur = pTimes.back() - sumTimes.back();
            double per = 100.* dur / pTimes.back();
            double perTotal = 100.*dur / durTotal;
            std::cerr << indent << rjust(pNames.back()+"/missing", 40) << ": " << 0 << ", " << dec(perTotal,1) << "%, " << dec(per,1) << '%' << std::endl;
            
            sumTimes.pop_back();
            pNames.pop_back();
            pTimes.pop_back();
            indent.resize(std::max((size_t)0,indent.size()-2));
          }
        }

        if (pNames.size() >= 1)
        {
          pName = pNames.back();
          pTime = pTimes.back();
        }
        
        std::string name2 = name;
        if (pName != "") name2 = replace(name2, pName+"_", pName+"/");
        name2 = replace(name2, "0", "::");
        name2 = rjust(name2, 40);
        double perTotal = 100.*dur / durTotal;

        if (startsWith(name,pName))
        {
          double per = 100.* dur / pTime;

          if (dur < 0 || loopCounts_max[i] > 0 || startsWith(name,"loop") || name2 == "TimeShower::pTnext/loopOverHardPartons")
          {
            std::cerr << indent << "~~" << name2 << ": " << dec(loopCounts_average[i],2) << " [" << loopCounts_min[i] << ", " << loopCounts_max[i] << "] iterations" << std::endl;
            indent += "  ";
          }
          else
          {
            if (sumTimes.empty())
            {
              std::cerr << "invalid" << std::endl;
              exit(0);
            }
            sumTimes.back() += dur;
                        double perTotal = 100.*dur / durTotal;
            std::cerr << indent << name2 << ": " << calls << ", " << dec(perTotal,1) << "%, " << dec(per,1) << '%' << std::endl;
          }
        }
        else
        {
          sumTimes.push_back(0.0);
          pNames.push_back(name);
          pTimes.push_back(dur);
          std::cerr << indent << name2 << ": " << calls << ", " << dec(perTotal,1) << "%, " << std::endl;
          indent += "  ";
        }

        if (pNames.size() == 1 && i == actualSize -1)
        {
          // before exiting scope, print the missing time
          double dur = pTimes.back() - sumTimes.back();
          double per = 100.* dur / pTimes.back();
          double perTotal = 100.*dur / durTotal;
          std::cerr << indent << rjust(pNames.back()+"/missing", 40, ' ') << ": " << 0 << ", " << dec(perTotal,1) << "%, " << dec(per,1) << '%' << std::endl;
        }        
      }
      std::cerr << "-------------------------" << std::endl;
    }
  
    // ---- class definitions

    void PlaceHolder::setUniqueIndex(int uniqueIndex_, const char* name_)
    {
      uniqueIndex = uniqueIndex_;
      Benchmark::names[uniqueIndex_] = name_;
    }

    Stopwatch::~Stopwatch()
    {
      if (uniqueIndex < 0) return;
      stop();
      Benchmark::durations[uniqueIndex] += totalDuration;
      Benchmark::callCounts[uniqueIndex] += 1;
    }

    void Stopwatch::reset()
    {
      totalDuration = 0.0;
      running = false;
    }

    void Stopwatch::start()
    {
      running = true;
      startTime = std::chrono::high_resolution_clock::now();
    }

    void Stopwatch::stop()
    {
      if (running)
      {
        running = false;
        auto endTime = std::chrono::high_resolution_clock::now();
        totalDuration += std::chrono::duration<double>(endTime-startTime).count();
      }
    }

    LoopCounter::~LoopCounter()
    {
      if (uniqueIndex < 0) return;
      Benchmark::durations[uniqueIndex] -= 1;
      Benchmark::callCounts[uniqueIndex] += totalCount;
    }

    void LoopCounter::count()
    {
      totalCount += 1;
    }

    void LoopCounter::finish(int uniqueIndex_)
    {
      Benchmark::durations[uniqueIndex_] = -0.1;
      Benchmark::callCounts[uniqueIndex_] = 0;
      Benchmark::names[uniqueIndex_] = Benchmark::names[uniqueIndex];
    }

  }

}

// ----------------------------------------------
