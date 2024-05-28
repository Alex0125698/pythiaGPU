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

    bool startsWith(const std::string& str, const std::string& prefix, bool case_sensitive)
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

    std::string replace(const std::string mainstring, const std::string substring, const std::string replacement, const int num_occurences)
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

    std::string ljust(const std::string to_be_filled, const size_t width, const char fill_character)
    {
      std::string tmp(std::max(to_be_filled.size(),width),fill_character);
      std::copy(to_be_filled.begin(),to_be_filled.end(),tmp.begin());
      return tmp;
    }

    std::string rjust(const std::string to_be_filled, const size_t width, const char fill_character)
    {
      std::string tmp(std::max(to_be_filled.size(),width),fill_character);
      std::copy(to_be_filled.begin(),to_be_filled.end(),tmp.begin()+tmp.size()-to_be_filled.size());
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

    constexpr int nPlaceholders = 10000;
    bool initialized = false;
    int nEvents = 0;
    Data data;

    // ---- class definitions

    void init()
    {
      if (!initialized)
      {
        initialized = true;
        data.nextIndex = 0;
        data.names.resize(nPlaceholders, nullptr);
        data.durations.resize(nPlaceholders);
        data.callCounts.resize(nPlaceholders);
        data.loopCounts_average.resize(nPlaceholders);
        data.loopCounts_min.resize(nPlaceholders);
        data.loopCounts_max.resize(nPlaceholders);
        
        clearData();
      }
    }

    void clearData()
    {
      data.nextIndex = 0;
      // std::fill(data.names.begin(), data.names.end(), nullptr);
      std::fill(data.durations.begin(), data.durations.end(), .0);
      std::fill(data.callCounts.begin(), data.callCounts.end(), 0);
      std::fill(data.loopCounts_average.begin(), data.loopCounts_average.end(), 0);
      std::fill(data.loopCounts_min.begin(), data.loopCounts_min.end(), +std::numeric_limits<double>::max());
      std::fill(data.loopCounts_max.begin(), data.loopCounts_max.end(), -std::numeric_limits<double>::max());
    }

    Data& getTimings()
    {
      return data;
    }

    int getNextUniqueIndex(const char* name)
    {
      if (!initialized)
      {
        std::cerr << "invalid call to getNextUniqueIndex before init" << std::endl;
        exit(0);
        return -1;
      }

      // if (data.nextIndex >= nPlaceholders)
      // {
      //   std::cerr << "too many calls to getNextUniqueIndex" << std::endl;
      //   exit(0);
      //   return -1;
      // }
      
      data.names[data.nextIndex] = name;
      int tmp = data.nextIndex;
      data.nextIndex++;
      return tmp;
    }

    void recordLoopMinMaxAverage()
    {
      int actualSize = data.nextIndex;
      for (int i=0; i<actualSize; i++)
      {
          if (data.durations[i] >= 0) continue; // skip stopwatch
          if (data.durations[i] == -0.1) continue;
          if (data.durations[i] == -0.2) continue;
          if (data.durations[i] == -0.3) continue;

          double ave = data.callCounts[i] / -data.durations[i];
          data.callCounts[i] = 0;
          data.durations[i] = -0.1;
          data.loopCounts_average[i] = (ave + data.loopCounts_average[i] * nEvents) / (nEvents+1);
          data.loopCounts_min[i] = std::min(ave, data.loopCounts_min[i]);
          data.loopCounts_max[i] = std::max(ave, data.loopCounts_max[i]);
      }
      nEvents += 1;
    }

    void cleanUpTimingData()
    {
      // merge duplicates

      int actualSize = data.nextIndex;
      std::vector<int> eraseIndices;

      for (int i=0; i<actualSize; ++i)
      {
        // don't merge endScopes!
        if (data.durations[i] < 0.) continue;
        
        for (int j=i+1; j<actualSize; ++j)
        {
          // don't merge endScopes!
          if (data.durations[i] < 0.) continue;

          if (std::strcmp(data.names[j],data.names[i]) == 0)
          {
            // ???
            // if (data.durations[i] == -0.1 || data.durations[j] == -0.1) continue;
            eraseIndices.push_back(j);
            data.durations[i] += data.durations[j];
            data.callCounts[i] += data.callCounts[j];
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
          data.names[i] = data.names[j];
          data.durations[i] = data.durations[j];
          data.callCounts[i] = data.callCounts[j];
          data.loopCounts_average[i] = data.loopCounts_average[j];
          data.loopCounts_min[i] = data.loopCounts_min[j];
          data.loopCounts_max[i] = data.loopCounts_max[j];
        }

        i++; j++;
      }
      data.nextIndex = i;
    }

    void printTimings()
    {
      std::vector<std::string> pNames;
      std::vector<double> pTimes;
      std::vector<double> sumTimes;
      std::string indent = "";

      double durTotal = data.durations[0];

      std::cerr << "------ task hierarchy ------" << std::endl;
      for (int i=0; i<data.nextIndex; i++)
      {
        auto dur = data.durations[i];
        auto name = data.names[i];
        auto calls = data.callCounts[i];
        
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
            std::cerr << indent << rjust(pNames.back()+"/missing", 40) << ": " << 0 << ", " << dec(perTotal,1) << "%, " << dec(per,1) << "%" << std::endl;
            
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

          if (dur < 0 || data.loopCounts_max[i] > 0 || startsWith(name,"loop") || name2 == "TimeShower::pTnext/loopOverHardPartons")
          {
            std::cerr << indent << "~~" << name2 << ": " << dec(data.loopCounts_average[i],2) << " [" << data.loopCounts_min[i] << ", " << data.loopCounts_max[i] << "] iterations" << std::endl;
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
            std::cerr << indent << name2 << ": " << calls << ", " << dec(perTotal,1) << "%, " << dec(1000.*dur,1) << std::endl;
            // std::cerr << indent << name2 << ": " << calls << ", " << dec(perTotal,1) << "%, " << dec(per,1) << "%" << std::endl;
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

        if (pNames.size() == 1 && i == data.nextIndex -1)
        {
          // before exiting scope, print the missing time
          double dur = pTimes.back() - sumTimes.back();
          double per = 100.* dur / pTimes.back();
          double perTotal = 100.*dur / durTotal;
          std::cerr << indent << rjust(pNames.back()+"/missing", 40, ' ') << ": " << 0 << ", " << dec(perTotal,1) << "%, " << dec(per,1) << "%" << std::endl;
        }        
      }
      std::cerr << "-------------------------" << std::endl;
    }

    bool started()
    {
      return initialized;
    }

    // ---- class definitions

    void PlaceHolder::setUniqueIndex(const int uniqueIndex_)
    {
      uniqueIndex = uniqueIndex_;
    }

    Stopwatch::Stopwatch()
    {
      start();
    }

    Stopwatch::~Stopwatch()
    {
      if (uniqueIndex == -1) return;
      stop();
      data.durations[uniqueIndex] += totalDuration;
      data.callCounts[uniqueIndex] += 1;
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

    LoopCounter::LoopCounter()
    {
      if(uniqueIndex == -1) return;

      if (data.durations[uniqueIndex] == -0.1)
      {
        data.durations[uniqueIndex] = 0.;
      }
    }

    LoopCounter::~LoopCounter()
    {
      if (uniqueIndex == -1) return;
      data.durations[uniqueIndex] -= 1;
      data.callCounts[uniqueIndex] += totalCount;
    }

    void LoopCounter::count()
    {
      totalCount += 1;
    }

    void LoopCounter::finish()
    {
      if (uniqueIndex == -1) return;
      // data.durations[uniqueIndex_] = -0.1;
      // data.callCounts[uniqueIndex_] = 0;
      // data.names[uniqueIndex_] = data.names[uniqueIndex];
    }

    Scope::~Scope()
    {
      if (uniqueIndex == -1) return;
      data.durations[uniqueIndex] = -0.2;
      if (!endScope) return;
      int uniqueIndex2 = getNextUniqueIndex("endScope");
      data.durations[uniqueIndex2] = -0.3;
    }
    
  }

}

// ----------------------------------------------
