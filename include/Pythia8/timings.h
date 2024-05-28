

#ifndef TIMINGS_H
#define TIMINGS_H

#include <chrono>
#include <vector>
#include <string>

namespace Pythia8 
{
  namespace Benchmark
  {
    bool startsWith(const std::string& str, const std::string& prefix, bool case_sensitive = true);
    std::string replace(const std::string mainstring, const std::string substring, const std::string replacement, const int num_occurences = -1);
    std::string ljust(const std::string to_be_filled, const size_t width, const char fill_character = ' ');
    std::string rjust(const std::string to_be_filled, const size_t width, const char fill_character = ' ');
    std::vector<size_t> argSort(const std::vector<std::string>& vec);
    void permutate(std::vector<std::string>& arr, const std::vector<size_t>& indices);
    void permutate(std::vector<int>& arr, const std::vector<size_t>& indices);
    bool in(const int x, const std::vector<int>& v);
    std::string dec(const double val, const int nDecimals);

    struct PlaceHolder
    {
      void setUniqueIndex(const int uniqueIndex);
      short uniqueIndex = -1;
    };

    struct Stopwatch : public PlaceHolder
    {
      Stopwatch();
      ~Stopwatch();
      void reset();
      void start();
      void stop();
      std::chrono::high_resolution_clock::time_point startTime;
      float totalDuration = 0.0;
      bool running = false;
    };

    struct LoopCounter : public PlaceHolder
    {
        LoopCounter();
        ~LoopCounter();
        void count();
        void finish();
        unsigned short totalCount = 0;
    };

    struct Scope : public PlaceHolder
    {
      bool endScope = false;
      ~Scope();
    };

    struct Data
    {
      int nextIndex;
      std::vector<const char*> names;
      std::vector<double> durations;
      std::vector<int> callCounts;
      std::vector<double> loopCounts_average;
      std::vector<double> loopCounts_min;
      std::vector<double> loopCounts_max;
    };
    
    bool started();
    int getNextUniqueIndex(const char* name);
    void init();
    void recordLoopMinMaxAverage();
    void clearData();
    Data& getTimings();
  }
}

#define ENABLE_BENCHMARKS

#ifdef ENABLE_BENCHMARKS

#define Benchmark_placeholder(name)             \
Benchmark::PlaceHolder pl_##name;          \
if (Benchmark::started())                       \
{                                               \
  static int ind = Benchmark::getNextUniqueIndex(#name);   \
  pl_##name.setUniqueIndex(ind);                \
}

#define Benchmark_start(name)                   \
Benchmark::Stopwatch sw_##name;            \
if (Benchmark::started())                       \
{                                               \
  static int ind = Benchmark::getNextUniqueIndex(#name);   \
  sw_##name.setUniqueIndex(ind);                \
}

#define Benchmark_stop(name)                    \
sw_##name.stop();

#define Benchmark_loopStart(name)               \
{Benchmark::Scope lsc_##name;                   \
Benchmark::LoopCounter lc_##name;               \
if (Benchmark::started())                       \
{                                               \
  static int ind = Benchmark::getNextUniqueIndex("LOOP");   \
  lsc_##name.setUniqueIndex(ind);                \
  static int ind2 = Benchmark::getNextUniqueIndex(#name);   \
  lc_##name.setUniqueIndex(ind2);                \
  static bool endScope = true;                  \
  if (endScope)                                 \
  {                                             \
    endScope = false;                           \
    lsc_##name.endScope = true;                           \
  }                                             \
}

#define Benchmark_loopCount(name)               \
lc_##name.count();                              

#define Benchmark_loopStop(name)                \
}                          

#define Benchmark_scope(name)                   \
Benchmark::Scope sc_##name;                \
if (Benchmark::started())                       \
{                                               \
  static int ind = Benchmark::getNextUniqueIndex(#name);   \
  sc_##name.setUniqueIndex(ind);                \
  static bool endScope = true;                  \
  if (endScope)                                 \
  {                                             \
    endScope = false;                           \
    sc_##name.endScope = true;                           \
  }                                             \
}

#define Benchmark_function(name)               \
Benchmark_scope(name)                          \
Benchmark_start(all)

#else

#define Benchmark_placeholder(name) ;
#define Benchmark_start(name) ;
#define Benchmark_stop(name) ;
#define Benchmark_loopStart(name) ;
#define Benchmark_loopCount(name) ;
#define Benchmark_loopStop(name) ;

#endif


#endif // TIMINGS_H
