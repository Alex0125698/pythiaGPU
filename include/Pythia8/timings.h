

#ifndef TIMINGS_H
#define TIMINGS_H

#include <chrono>
#include <vector>
#include <string>

namespace Pythia8 
{
  namespace Benchmark
  {
    class PlaceHolder
    {
      public:

        int uniqueIndex = -1;

      public:
      
        void setUniqueIndex(int uniqueIndex_, const char* name_);
    };

    class Stopwatch : public PlaceHolder
    {
      public:

        std::chrono::high_resolution_clock::time_point startTime;
        double totalDuration = 0.0;
        bool running = false;

      public:

        ~Stopwatch();
        void reset();
        void start();
        void stop();
    };

    class LoopCounter : public PlaceHolder
    {
      public:

        int totalCount = 0;
      
      public:

        ~LoopCounter();
        void count();
        void finish(int uniqueIndex_);
    };

    bool started();
    int getUniqueIndex();
    void init();
    void printTimings();
    void recordLoopMinMaxAverage();

  }

}

#define Benchmark_placeholder(name)             \
Benchmark::PlaceHolder pl_##name;               \
if (Benchmark::started())                       \
{                                               \
  static int ind = Benchmark::getUniqueIndex(); \
  pl_##name.setUniqueIndex(ind,#name);          \
}

#define Benchmark_start(name)                   \
Benchmark::Stopwatch sw_##name;                 \
if (Benchmark::started())                       \
{                                               \
  static int ind = Benchmark::getUniqueIndex(); \
  sw_##name.setUniqueIndex(ind,#name);          \
  sw_##name.start();                            \
}

#define Benchmark_stop(name) \
sw_##name.stop();

#define Benchmark_loopStart(name)               \
Benchmark::LoopCounter lc_##name;               \
if (Benchmark::started())                       \
{                                               \
  static int ind = Benchmark::getUniqueIndex(); \
  lc_##name.setUniqueIndex(ind,#name);          \
}

#define Benchmark_loopCount(name) \
lc_##name.count();

#define Benchmark_loopStop(name)                \
if (Benchmark::started())                       \
{                                               \
  static int ind = Benchmark::getUniqueIndex(); \
  lc_##name.finish(ind);                        \
}

#endif // TIMINGS_H
