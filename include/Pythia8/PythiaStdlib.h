// PythiaStdlib.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for functionality pulled in from Stdlib,
// plus a few useful utilities (small powers; positive square root,
// Gamma function).

#ifndef Pythia8_PythiaStdlib_H
#define Pythia8_PythiaStdlib_H

// make timings avaliable everywhere
#include "Pythia8/timings.h"

// Stdlib header files for mathematics.
#include <cmath>
#include <cstdlib>
#include <algorithm>

// Stdlib header files for strings and containers.
#include <string>
#include <vector>
#include <map>
#include <deque>
#include <set>
#include <list>

// Stdlib header file for dynamic library loading.
#ifndef _MSC_VER
#define dlsym __
#include <dlfcn.h>
#undef dlsym
#endif

// Redefine dlsym to suppress compiler warnings.
extern "C" void *(*dlsym(void *handle, const char *symbol))();

// Stdlib header file for input and output.
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "Pythia8/SIMDString.h"

// Define pi if not yet done.
#ifndef M_PI
#define M_PI 3.1415926535897932385
#endif

// By this declaration you do not need to use std:: qualifier everywhere.
//using namespace std;

// Alternatively you can specify exactly which std:: methods will be used.
// Now made default so std does not spill outside namespace Pythia8.
namespace Pythia8 {


// Generic utilities and mathematical functions.
using std::swap;
using std::max;
using std::min;
using std::abs;
using std::sort;

// Strings and containers.
using std::pair;
using std::make_pair;
using std::string;
using std::vector;
using std::map;
using std::multimap;
using std::deque;
using std::set;

// using cstring = const char*;
// using string = SIMDString<64>;

// Input/output streams.
using std::cin;
using std::cout;
using std::cerr;
using std::istream;
using std::ostream;
using std::fstream;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::istringstream;
using std::ostringstream;
using std::ios;

// Input/output formatting.
using std::endl;
using std::fixed;
using std::scientific;
using std::left;
using std::right;
using std::setw;
using std::setprecision;

} // end namespace Pythia8

namespace Pythia8 {

inline void toLowerInPlace(string& name)
{
  for (int i=0; i<(int)name.size(); ++i)
  {
    constexpr char delta = 'a' - 'A';
    bool swap = name[i] >= 'A' && name[i] <= 'Z';
    name[i] = swap ? name[i] + delta : name[i];
  }
}

const string& toLower(const string& name)
{
  // for some reason we also trim whitespace characters
  // will assume no non-printable characters are being used at all

  int tstart = 0, tend = name.size()-1;
  while (name[tstart] <= ' ' && tstart < name.size()) ++tstart;
  while (name[tend] <= ' ' && tend >= 0) --tend;

  // figure out where the non-lowercase begins
  int start = tstart;
  while ((name[start] < 'A' || name[start] > 'Z') && start <= tend) ++start;

  // if no whitespace or uppercase then return original string
  if (start >= name.size()) return name;

  // prepare buffer to write new string
  thread_local string tmp;
  if (tmp.size() < name.size()) tmp = name;
  tmp.resize(tend-tstart+1);

  // write everything after whitespace and before first uppercase
  for (int i=tstart; i<start; ++i)
  {
    tmp[i-tstart] = name[i];
  }

  // write rest up to start of trailing whitespace
  for (int i=start; i<=tend; ++i)
  {
    constexpr char delta = 'a' - 'A';
    bool swap = name[i] >= 'A' && name[i] <= 'Z';
    tmp[i-tstart] = swap ? name[i] + delta : name[i];
  }

  return tmp;
}

// Powers of small integers - for balance speed/code clarity.
inline double pow2(const double& x) {return x*x;}
inline double pow3(const double& x) {return x*x*x;}
inline double pow4(const double& x) {return x*x*x*x;}
inline double pow5(const double& x) {return x*x*x*x*x;}
inline double pow6(const double& x) {return x*x*x*x*x*x;}

// Avoid problem with negative square root argument (from roundoff).
inline double sqrtpos(const double& x) {return sqrt( max( 0., x));}

// The Gamma function for real argument.
double GammaReal(double x);

  using cstring = string;
  // using cstring = std::string_view;
  // using cstring = const char*;
  // using stringref = const std::string_view;
  using stringref = const std::string&;
  using stringref2 = const std::string&;
  using stringbuf = std::string;
  // using stringbuf = std::string;

} // end namespace Pythia8

#endif // Pythia8_PythiaStdlib_H
