// BeamShape.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the BeamShape class.

#include "Pythia8/BeamShape.h"

namespace Pythia8 {

//==========================================================================

// The BeamShape class.

//--------------------------------------------------------------------------

// Initialize beam parameters.

  void BeamShape::init( Settings& settings, Rndm* rndmPtrIn) {

  // Save pointer.
  rndmPtr             = rndmPtrIn;

  // Main flags.
  allowMomentumSpread = settings.get(Flag::Beams_allowMomentumSpread);
  allowVertexSpread   = settings.get(Flag::Beams_allowVertexSpread);

  // Parameters for beam A momentum spread.
  sigmaPxA            = settings.get(Param::Beams_sigmaPxA);
  sigmaPyA            = settings.get(Param::Beams_sigmaPyA);
  sigmaPzA            = settings.get(Param::Beams_sigmaPzA);
  maxDevA             = settings.get(Param::Beams_maxDevA);

  // Parameters for beam B momentum spread.
  sigmaPxB            = settings.get(Param::Beams_sigmaPxB);
  sigmaPyB            = settings.get(Param::Beams_sigmaPyB);
  sigmaPzB            = settings.get(Param::Beams_sigmaPzB);
  maxDevB             = settings.get(Param::Beams_maxDevB);

  // Parameters for beam vertex spread.
  sigmaVertexX        = settings.get(Param::Beams_sigmaVertexX);
  sigmaVertexY        = settings.get(Param::Beams_sigmaVertexY);
  sigmaVertexZ        = settings.get(Param::Beams_sigmaVertexZ);
  maxDevVertex        = settings.get(Param::Beams_maxDevVertex);
  sigmaTime           = settings.get(Param::Beams_sigmaTime);
  maxDevTime          = settings.get(Param::Beams_maxDevTime);

  // Parameters for beam vertex offset.
  offsetX             = settings.get(Param::Beams_offsetVertexX);
  offsetY             = settings.get(Param::Beams_offsetVertexY);
  offsetZ             = settings.get(Param::Beams_offsetVertexZ);
  offsetT             = settings.get(Param::Beams_offsetTime);

}

//--------------------------------------------------------------------------

// Set the two beam momentum deviations and the beam vertex.

void BeamShape::pick() {

  // Reset all values.
  deltaPxA = deltaPyA = deltaPzA = deltaPxB = deltaPyB = deltaPzB
    = vertexX = vertexY = vertexZ = vertexT = 0.;

  // Set beam A momentum deviation by a three-dimensional Gaussian.
  if (allowMomentumSpread) {
    double totalDev, gauss;
    do {
      totalDev = 0.;
      if (sigmaPxA > 0.) {
        gauss     = rndmPtr->gauss();
        deltaPxA  = sigmaPxA * gauss;
        totalDev += gauss * gauss;
      }
      if (sigmaPyA > 0.) {
        gauss     = rndmPtr->gauss();
        deltaPyA  = sigmaPyA * gauss;
        totalDev += gauss * gauss;
      }
      if (sigmaPzA > 0.) {
        gauss     = rndmPtr->gauss();
        deltaPzA  = sigmaPzA * gauss;
        totalDev += gauss * gauss;
      }
    } while (totalDev > maxDevA * maxDevA);

    // Set beam B momentum deviation by a three-dimensional Gaussian.
    do {
      totalDev = 0.;
      if (sigmaPxB > 0.) {
        gauss     = rndmPtr->gauss();
        deltaPxB  = sigmaPxB * gauss;
        totalDev += gauss * gauss;
      }
      if (sigmaPyB > 0.) {
        gauss     = rndmPtr->gauss();
        deltaPyB  = sigmaPyB * gauss;
        totalDev += gauss * gauss;
      }
      if (sigmaPzB > 0.) {
        gauss     = rndmPtr->gauss();
        deltaPzB  = sigmaPzB * gauss;
        totalDev += gauss * gauss;
      }
    } while (totalDev > maxDevB * maxDevB);
  }

  // Set beam vertex location by a three-dimensional Gaussian.
  if (allowVertexSpread) {
    double totalDev, gauss;
    do {
      totalDev = 0.;
      if (sigmaVertexX > 0.) {
        gauss     = rndmPtr->gauss();
        vertexX   = sigmaVertexX * gauss;
        totalDev += gauss * gauss;
      }
      if (sigmaVertexY > 0.) {
        gauss     = rndmPtr->gauss();
        vertexY   = sigmaVertexY * gauss;
        totalDev += gauss * gauss;
      }
      if (sigmaVertexZ > 0.) {
        gauss     = rndmPtr->gauss();
        vertexZ   = sigmaVertexZ * gauss;
        totalDev += gauss * gauss;
      }
    } while (totalDev > maxDevVertex * maxDevVertex);

    // Set beam collision time by a Gaussian.
    if (sigmaTime > 0.) {
      do gauss    = rndmPtr->gauss();
      while (abs(gauss) > maxDevTime);
      vertexT     = sigmaTime * gauss;
    }

    // Add offset to beam vertex.
    vertexX      += offsetX;
    vertexY      += offsetY;
    vertexZ      += offsetZ;
    vertexT      += offsetT;
  }

}

//==========================================================================

} // end namespace Pythia8
