// SigmaSUSY.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// Main authors of this file: N. Desai, P. Skands
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for Supersymmetric process differential cross sections.
// Contains structes derived from SigmaProcess via Sigma2Process.

#ifndef Pythia8_SigmaSUSY_H
#define Pythia8_SigmaSUSY_H

#include "SigmaProcessGPU.h"

#include "Pythia8/EventGPU.h"
// #include "Pythia8/PhaseSpace.h"
#include "Pythia8/PythiaComplexGPU.h"
#include "Pythia8/SigmaProcessGPU.h"
#include "Pythia8/SusyCouplingsGPU.h"

namespace Pythia8 {

// q qbar -> neutralino_i neutralino_j.
struct Sigma2qqbar2chi0chi0Data : public SigmaProcessData
{
  int_star     id3chi, id4chi, codeSave;
  string_star  nameSave;
  double_star  sigma0, ui, uj, ti, tj, openFracPair;
  complex*     propZ;
};

// q qbar -> neutralino_i chargino_j.
struct Sigma2qqbar2charchi0Data : public Sigma2qqbar2chi0chi0Data
{
  complex propW;
};

// q qbar -> chargino+_i chargino-_j.
struct Sigma2qqbar2charcharData : public Sigma2qqbar2chi0chi0Data
{
};

// q g -> neutralino_i squark_j (and cc)
struct Sigma2qg2chi0squarkData : public SigmaProcessData
{
  int_star     id3chi, id4sq, codeSave;
  string_star  nameSave;
  double_star  sigma0, ui, uj, ti, tj, openFracPair;
};

// q g -> chargino_i squark_j (incl cc)
struct Sigma2qg2charsquarkData : public Sigma2qg2chi0squarkData
{
  int_star id3Sav, id4Sav;
};

// q q' -> ~q_i ~q_j
struct Sigma2qq2squarksquarkData : public SigmaProcessData
{
  int_star id3Sav, id4Sav, codeSave, iGen3, iGen4, nNeut;
  string_star   nameSave;
  bool_star     isUD, onlyQCD;
  double_star   m2Glu;
  double** m2Neut;
  double** m2Char;
  double_star sigmaChar, sigmaNeut, sigmaGlu;
  double_star sigmaCharNeut, sigmaCharGlu, sigmaNeutGlu;
  double_star openFracPair;
  double_star tGlu, uGlu;
  double** tNeut, uNeut, tChar, uChar;
  double_star sumCt, sumCu, sumNt, sumNu, sumGt, sumGu, sumInterference;
};

// q qbar' -> ~q_i ~q*_j
struct Sigma2qqbar2squarkantisquarkData : public SigmaProcessData
{
  int_star     id3Sav, id4Sav, codeSave, iGen3, iGen4, nNeut;
  string_star  nameSave;
  bool_star    isUD, isCC, onlyQCD;
  double_star m2Glu;
  double** m2Neut;
  double_star xW;
  double_star openFracPair;
  double_star sigmaEW, sigmaGlu, sigmaEWG;
  double_star tGlu, uGlu;
  double** tNeut;
  double** uNeut;
  complex* propZW;
  double_star sumColS, sumColT, sumInterference;
};

// g g -> ~q ~q*
struct Sigma2gg2squarkantisquarkData : public SigmaProcessData
{

  int_star     id3Sav, id4Sav, codeSave;
  string_star  nameSave;
  double_star sigma, m2Sq, openFracPair;
};

// q g -> ~q ~g
struct Sigma2qg2squarkgluinoData : public SigmaProcessData
{
  int_star     id3Sav, codeSave;
  string_star  nameSave;
  double_star sigmaA, sigmaB, comFacHat, m2Glu, m2Sq, openFracPair;
};

// g g -> gluino gluino.
struct Sigma2gg2gluinogluinoData : public SigmaProcessData
{
  double_star sigTS, sigUS, sigTU, sigSum, sigma, openFracPair;
};

// q qbar -> gluino gluino.
struct Sigma2qqbar2gluinogluinoData : public SigmaProcessData
{
  double_star openFracPair, s34Avg, sigS, tHG, uHG, tHG2, uHG2;
};

// q q -> ~q
struct Sigma1qq2antisquarkData : public SigmaProcessData
{
  double_star mRes, GammaRes, m2Res, sigBW, widthOut;
  int_star    codeSave, idRes;
  string_star nameSave;
};

// q qbar -> neutralino_i gluino.
struct Sigma2qqbar2chi0gluinoData : public SigmaProcessData
{
  int_star     id3chi, id4chi, codeSave;
  string_star  nameSave;
  double_star  sigma0, ui, uj, ti, tj, openFracPair;
};

// q qbar -> neutralino_i chargino_j.
struct Sigma2qqbar2chargluinoData : public Sigma2qqbar2chi0gluinoData
{
  complex* propW;
};

// q qbar' -> ~q_i ~q*_j
struct Sigma2qqbar2sleptonantisleptonData : public Sigma2qqbar2squarkantisquarkData
{
  int_star     id3Sav, id4Sav, codeSave, iGen3, iGen4, nNeut;
  string_star  nameSave;
  bool_star    isUD;
  double**    vec_m2Neut;
  double_star xW;
  double_star openFracPair;
  double_star sigmaEW;
  double** vec_tNeut;
  double** vec_uNeut;
  complex* propZW;
  double_star sumColS, sumColT, sumInterference;
};

} // end namespace Pythia8

#endif // Pythia8_SigmaSUSY_H
