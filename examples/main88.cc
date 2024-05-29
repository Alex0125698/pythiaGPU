// main88.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This program is written by Stefan Prestel.
// It illustrates how to do UNLOPS merging, see the Matrix Element
// Merging page in the online manual. An example command is
//     ./main88 main88.cmnd w_production hepmcout88.dat
// where main88.cmnd supplies the commands, w_production provides the
// input LHE events, and hepmcout88.dat is the output file. This
// example requires HepMC.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include <unistd.h>

using namespace Pythia8;

//==========================================================================

// Example main programm to illustrate UNLOPS merging

int main( int argc, char* argv[] ){

  // Check that correct number of command-line arguments
  if (argc != 4) {
    cerr << " Unexpected number of command-line arguments ("<<argc<<"). \n"
         << " You are expected to provide the arguments" << endl
         << " 1. Input file for settings" << endl
         << " 2. Name of the input LHE file (with path), up to the '_tree'"
         << " or '_powheg' identifiers" << endl
         << " 3. Output hepmc file name" << endl
         << " Program stopped. " << endl;
    return 1;
  }

  Pythia pythia;

  // Input parameters:
  //  1. Input file for settings
  //  2. Path to input LHE file
  //  3. Output histogram path
  pythia.readFile(argv[1]);

  // Interface for conversion from Pythia8::Event to HepMC one.
  HepMC::Pythia8ToHepMC ToHepMC;
  // Specify file where HepMC events will be stored.
  HepMC::IO_GenEvent ascii_io(argv[3], std::ios::out);
  // Switch off warnings for parton-level events.
  ToHepMC.set_print_inconsistency(false);
  ToHepMC.set_free_parton_warnings(false);
  // Do not store cross section information, as this will be done manually.
  ToHepMC.set_store_pdf(false);
  ToHepMC.set_store_proc(false);
  ToHepMC.set_store_xsec(false);

  // Path to input events, with name up to the "_tree", "_powheg" identifier
  // included.
  string iPath = string(argv[2]);

  // Number of events
  int nEvent   = pythia.get(Mode::Main_numberOfEvents);
  // Maximal number of additional LO jets.
  int nMaxLO   = pythia.get(Mode::Merging_nJetMax);
  // maximal number of additional NLO jets.
  int nMaxNLO  = pythia.get(Mode::Merging_nJetMaxNLO);

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

  // Switch off all showering and MPI when extimating the cross section after
  // the merging scale cut.
  bool fsr = pythia.get(Flag::PartonLevel_FSR);
  bool isr = pythia.get(Flag::PartonLevel_ISR);
  bool mpi = pythia.get(Flag::PartonLevel_MPI);
  bool had = pythia.get(Flag::HadronLevel_all);
  pythia.settings.set(Flag::PartonLevel_FSR, false);
  pythia.settings.set(Flag::PartonLevel_ISR, false);
  pythia.settings.set(Flag::HadronLevel_all, false);
  pythia.settings.set(Flag::PartonLevel_MPI, false);

  // Switch on cross section estimation procedure.
  pythia.settings.set(Flag::Merging_doXSectionEstimate, true);
  pythia.settings.set(Flag::Merging_doUNLOPSTree, true);

  int njetcounterLO  = nMaxLO;
  string iPathTree   = iPath + "_tree";

  // Save estimates in vectors.
  vector<double> xsecLO;
  vector<double> nSelectedLO;
  vector<double> nAcceptLO;
  vector<int> strategyLO;

  cout << endl << endl << endl;
  cout << "Start estimating unlops tree level cross section" << endl;

  while(njetcounterLO >= 0){

    // From njetcounter, choose LHE file
    stringstream in;
    in   << "_" << njetcounterLO << ".lhe";
#ifdef GZIPSUPPORT
    if(access( (iPathTree+in.str()+".gz").c_str(), F_OK) != -1) in << ".gz";
#endif
    string LHEfile = iPathTree + in.str();
    pythia.settings.set(Mode::Merging_nRequested, njetcounterLO);
    pythia.settings.set(Mode::Beams_frameType, 4);
    pythia.settings.set(Word::Beams_LHEF, LHEfile);
    pythia.init();

    // Start generation loop
    for( int iEvent=0; iEvent<nEvent; ++iEvent ){
      // Generate next event
      if( !pythia.next() ) {
        if( pythia.info.atEndOfFile() ){
          break;
        }
        else continue;
      }
    } // end loop over events to generate

    // print cross section, errors
    pythia.stat();

    xsecLO.push_back(pythia.info.sigmaGen());
    nSelectedLO.push_back(pythia.info.nSelected());
    nAcceptLO.push_back(pythia.info.nAccepted());
    strategyLO.push_back(pythia.info.lhaStrategy());

    // Restart with ME of a reduced the number of jets
    if( njetcounterLO > 0 )
      njetcounterLO--;
    else
      break;

  }

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

  cout << endl << endl << endl;
  cout << "Start estimating unlops virtual corrections cross section" << endl;

  pythia.settings.set(Flag::Merging_doUNLOPSTree, false);
  pythia.settings.set(Flag::Merging_doUNLOPSLoop, true);

  int njetcounterNLO  = nMaxNLO;
  string iPathLoop    = iPath + "_powheg";

  // Save estimates in vectors.
  vector<double> xsecNLO;
  vector<double> nSelectedNLO;
  vector<double> nAcceptNLO;
  vector<int> strategyNLO;

  while(njetcounterNLO >= 0){

    // From njetcounter, choose LHE file
    stringstream in;
    in   << "_" << njetcounterNLO << ".lhe";
#ifdef GZIPSUPPORT
    if(access( (iPathLoop+in.str()+".gz").c_str(), F_OK) != -1) in << ".gz";
#endif
    string LHEfile = iPathLoop + in.str();
    pythia.settings.set(Mode::Merging_nRequested, njetcounterNLO);
    pythia.settings.set(Mode::Beams_frameType, 4);
    pythia.settings.set(Word::Beams_LHEF, LHEfile);
    pythia.init();

    // Start generation loop
    for( int iEvent=0; iEvent<nEvent; ++iEvent ){
      // Generate next event
      if( !pythia.next() ) {
        if( pythia.info.atEndOfFile() ){
          break;
        }
        else continue;
      }
    } // end loop over events to generate

    // print cross section, errors
    pythia.stat();

    xsecNLO.push_back(pythia.info.sigmaGen());
    nSelectedNLO.push_back(pythia.info.nSelected());
    nAcceptNLO.push_back(pythia.info.nAccepted());
    strategyNLO.push_back(pythia.info.lhaStrategy());

    // Restart with ME of a reduced the number of jets
    if( njetcounterNLO > 0 )
      njetcounterNLO--;
    else
      break;

  }

  int sizeLO   = int(xsecLO.size());
  int sizeNLO  = int(xsecNLO.size());

  // Switch off cross section estimation.
  pythia.settings.set(Flag::Merging_doXSectionEstimate, false);

  // Switch showering and multiple interaction back on.
  pythia.settings.set(Flag::PartonLevel_FSR, fsr);
  pythia.settings.set(Flag::PartonLevel_ISR, isr);
  pythia.settings.set(Flag::HadronLevel_all, had);
  pythia.settings.set(Flag::PartonLevel_MPI, mpi);

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

  // Declare sample cross section for output.
  double sigmaTemp  = 0.;
  vector<double> sampleXStree;
  vector<double> sampleXSvirt;
  vector<double> sampleXSsubtTree;
  vector<double> sampleXSsubtVirt;
  // Cross section an error.
  double sigmaTotal  = 0.;
  double errorTotal  = 0.;

  // Switch on tree-level processing.
  pythia.settings.set(Flag::Merging_doUNLOPSTree, true);
  pythia.settings.set(Flag::Merging_doUNLOPSLoop, false);
  pythia.settings.set(Flag::Merging_doUNLOPSSubt, false);
  pythia.settings.set(Flag::Merging_doUNLOPSSubtNLO, false);
  pythia.settings.set(Mode::Merging_nRecluster, 0);

  // Start looping through input event files.
  njetcounterLO = nMaxLO;
  iPathTree     = iPath + "_tree";

  while(njetcounterLO >= 0){

    // From njetcounter, choose LHE file
    stringstream in;
    in   << "_" << njetcounterLO << ".lhe";
#ifdef GZIPSUPPORT
    if(access( (iPathTree+in.str()+".gz").c_str(), F_OK) != -1) in << ".gz";
#endif
    string LHEfile = iPathTree + in.str();

    cout << endl << endl << endl
         << "Start tree level treatment for " << njetcounterLO << " jets"
         << endl;

    // UNLOPS does not contain a zero-jet tree-level sample.
    if ( njetcounterLO == 0 ) break;
    pythia.settings.set(Mode::Merging_nRequested, njetcounterLO);
    pythia.settings.set(Mode::Beams_frameType, 4);
    pythia.settings.set(Word::Beams_LHEF, LHEfile);
    pythia.init();

    // Remember position in vector of cross section estimates.
    int iNow = sizeLO-1-njetcounterLO;

    // Start generation loop
    for( int iEvent=0; iEvent<nEvent; ++iEvent ){

      // Generate next event
      if( !pythia.next() ) {
        if( pythia.info.atEndOfFile() ) break;
        else continue;
      }

      // Get event weight(s).
      double weightNLO  = pythia.info.mergingWeightNLO();
      double evtweight  = pythia.info.weight();
      weightNLO        *= evtweight;
      // Do not print zero-weight events.
      if ( weightNLO == 0. ) continue;

      // Construct new empty HepMC event.
      HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
      // Get correct cross section from previous estimate.
      double normhepmc = xsecLO[iNow] / nAcceptLO[iNow];
      // Set hepmc event weight.
      hepmcevt->weights().push_back(weightNLO*normhepmc);
      // Fill HepMC event.
      ToHepMC.fill_next_event( pythia, hepmcevt );
      // Add the weight of the current event to the cross section.
      sigmaTotal += weightNLO*normhepmc;
      sigmaTemp  += weightNLO*normhepmc;
      errorTotal += pow2(weightNLO*normhepmc);
      // Report cross section to hepmc.
      HepMC::GenCrossSection xsec;
      xsec.set_cross_section( sigmaTotal*1e9, pythia.info.sigmaErr()*1e9 );
      hepmcevt->set_cross_section( xsec );
      // Write the HepMC event to file. Done with it.
      ascii_io << hepmcevt;
      delete hepmcevt;
    } // end loop over events to generate

    // print cross section, errors
    pythia.stat();

    // Restart with ME of a reduced the number of jets
    if( njetcounterLO > 0 )
      njetcounterLO--;
    else
      break;

    // Save sample cross section for output.
    sampleXStree.push_back(sigmaTemp);
    sigmaTemp = 0.;

  }

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

  cout << endl << endl << endl;
  cout << "Start unlops virtual corrections part" << endl;

  // Switch on loop-level processing.
  pythia.settings.set(Flag::Merging_doUNLOPSTree, false);
  pythia.settings.set(Flag::Merging_doUNLOPSLoop, true);
  pythia.settings.set(Flag::Merging_doUNLOPSSubt, false);
  pythia.settings.set(Flag::Merging_doUNLOPSSubtNLO, false);
  pythia.settings.set(Mode::Merging_nRecluster, 0);

  njetcounterNLO = nMaxNLO;
  iPathLoop= iPath + "_powheg";

  while(njetcounterNLO >= 0){

    // From njetcounter, choose LHE file
    stringstream in;
    in   << "_" << njetcounterNLO << ".lhe";
#ifdef GZIPSUPPORT
    if(access( (iPathLoop+in.str()+".gz").c_str(), F_OK) != -1) in << ".gz";
#endif
    string LHEfile = iPathLoop + in.str();

    cout << endl << endl << endl
         << "Start loop level treatment for " << njetcounterNLO << " jets"
         << endl;

    pythia.settings.set(Mode::Merging_nRequested, njetcounterNLO);
    pythia.settings.set(Mode::Beams_frameType, 4);
    pythia.settings.set(Word::Beams_LHEF, LHEfile);
    pythia.init();

    // Remember position in vector of cross section estimates.
    int iNow = sizeNLO-1-njetcounterNLO;

    // Start generation loop
    for( int iEvent=0; iEvent<nEvent; ++iEvent ){

      // Generate next event
      if( !pythia.next() ) {
        if( pythia.info.atEndOfFile() ) break;
        else continue;
      }

      // Get event weight(s).
      double weightNLO  = pythia.info.mergingWeightNLO();
      double evtweight  = pythia.info.weight();
      weightNLO        *= evtweight;
      // Do not print zero-weight events.
      if ( weightNLO == 0. ) continue;

      // Construct new empty HepMC event.
      HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
      // Get correct cross section from previous estimate.
      double normhepmc = xsecNLO[iNow] / nAcceptNLO[iNow];
      // powheg weighted events
      if( abs(strategyNLO[iNow]) == 4)
        normhepmc = 1. / (1e9*nSelectedNLO[iNow]);
      // Set hepmc event weight.
      hepmcevt->weights().push_back(weightNLO*normhepmc);
      // Fill HepMC event.
      ToHepMC.fill_next_event( pythia, hepmcevt );
      // Add the weight of the current event to the cross section.
      sigmaTotal += weightNLO*normhepmc;
      sigmaTemp  += weightNLO*normhepmc;
      errorTotal += pow2(weightNLO*normhepmc);
      // Report cross section to hepmc
      HepMC::GenCrossSection xsec;
      xsec.set_cross_section( sigmaTotal*1e9, pythia.info.sigmaErr()*1e9 );
      hepmcevt->set_cross_section( xsec );
      // Write the HepMC event to file. Done with it.
      ascii_io << hepmcevt;
      delete hepmcevt;

    } // end loop over events to generate

    // print cross section, errors
    pythia.stat();
    // Save sample cross section for output.
    sampleXSvirt.push_back(sigmaTemp);
    sigmaTemp = 0.;

    // Restart with ME of a reduced the number of jets
    if( njetcounterNLO > 0)
      njetcounterNLO--;
    else
      break;

  }

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

  cout << endl << endl << endl;
  cout << "Shower subtractive events" << endl;

  // Switch on processing of counter-events.
  pythia.settings.set(Flag::Merging_doUNLOPSTree, false);
  pythia.settings.set(Flag::Merging_doUNLOPSLoop, false);
  pythia.settings.set(Flag::Merging_doUNLOPSSubt, true);
  pythia.settings.set(Flag::Merging_doUNLOPSSubtNLO, false);
  pythia.settings.set(Mode::Merging_nRecluster, 1);

  int nMaxCT = nMaxLO;
  int njetcounterCT = nMaxCT;
  string iPathSubt= iPath + "_tree";

  while(njetcounterCT >= 1){

    // From njetcounter, choose LHE file
    stringstream in;
    in   << "_" << njetcounterCT << ".lhe";
#ifdef GZIPSUPPORT
    if(access( (iPathSubt+in.str()+".gz").c_str(), F_OK) != -1) in << ".gz";
#endif
    string LHEfile = iPathSubt + in.str();

    cout << endl << endl << endl
         << "Start subtractive treatment for " << njetcounterCT << " jets"
         << endl;

    pythia.settings.set(Mode::Merging_nRequested, njetcounterCT);
    pythia.settings.set(Mode::Beams_frameType, 4);
    pythia.settings.set(Word::Beams_LHEF, LHEfile);
    pythia.init();

    // Remember position in vector of cross section estimates.
    int iNow = sizeLO-1-njetcounterCT;

    // Start generation loop
    for( int iEvent=0; iEvent<nEvent; ++iEvent ){

      // Generate next event
      if( !pythia.next() ) {
        if( pythia.info.atEndOfFile() ) break;
        else continue;
      }

      // Get event weight(s).
      double weightNLO  = pythia.info.mergingWeightNLO();
      double evtweight  = pythia.info.weight();
      weightNLO        *= evtweight;
      // Do not print zero-weight events.
      if ( weightNLO == 0. ) continue;

      // Construct new empty HepMC event.
      HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
      // Get correct cross section from previous estimate.
      double normhepmc = -1*xsecLO[iNow] / nAcceptLO[iNow];
      // Set hepmc event weight.
      hepmcevt->weights().push_back(weightNLO*normhepmc);
      // Fill HepMC event.
      ToHepMC.fill_next_event( pythia, hepmcevt );
      // Add the weight of the current event to the cross section.
      sigmaTotal += weightNLO*normhepmc;
      sigmaTemp  += weightNLO*normhepmc;
      errorTotal += pow2(weightNLO*normhepmc);
      // Report cross section to hepmc.
      HepMC::GenCrossSection xsec;
      xsec.set_cross_section( sigmaTotal*1e9, pythia.info.sigmaErr()*1e9 );
      hepmcevt->set_cross_section( xsec );
      // Write the HepMC event to file. Done with it.
      ascii_io << hepmcevt;
      delete hepmcevt;

    } // end loop over events to generate

    // print cross section, errors
    pythia.stat();
    // Save sample cross section for output.
    sampleXSsubtTree.push_back(sigmaTemp);
    sigmaTemp = 0.;

    // Restart with ME of a reduced the number of jets
    if( njetcounterCT > 1 )
      njetcounterCT--;
    else
      break;

  }

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

  cout << endl << endl << endl;
  cout << "Shower subtractive events" << endl;

  pythia.settings.set(Flag::Merging_doUNLOPSTree, false);
  pythia.settings.set(Flag::Merging_doUNLOPSLoop, false);
  pythia.settings.set(Flag::Merging_doUNLOPSSubt, false);
  pythia.settings.set(Flag::Merging_doUNLOPSSubtNLO, true);
  pythia.settings.set(Mode::Merging_nRecluster, 1);

  nMaxCT = nMaxNLO;
  njetcounterCT = nMaxCT;
  iPathSubt= iPath + "_powheg";

  while(njetcounterCT >= 1){

    // From njetcounter, choose LHE file
    stringstream in;
    in   << "_" << njetcounterCT << ".lhe";
#ifdef GZIPSUPPORT
    if(access( (iPathSubt+in.str()+".gz").c_str(), F_OK) != -1) in << ".gz";
#endif
    string LHEfile = iPathSubt + in.str();

    cout << endl << endl << endl
         << "Start subtractive treatment for " << njetcounterCT << " nlo jets"
         << endl;

    pythia.settings.set(Mode::Merging_nRequested, njetcounterCT);
    pythia.settings.set(Mode::Beams_frameType, 4);
    pythia.settings.set(Word::Beams_LHEF, LHEfile);
    pythia.init();

    // Remember position in vector of cross section estimates.
    int iNow = sizeNLO-1-njetcounterCT;

    // Start generation loop
    for( int iEvent=0; iEvent<nEvent; ++iEvent ){

      // Generate next event
      if( !pythia.next() ) {
        if( pythia.info.atEndOfFile() ) break;
        else continue;
      }

      // Get event weight(s).
      double weightNLO  = pythia.info.mergingWeightNLO();
      double evtweight  = pythia.info.weight();
      weightNLO        *= evtweight;
      // Do not print zero-weight events.
      if ( weightNLO == 0. ) continue;

      // Construct new empty HepMC event.
      HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
      // Get correct cross section from previous estimate.
      double normhepmc = -1*xsecNLO[iNow] / nAcceptNLO[iNow];
      // powheg weighted events
      if( abs(strategyNLO[iNow]) == 4)
        normhepmc = -1. / (1e9*nSelectedNLO[iNow]);
      // Set hepmc event weight.
      hepmcevt->weights().push_back(weightNLO*normhepmc);
      // Fill HepMC event.
      ToHepMC.fill_next_event( pythia, hepmcevt );
      // Add the weight of the current event to the cross section.
      sigmaTotal += weightNLO*normhepmc;
      sigmaTemp  += weightNLO*normhepmc;
      errorTotal += pow2(weightNLO*normhepmc);
      // Report cross section to hepmc.
      HepMC::GenCrossSection xsec;
      xsec.set_cross_section( sigmaTotal*1e9, pythia.info.sigmaErr()*1e9 );
      hepmcevt->set_cross_section( xsec );
      // Write the HepMC event to file. Done with it.
      ascii_io << hepmcevt;
      delete hepmcevt;

    } // end loop over events to generate

    // print cross section, errors
    pythia.stat();
    // Save sample cross section for output.
    sampleXSsubtVirt.push_back(sigmaTemp);
    sigmaTemp = 0.;

    // Restart with ME of a reduced the number of jets
    if( njetcounterCT > 1 )
      njetcounterCT--;
    else
      break;

  }

  // Print cross section information.
  cout << endl << endl;
  cout << " *---------------------------------------------------*" << endl;
  cout << " |                                                   |" << endl;
  cout << " | Sample cross sections after UNLOPS merging        |" << endl;
  cout << " |                                                   |" << endl;
  cout << " | Leading order cross sections (mb):                |" << endl;
  for (int i = 0; i < int(sampleXStree.size()); ++i)
    cout << " |     " << sampleXStree.size()-1-i+1 << "-jet:  "
         << setw(17) << scientific << setprecision(6)
         << sampleXStree[i] << "                     |" << endl;
  cout << " |     (No 0-jet tree-level sample in UNLOPS)        |" << endl;
  cout << " |                                                   |" << endl;
  cout << " | NLO order cross sections (mb):                    |" << endl;
  for (int i = 0; i < int(sampleXSvirt.size()); ++i)
    cout << " |     " << sampleXSvirt.size()-1-i << "-jet:  "
         << setw(17) << scientific << setprecision(6)
         << sampleXSvirt[i] << "                     |" << endl;
  cout << " |                                                   |" << endl;
  cout << " | Leading-order subtractive cross sections (mb):    |" << endl;
  for (int i = 0; i < int(sampleXSsubtTree.size()); ++i)
    cout << " |     " << sampleXSsubtTree.size()-1-i+1 << "-jet:  "
         << setw(17) << scientific << setprecision(6)
         << sampleXSsubtTree[i] << "                     |" << endl;
  cout << " |                                                   |" << endl;
  if ( sampleXSsubtVirt.size() > 0) {
  cout << " | NLO subtractive cross sections (mb):              |" << endl;
  for (int i = 0; i < int(sampleXSsubtVirt.size()); ++i)
    cout << " |     " << sampleXSsubtVirt.size()-1-i+1 << "-jet:  "
         << setw(17) << scientific << setprecision(6)
         << sampleXSsubtVirt[i] << "                     |" << endl;
  cout << " |                                                   |" << endl;
  }
  cout << " |---------------------------------------------------|" << endl;
  cout << " |---------------------------------------------------|" << endl;
  cout << " | Inclusive cross sections:                         |" << endl;
  cout << " |                                                   |" << endl;
  cout << " | UNLOPS merged inclusive cross section:            |" << endl;
  cout << " |    " << setw(17) << scientific << setprecision(6)
       << sigmaTotal << "  +-  " << setw(17) << sqrt(errorTotal) << " mb "
       << "   |" << endl;
  cout << " |                                                   |" << endl;
  cout << " | NLO inclusive cross section:                      |" << endl;
  cout << " |    " << setw(17) << scientific << setprecision(6)
       << xsecNLO.back() << " mb                           |"  << endl;
  cout << " |                                                   |" << endl;
  cout << " | LO inclusive cross section:                       |" << endl;
  cout << " |    " << setw(17) << scientific << setprecision(6)
       << xsecLO.back() << " mb                           |" << endl;
  cout << " |                                                   |" << endl;
  cout << " *---------------------------------------------------*" << endl;
  cout << endl << endl;

  // Done
  return 0;

}
