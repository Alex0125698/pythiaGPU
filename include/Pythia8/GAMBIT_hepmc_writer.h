#ifndef Pythia8_GAMBIT_hepmc_writer_H
#define Pythia8_GAMBIT_hepmc_writer_H

#include "HepMC3/GenEvent.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "HepMC3/Print.h"

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"


namespace Pythia8
{

  class GAMBIT_hepmc_writer
  {

  private:

    std::string filename2;
    std::string filename3;

    bool HepMC2_ON = false;
    bool HepMC3_ON = false;
 
    HepMC3::Pythia8ToHepMC3 *pythiaToHepMC;
    HepMC3::WriterAsciiHepMC2 *file2;
    HepMC3::WriterAscii *file3;
 

  public:

    // Destructor
    ~GAMBIT_hepmc_writer()
    {
      if(HepMC2_ON)
      {
        file2->close();
        delete file2;
      }

      if(HepMC3_ON)
      {
        file3->close();
        delete file3;
      }

      delete pythiaToHepMC;
    }

    void init(std::string filename_in, bool HepMC2, bool HepMC3)
    {
      if (HepMC2)
      {
        filename2 = filename_in+"2";
        file2 = new HepMC3::WriterAsciiHepMC2(filename2.c_str());
        HepMC2_ON = true;
      }

      if (HepMC3)
      {
        filename3 = filename_in+"3";
        file3 = new HepMC3::WriterAscii(filename3.c_str());
        HepMC3_ON = true;
      }

      pythiaToHepMC = new HepMC3::Pythia8ToHepMC3();
    }


    // Write current event to file in HepMC3 format
    void write_event_HepMC3(Pythia* pythia)
    {
      // Construct new empty HepMC event and fill it.
      HepMC3::GenEvent hepmc( HepMC3::Units::GEV, HepMC3::Units::MM );

      pythiaToHepMC->fill_next_event(pythia->event, &hepmc, -1, &(pythia->info));

      // Write the HepMC event to file. Done with it.
      file3->write_event(hepmc);
    }
  
    // Write current event to file in HepMC2 format
    void write_event_HepMC2(Pythia* pythia)
    {
      // Construct new empty HepMC event and fill it.
      HepMC3::GenEvent hepmc( HepMC3::Units::GEV, HepMC3::Units::MM );

      pythiaToHepMC->fill_next_event(pythia->event, &hepmc, -1, &(pythia->info));

      // Write the HepMC event to file. Done with it.
      file2->write_event(hepmc);
    }

    // Returns current event in HepMC format
    void convert_to_HepMC_event(Pythia *pythia, HepMC3::GenEvent &event)
    {
      pythiaToHepMC->fill_next_event(pythia->event, &event, -1, &(pythia->info));
    }
  
};

}

#endif // Pythia8_GAMBIT_hepmc_writer_H