// -------------------------------------------------------------------------------------------
//    Example of usage of the GenWF class for photons read at both ends of the bar. 
//    
//        Author: F.Sanchez (Universite de Geneve) 
//        V1 :  July 2025 
// --------------------------------------------------------------------------------------------

#include "GenWF.h"

int main(int argc, char **argv){

  GenWF gwf("config.json");

  double L = gwf.GetBarLength();

  double Ll = 0.5; // Set the energy in the center of the bar. 
  

  TRandom r; 
  
  double energy = r.Landau(2.0,0.22);

  double time0 = 0; // This is the time of the first sample in the sample array 

  // Generate photons at a distance Ll from the sensor, for an "energy", returning the time of the first cell and the content of the cell.  

  wf.GeneratefromdEdx(Ll,energy);


  std::vector<double> wfdata_barend0 = gwf.GetWaveForm(0); // Returns the signal in each of the samples.
  double time0 = gwf.GetCell0Time(0);
  std::vector<double> wfdata_barend1 = gwf.GetWaveForm(1); // Returns the signal in each of the samples.
  double time1 = gwf.GetCell0Time(1);

  std::cout << " Time 0 " << time0 << "  " << time1 << std::endl;
  
  for(int i =0 ; i < wfdata.size(); i++ ) {
    
    std::cout << "Sample["<<i<<"]= \t" << wfdata_barend0[i] << "   \t" << wfdata_barend1[i] << std::endl; 
    
  }


  
  
  return 0;
}

