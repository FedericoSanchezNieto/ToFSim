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
  
  std::vector<double> wfdata = gwf.GeneratefromdEdx(Ll,energy,time0); // Returns the signal in each of the samples.

  double time1 =0; 
  // This will be the other end of the bar 
  std::vector<double> wfdata2 = gwf.GeneratefromdEdx(L-Ll,energy,time1); // Returns the signal in each of the samples. 

  std::cout << " Time 0 " << time0 << "  " << time1 << std::endl;
  
  for(int i =0 ; i < wfdata.size(); i++ ) {
    
    std::cout << "Sample["<<i<<"]= \t" << wfdata[i] << "   \t" << wfdata2[i] << std::endl; 
    
  }


  
  
  return 0;
}

