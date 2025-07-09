// -------------------------------------------------------------------------------------------
//    Photon production, transport, detection electronics simulation including sampling 
//    
//        Author: F.Sanchez (Universite de Geneve) 
//        V1 :  July 2025 
// --------------------------------------------------------------------------------------------

#ifndef GenWF_h
#define GenWF_h 

#include <TROOT.h>
#include <TRandom.h>
#include <TH1D.h>
#include <TF1.h> 

#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::istringstream

#include <fstream>
#include <map>

Double_t WFResponsefunction(Double_t *x, Double_t *par)
{
   Float_t xx =x[0];
   Double_t f = par[0]*TMath::Exp(-(xx-par[1])*(xx-par[1]) / ( 2 * TMath::Power( par[2] * ( xx - par[1] ) + par[3] , 2 ) ) ) ;
   return f;
}

class  GenWF{
 public :

  GenWF();
  GenWF(std::string jsonfile );
  virtual ~GenWF() {delete r; }
  void PrintSettings(void);

  std::vector<double> GeneratefromdEdx(double L,double energy,double &Cell0Time); 
  std::vector<double> Generate(double L,double nphotons,double &Cell0Time);

  double GetVelocity(void){ return c;}
  double GetScintFallTime(void) { return FallTime;}
  double GetBarLength(void) { return BarLength; } 

  // Counter getters
  double NPhotonsProduced(void) { return NPhotons; }
  double NPhotonsInAcceptance(void) { return NPhevt; }
  double NPhotonsSimultated(void) { return NPhcnt; }
  

 private :

  // Geometry 
  double BarLength; // in m
  double ReflectedFraction; // % are reflected in the other end.

  // Scintillator 
  
  double n; //index of refraction;
  double c; // velocity m/ns   
  double RaiseTime; // ns
  double FallTime; // ns 
  double Attenuation; // m
  double AttenuationShort; // m
  double FracShort; 
  double photonsperdedx;

  // Time smearing from Mie Scattering, wrong internal reflection, etc... 
  double TimeSmear; // probability, 0 if disabled  

  // MPPC's 
  double quantumefficiency;
  double fillfactor;
  double photocoverage; // Surface covered by the SiPM
  double PhotonEff;
  double Amplitude; 

  // Electronics
  
  int TotSamples;
  int PreSamples;
  int PostSamples; 
  double SamplingTime; // ns
  double SamplingCell0Time; // ns
  double threshold;
  double noise;

  // Counters
  double NPhotons;
  double NPhevt;
  double NPhcnt; 

  // WF
  double WF0;
  double WF1;
  double WF2;
  double WF3;

  TRandom *r;

  TF1 *fsc;
  TF1 *wfrf;
  TH1F *twf;

};


GenWF::GenWF(){

  // Some default configuration. 
  
  BarLength = 2.2; // in m 
  n = 1.58 ; //index of refraction;
  c= 300000000./n*1e-9; // m/ns
  RaiseTime = 0.9; // ns
  FallTime = 2.1; // ns 
  Attenuation = 3.8; // m
  AttenuationShort = 3.8; // m
  FracShort = 0.21; 
  
  SamplingTime = 0.3125; // ns 
  ReflectedFraction = 0.20; // % are reflected in the other end.

  quantumefficiency = 0.40;
  fillfactor = 0.74;
  photocoverage = (6.*6.)*8./10./120.;

  photonsperdedx = 10000./2.2;

  PhotonEff = quantumefficiency*fillfactor*photocoverage;

  noise = 0.2; 

  threshold = 1.;

  TotSamples = 60;
  PreSamples = 15;

  PostSamples = TotSamples-PreSamples; 
  
  PrintSettings(); 
  
  r = new TRandom();

  fsc  = new TF1("fsc","(1.-TMath::Exp(-x/[0]))*TMath::Exp(-x/[1])",0.,30.);
  fsc->SetParameters(RaiseTime,FallTime); 

  wfrf = new TF1("wfrf",WFResponsefunction,0.,1000.,4);
  wfrf->SetParameters(0.01,100.,0.29,7.0);
 
}


std::string getSubstringBetween(const std::string& str, char startChar, char endChar) {
    size_t startPos = str.find(startChar);
    if (startPos == std::string::npos) {
        return "";  // Start character not found
    }

    size_t endPos = str.find(endChar, startPos + 1);
    if (endPos == std::string::npos) {
        return "";  // End character not found
    }

    return str.substr(startPos + 1, endPos - startPos - 1);
}

GenWF::GenWF(std::string jsonfile ){

  std::map<std::string,double> config;

  std::ifstream file(jsonfile); // Open the file "data.txt"

  if (!file.is_open()) {
    std::cerr << "Error opening the file." << std::endl;
    return;
  }
  
  std::string line;
 
  while (std::getline(file, line)) {
    // Process each line as needed
    std::istringstream iss(line);
    std::string str;
    double number;
    iss >> str >> number;

    if( str == "{" || str == "}" || str == "" ) continue; 
   
    config[ getSubstringBetween(str,'"','"')] = number; 
  }
 
  BarLength = config["BarLength"]; // in m

  n = config["IndexRefraction"]; //index of refraction;
  c= 300000000/n*1e-9; // m/ns
  RaiseTime = config["RaiseTime"]; // ns
  FallTime = config["FallTime"]; // ns 
  Attenuation = config["Attenuation"]; // m
  AttenuationShort = config["AttenuationShort"]; // m
  FracShort = config["FracShort"]; // m
  SamplingTime = config["SamplingTime"]; // ns

  SamplingCell0Time = config["SamplingCell0Time"]; // ns
    
  ReflectedFraction = config["ReflectedFraction"]; // % are reflected in the other end.

  TimeSmear = config["TimeSmear"]; // Time Smear 
  
  quantumefficiency = config["quantumefficiency"];
  fillfactor = config["fillfactor"];
  photocoverage = config["photocoverage"];

  photonsperdedx = config["photonsperdedx"];

  PhotonEff = quantumefficiency*fillfactor*photocoverage;

  noise = config["noise"]; 

  threshold = config["threshold"];

  TotSamples = config["TotSamples"];
  PreSamples = config["PreSamples"];

  Amplitude = config["Amplitude"];

  PostSamples = TotSamples-PreSamples; 
 
  r = new TRandom();

  fsc  = new TF1("fsc","(1.-TMath::Exp(-x/[0]))*TMath::Exp(-x/[1])",0.,30.);
  fsc->SetParameters(RaiseTime,FallTime); 

  wfrf = new TF1("wfrf",WFResponsefunction,0.,1000.,4);
  WF0 = config["WF0"];
  WF1 = config["WF1"];
  WF2 = config["WF2"];
  WF3 = config["WF3"];
  //  std::cout << WF0 << std::endl;
  wfrf->SetParameters(WF0,WF1,WF2,WF3);

  //  TObject *a = NULL; 
  //  if( (a = gROOT->FindObject("twf")) ) delete a;
  twf = new TH1F("twf","  ",301,-0.5,300.5);

  file.close();

  PrintSettings(); 
  
}



void GenWF::PrintSettings(void) {

  std::cout << std::endl<< std::endl;
  std::cout << " MC photon generation settings " << std::endl<<std::endl;
  std::cout << " ----------------------------- " << std::endl;
  std::cout << " -------- Geometry  ---------- " << std::endl;
  std::cout << " ----------------------------- " << std::endl;
  std::cout << " Bar Length " << BarLength << " m " << std::endl;
  std::cout << " Bar photo coverage " << photocoverage << std::endl;
  std::cout << " ----------------------------- " << std::endl;
  std::cout << " --------   MPPC    ---------- " << std::endl;
  std::cout << " ----------------------------- " << std::endl;
  std::cout << " Quantum efficiency " << quantumefficiency << std::endl;
  std::cout << " Pixel fill factor " << fillfactor << std::endl;
  std::cout << " Total photon detection efficiency " << PhotonEff << std::endl;
  std::cout << " Amplitude " << Amplitude << std::endl;
  std::cout << " ----------------------------- " << std::endl;
  std::cout << " ------- Scintillator -------- " << std::endl;
  std::cout << " ----------------------------- " << std::endl;
  std::cout << " Index refraction " << n << std::endl; 
  std::cout << " Raise Time " << RaiseTime << " ns" << std::endl; 
  std::cout << " Fall Time " << FallTime << " ns" << std::endl;
  std::cout << " Photon speed " << c <<  " m/ns" << std::endl;
  std::cout << " Reflection fraction " << ReflectedFraction << std::endl;
  std::cout << " Attenuation " << Attenuation << " m" << std::endl;
  std::cout << " AttenuationShort " << AttenuationShort << " m" << std::endl;
  std::cout << " FracShort " << FracShort << " m" << std::endl;
  std::cout << " Photonsperdedx " <<  photonsperdedx << " photons/MeV " << std::endl;
  std::cout << " TimeSmear " << TimeSmear << " ns " << std::endl;
  std::cout << " ----------------------------- " << std::endl;
  std::cout << " ------ Electronics  --------- " << std::endl;
  std::cout << " ----------------------------- " << std::endl;
  std::cout << " Sampling time " << SamplingTime << std::endl;
  std::cout << " Sampling Cell0 time " << SamplingCell0Time << std::endl;
  std::cout << " Electronics noise " << noise << std::endl;
  std::cout << " Trigger threshold " << threshold << std::endl;
  std::cout << " Waveform total samples " << TotSamples << std::endl;
  std::cout << " Waveform pre-samples " << PreSamples << std::endl;
  std::cout << " ----------------------------- " << std::endl;
  std::cout << " ------ WF parameters  --------- " << std::endl;
  std::cout << " ----------------------------- " << std::endl;
  printf("%3.2f*TMath::Exp(-(t-%3.2f)^2)/(2*(%3.2f*(t-%3.2f)+%3.2f)^2)\n",WF0,WF1,WF2,WF1,WF3) ;
  std::cout << " ----------------------------- " << std::endl;
}

std::vector<double> GenWF::GeneratefromdEdx(double L, double E,double &Cell0Time){

  double nphotons = photonsperdedx*E;

  //  std::cout << ReflectedFraction << "   " << nphotons << " --- > " ; 

  
  nphotons = nphotons/2.*(1.+ReflectedFraction); // Number of photons in one direction is half plus the reflected. 

  //  std::cout << nphotons << std::endl;
  
  return Generate(L,nphotons,Cell0Time); 
}

std::vector<double> GenWF::Generate(double L, double nphotons,double &Cell0Time) {

  NPhotons = nphotons; 
  
  // Photon detection efficiency 
  nphotons *= PhotonEff;

  twf->Reset();
  
  double mincos =  TMath::Sqrt(1.-1/n*1/n); // there is no photons below this value. Used for MC generation. 

  // mincos = 0.;
  
  double nphotonseff = nphotons*(1.-mincos);   // Effective number of photons. 

  //  std::cout << (1.-mincos) << std::endl; 
  
  int nphevt = r->Poisson(nphotonseff);

  //  std::cout << nphevt << "  " << nphotons << std::endl; 

  int nphcnt = 0;

  double timesamplingphase = r->Uniform(0.,SamplingCell0Time); // Constant for the event

  double tphotonmin =  L/c+timesamplingphase; // This is the first time arrival of a photon in the event. 

  for(int i = 0;i < nphevt; i++ ) {
    
    double cos = r->Uniform(mincos,1.); 
    double z = cos;

    // Simulate the returning photons.
    double d = 0;
    
    if( r->Uniform(0.,1.) > ReflectedFraction/(1.+ReflectedFraction) ) { // Direct photons 
      d = L/z;
    }
    else {  //Reflected photons, 
      d = (2.*BarLength-L)/z;
    }
      
    // Include the attenuation at the beginning to reduce computation.
    double u = r->Uniform(0.,1.);
    double v = r->Uniform(0.,1.);

    // Select between short and long Attenuation. 
    double Att = Attenuation; 
    if( v < FracShort ) Att = AttenuationShort;
    
    if( u > TMath::Exp(-d/Att) ) continue; // Attenuation 

    double sin = TMath::Sqrt(1.-cos*cos); 
    
    double phi = r->Uniform(-3.141592,3.141592);
    double x = sin*TMath::Cos(phi);
    double y = sin*TMath::Sin(phi);
    
    double siny = TMath::Sqrt(1.-y*y); //refraction angle in the horizontal boundaries
    double sinx = TMath::Sqrt(1.-x*x); //refraction angle in the vertical boundaries 
    
    bool directphoton = false;
    
    double sinymax = L/TMath::Sqrt(L*L+0.005*0.005);
    double sinxmax = L/TMath::Sqrt(L*L+0.06*0.06); 
    
    if( siny > sinymax  || sinx > sinxmax  ) directphoton = true; // Check for direct photons entering.
    
    bool internalreflection = (siny  > 1/n && sinx > 1/n ) ;
    
    if(  internalreflection || directphoton  ) {  // internal reflection criteria for all bars
      if( sin < 1/n ) { // photons leave the bar at the end. 
	nphcnt++; 
  
	//Randomise the scintillator de-excitation 
	double time = d/c+fsc->GetRandom();

	// Add gaussian noise from Mie scattering and bad reflections. This really changes the arrival time of a photon. 
	time +=  r->Gaus(0,TimeSmear*(d/L));

	//And the sampling phase
	time +=  timesamplingphase;

	// And the electronics shaper 
	for(int j = 0; j < twf->GetNbinsX();j++ ) {
	  // The WF histogram contains the time with respect the particle t0 
	  double tj = time/SamplingTime+j;    
	  double ftj = wfrf->Eval((double)j-0.5); // wfrf->Integral((double)j-0.5,(double)j+0.5); // WF Shape
	  ftj += wfrf->Eval((double)j) ;
	  ftj += wfrf->Eval((double)j+0.5) ;
	  ftj /= 3.;
	  twf->Fill(tj,ftj);
	}
      }
    }
    
  } // end loop one event

  NPhevt = nphevt;
  NPhcnt = nphcnt; 
  
  std::vector<double> wf;

  int ith = -10;

  // Add noise and get the threshold 
  
  for(int i = 0; i < twf->GetNbinsX();i++) {

    double val = Amplitude*twf->GetBinContent(i+1)+r->Gaus(0.,noise);

    if( val > threshold && ith == -10 ) {
      ith = i;
    }
  
    wf.push_back(val); 
  }

  //
  // Find the fine time of the thrshold crossing.
  // Simple linear interpolation might be sufficient given the
  // precision in the SamplingCell0Time.
  //
  double dw = (wf[ith]-wf[ith-1]);
  double dt = SamplingTime; 
  
  double ftime = ((threshold-wf[ith-1])/dw)*dt;
  
  ftime = ((int)(ftime/SamplingCell0Time))*SamplingCell0Time; 


  // Copy the valid section of the WF 
  std::vector<double> wfpub(&wf[ith-PreSamples],&wf[ith+PostSamples]); 


  // Recompute the Cell0Time to the beginning of the WF 
  Cell0Time = (ith-PreSamples)*SamplingTime+ftime;
  
  return wfpub; 
}


#endif 
