// -------------------------------------------------------------------------------------------
//    Photon production, transport, detection electronics simulation including sampling 
//    
//        Author: F.Sanchez (Universite de Geneve) 
//          V1 :  July 2025
//          V2 :  21 August 2025
// --------------------------------------------------------------------------------------------


#ifndef GenWF_h
#define GenWF_h 

#include <TROOT.h>
#include <TRandom.h>
#include <TH1D.h>
#include <TF1.h>
#include <TSpline.h> 

#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::istringstream

#include <fstream>
#include <map>

//#define  __WFfunc__

Double_t WFResponsefunction(Double_t *x, Double_t *par)
{
   Float_t xx =x[0];
   Double_t delta;
   
   if( xx <= par[1] ) 
     delta =  (xx-par[1])/(par[2] * ( xx - par[1] ) + par[3]) ;
   else 
     delta =  (xx-par[1])/(par[4] * ( xx - par[1] ) + par[5]) ;
   
   Double_t f = par[0]*TMath::Exp(-0.5*delta*delta) ;
   return f;
}

Double_t ScintillationResponse(Double_t *x, Double_t *par)
{
   Float_t xx =x[0];
   Double_t f = (1.-TMath::Exp(-xx/par[0]))*TMath::Exp(-xx/par[1]);
   // Double_t f = 1.4*(TMath::Exp(-xx/par[1])-TMath::Exp(-xx/par[0]));
   return f;
}

class  GenWF{
 public :

  GenWF();
  GenWF(std::string jsonfile );
  virtual ~GenWF() {delete r; }
  void PrintSettings(void);

  void InitiateSupportFunctions(void);

  bool GeneratefromdEdx(double L,double energy); 
  bool Generate(double L,double nphotons);

  double GetVelocity(void){ return speed;}
  double GetScintFallTime(void) { return FallTime;}
  double GetBarLength(void) { return BarLength; } 

  // Counter getters
  double NPhotonsProduced(void) { return NPhotons; }
  double NPhotonsSimulated(void) { return NPhevt; }
  double NPhotonsAccepted(int endcode) { if( endcode == 0 ) return NPhcnt0; else return NPhcnt1; }


  // Reflection 
  double ReflFactor2surf(double cos0, double n0, double n1, double &cos1 );
  double ReflFactorScintAir(double cos0,double &cos1 );
  double ReflFactorAirSiPM(double cos0,double &cos1 );

  void computeMieScatteredDirection(double& z1, double& z2, double& z3 );

  void computeRayleighScatteredDirection(double& z1, double& z2, double& z3 );

  double sampleCosThetaHG(void);

  TH1F *GetPhotonTimeArrivalHist(int Endcode ) { if(Endcode == 0 ) return hpt0; else return hpt1; } 

  int TransportPhoton(double L, double timesamplingphase); 

  bool ElectronicsSimulation(int endcode );
  
  // This is the output
  std::vector<double> wfpub0;
  double Cell0Time_0; 

  std::vector<double> wfpub1;
  double Cell0Time_1; 

  // and the access method
  std::vector<double> GetWaveForm(int EndCode) {
    if( EndCode == 0 ) 
      return wfpub0;
    else if( EndCode == 1 )
      return wfpub1; 
    else
      {std::cout << " Invalid reaoud end code " << EndCode << std::endl; exit(0);}
  }
  
  double GetCell0Time(int EndCode) {
    if( EndCode == 0 ) 
      return Cell0Time_0;
    else if( EndCode == 1 ) 
      return Cell0Time_1; 
    else 
      { std::cout << " Invalid reaoud end code " << EndCode << std::endl; exit(0);}
  }
  
 private :

  TH1F *hpt0;
  TH1F *hpt1;
  
  // Geometry 
  double BarLength; // in m

  double Reflectivity3DMat; 
  
  // Scintillator 
  
  double n; //index of refraction;
  double nSi; //index of refraction of Silicon
  double nSiCover; //index of refraction of Silicon Cover
  double speed; // light velocity m/ns   
  double RaiseTime; // ns
  double FallTime; // ns 
  double Attenuation; // m
  double photonsperMeV;

  double gMie;
  double LambdaMie;
  double LambdaRayleigh;

  
  // MPPC's 
  double quantumefficiency;
  double fillfactor;
  double photocoverage; // Surface covered by the photosensitive are of SiPM
  double Sicoverage; // Surface covered by the Silicium in the SiPM.
  double PhotonEff;
  double Amplitude; 

  // Electronics
  
  int TotSamples;
  int PreSamples;
  int PostSamples; 
  double SamplingTime; // ns
  double threshold;
  double noise;
  double delaybtwSiPM;
  
  // Counters
  double NPhotons;
  double NPhevt;
  double NPhcnt0; // Photons detected in one end. 
  double NPhcnt1; // Photons detected in the other end. 

  // Mirror or internal refleciton. 
  bool Mirror;
  
  // WF
  double WF0;
  double WF1;
  double WF2;
  double WF3;
  double WF4;
  double WF5;
  
  double BallShapeTime; 
  

  TRandom *r;

  TF1 *fsc;
  TF1 *wfrf;
  double wfnorm;
  TF1 *Rayf; 
  
  TH1F *twf0; // Photon Arrival time for the z=0 end.
  TH1F *twfL; // Photon Arrival time for the z=L end.
  
};


double GenWF::ReflFactor2surf(double cos0, double n0, double n1, double &cos1 ) {
  double sin0 =sqrt(1-cos0*cos0);
  double sin1 = n0/n1*sin0;

  if( sin1 >=  1. ) return 1.; 

  cos1 = sqrt(1-sin1*sin1); 
  
  double rs = (n0*cos0-n1*cos1)/(n0*cos0+n1*cos1);
  double rp = (n1*cos0-n0*cos1)/(n1*cos0+n0*cos1);
  
  double R = (rs*rs+rp*rp)/2.; // reflection scint-air

  return R;
}


double GenWF::ReflFactorScintAir(double cos0,double &cos1 ) {
  double n0 =  n ; // Scintillator index of refraction 
  double n1 = 1.0; // Air index of refraction

   // Plastic-Air 
  double R = ReflFactor2surf(cos0,n0,n1,cos1); 

  return R; 
}

double GenWF::ReflFactorAirSiPM(double cos0,double &cos1 ) {
  double n0 = 1.0; // Air index of refraction
  double n1 = nSiCover; // Silicon Cover 
  double n2 = nSi; // Silicon index of refraction
  
  double cosloc;
  
   // Air-SiPMCover 
  double R1 = ReflFactor2surf(cos0,n0,n1,cosloc); 

  // SiPMCover-Si
  double R2 = ReflFactor2surf(cosloc,n1,n2,cos1); 

  
  return R1+(1.-R1)*R2; 
}


GenWF::GenWF(){

  // Some default configuration. 
  
  BarLength = 2.2; // in m 
  n = 1.58 ; //index of refraction Scint
  nSi = 4.0; // index of refraction Si
  nSiCover = 1.56; 
  speed = 299792458.e-9/n; // m/ns
  RaiseTime = 0.9; // ns
  FallTime = 2.1; // ns 
  

  gMie = 0.9;
  LambdaMie = 4.0;
  LambdaRayleigh = 3.8;
  Attenuation = 3.8; // m
  
  SamplingTime = 0.3125; // ns 
  Mirror = true;

  Reflectivity3DMat = 0.20; 
  
  quantumefficiency = 0.40;
  fillfactor = 0.74;
  photocoverage = (6.*6.)*8./10./120.;
  Sicoverage = (6.85*6.85)*8./10./120.;
  
  photonsperMeV = 10000.;

  PhotonEff = quantumefficiency*fillfactor*photocoverage/Sicoverage; // Photon efficiency within a SiPM sensor. 

  noise = 0.2; 

  delaybtwSiPM = 0.210; 
  
  threshold = 1.;

  TotSamples = 60;
  PreSamples = 15;

  PostSamples = TotSamples-PreSamples; 
  
  PrintSettings(); 
  
  r = new TRandom();

  InitiateSupportFunctions(); 
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
  speed = 299792458.e-9/n; // m/ns

  gMie = config["gMie"];
  LambdaMie = config["LambdaMie"];
  LambdaRayleigh = config["LambdaRayleigh"];

  nSiCover = config["IndexRefractionSiCover"]; //index of refraction;
  nSi = config["IndexRefractionSi"]; //index of refraction;
  RaiseTime = config["RaiseTime"]; // ns
  FallTime = config["FallTime"]; // ns 
  Attenuation = config["Attenuation"]; // m
  SamplingTime = config["SamplingTime"]; // ns
    
  Mirror = (bool) config["Mirror"]; // Mirror propagation

  Reflectivity3DMat = config["Reflectivity3DMat"]; // Reflectivity of the 3DPrinting material. 
    
  quantumefficiency = config["quantumefficiency"];
  fillfactor = config["fillfactor"];
  photocoverage = config["photocoverage"];

  Sicoverage = config["Sicoverage"]; // this is the fraction of the area covered by Si, larger than photo-coverage. 
  
  photonsperMeV = config["photonsperMeV"];

  PhotonEff = quantumefficiency*fillfactor*photocoverage/Sicoverage; // Photon efficiency within a SiPM sensor. 
  
  noise = config["noise"]; 

  delaybtwSiPM = config["delaybtwSiPM"]; 
  
  threshold = config["threshold"];

  TotSamples = config["TotSamples"];
  PreSamples = config["PreSamples"];

  Amplitude = config["Amplitude"];

  PostSamples = TotSamples-PreSamples; 

  WF0 = config["WF0"];
  WF1 = config["WF1"];
  WF2 = config["WF2"];
  WF3 = config["WF3"];
  WF4 = config["WF4"];
  WF5 = config["WF5"];

  
  BallShapeTime = config["BallShapeTime"];
  
  file.close();
  
  r = new TRandom();
  
  // This histogram will contain the waveform... units are samples not time units. 
  twf0 = new TH1F("twf0","  ",401,-0.5,400.5);
  twfL = new TH1F("twfL","  ",401,-0.5,400.5);
  
  // This histogram will contain the photon arrival times.
  hpt0 = new TH1F("hpt0","   ",1000,0.,100.);
  hpt1 = new TH1F("hpt1","   ",1000,0.,100.); 
  
  PrintSettings(); 

  InitiateSupportFunctions();

}


// Sample cos(theta) from Henyey-Greenstein distribution
double GenWF::sampleCosThetaHG(void) {
  double  xi = r->Uniform(0.,1.);

  if (std::abs(gMie) < 1e-6) {
    return 2.0 * xi - 1.0;  // Isotropic case
  }
  
  double term = (1. - gMie * gMie) / (1 - gMie + 2 * gMie * xi);

  double val = (1. + gMie * gMie - term * term) / (2 * gMie);
  
  return val; 
}

// Convert (theta, phi) to new direction vector assuming initial vector along z
void GenWF::computeMieScatteredDirection(double& z1, double& z2, double& z3 ) {

  double cosTheta = sampleCosThetaHG();
  double phi = r->Uniform(0.,2.*3.141592);
  
  double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);

  // some random y direction 
  double y1 = 0; double y2 = 1.; double y3 = 0;

  // get norm to d direction.
  double x1 = y2*z3 -z2*y3;
  double x2 = y3*z1 -z3*y1;
  double x3 = y1*z2 -z1*y2; 

  double nx = sqrt(x1*x1+y1*y1+z1*z1);
  x1 /= nx;x2 /= nx;x3 /= nx;
  
  y1 = z2*x3 - x2*z3;
  y2 = z3*x1 - x3*z1;
  y3 = z1*x2 - x1*z2;

  double sinphi = sin(phi);
  double cosphi = cos(phi); 
  
  z1 += x1*sinTheta*cosphi+y1*sinTheta*sinphi+z1*cosTheta;
  z2 += x2*sinTheta*cosphi+y2*sinTheta*sinphi+z2*cosTheta;
  z3 += x3*sinTheta*cosphi+y3*sinTheta*sinphi+z3*cosTheta;
   
  double nz = sqrt(z1*z1+z2*z2+z3*z3);
  z1/=nz;z2/=nz;z3/=nz;

  return;   
}


// Convert (theta, phi) to new direction vector assuming initial vector along z
void GenWF::computeRayleighScatteredDirection(double& z1, double& z2, double& z3 ) {

  double cosTheta = Rayf->GetRandom(-1.,1.); // r->Uniform(-1.,1.);
  double phi = r->Uniform(0.,2.*3.141592);
  double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
  double sinphi = sin(phi);
  double cosphi = cos(phi); 
  
    // some random y direction 
  double y1 = 0; double y2 = 1.; double y3 = 0;

  // get norm to d direction.
  double x1 = y2*z3 -z2*y3;
  double x2 = y3*z1 -z3*y1;
  double x3 = y1*z2 -z1*y2; 

  double nx = sqrt(x1*x1+y1*y1+z1*z1);
  x1 /= nx;x2 /= nx;x3 /= nx;
  
  y1 = z2*x3 - x2*z3;
  y2 = z3*x1 - x3*z1;
  y3 = z1*x2 - x1*z2;

  // Rotate the vector in the new direction.
  
  z1 += x1*sinTheta*cosphi+y1*sinTheta*sinphi+z1*cosTheta;
  z2 += x2*sinTheta*cosphi+y2*sinTheta*sinphi+z2*cosTheta;
  z3 += x3*sinTheta*cosphi+y3*sinTheta*sinphi+z3*cosTheta;
   
  double nz = sqrt(z1*z1+z2*z2+z3*z3);
  z1/=nz;z2/=nz;z3/=nz;
  
  return;   
}


void GenWF::PrintSettings(void) {

  std::cout << std::endl<< std::endl;
  std::cout << " MC photon generation settings " << std::endl<<std::endl;
  std::cout << " ----------------------------- " << std::endl;
  std::cout << " -------- Geometry  ---------- " << std::endl;
  std::cout << " ----------------------------- " << std::endl;
  std::cout << " Bar Length " << BarLength << " m " << std::endl;
  std::cout << " Bar active photo area coverage " << photocoverage << std::endl;
  std::cout << " Bar total Silicium area coverage " << Sicoverage << std::endl;
  std::cout << " ----------------------------- " << std::endl;
  std::cout << " --------   MPPC    ---------- " << std::endl;
  std::cout << " ----------------------------- " << std::endl;
  std::cout << " Quantum efficiency " << quantumefficiency << std::endl;
  std::cout << " Pixel fill factor " << fillfactor << std::endl;
  std::cout << " Total photon detection efficiency " << PhotonEff << std::endl;
  std::cout << " Amplitude " << Amplitude << std::endl;
  std::cout << " Index refraction Silicon " << nSi << std::endl;
  std::cout << " Index refraction Silicon Cover " << nSiCover << std::endl; 
  std::cout << " ----------------------------- " << std::endl;
  std::cout << " ------- Scintillator -------- " << std::endl;
  std::cout << " ----------------------------- " << std::endl;
  std::cout << " Index refraction " << n << std::endl; 
  std::cout << " Raise Time " << RaiseTime << " ns" << std::endl; 
  std::cout << " Fall Time " << FallTime << " ns" << std::endl;
  std::cout << " Photon speed " << speed <<  " m/ns" << std::endl;
  std::cout << " Mie g " << gMie << std::endl; 
  std::cout << " Mie Lambda " << LambdaMie << std::endl;
  std::cout << " Rayleigh Lambda " << LambdaRayleigh << std::endl;
  std::cout << " 3D Printing material Reflectivity " << Reflectivity3DMat << std::endl;
  std::cout << " Attenuation " << Attenuation << " m" << std::endl;
  std::cout << " # Photons " <<  photonsperMeV << " photons/MeV " << std::endl;
  std::cout << " ----------------------------- " << std::endl;
  std::cout << " ------ Electronics  --------- " << std::endl;
  std::cout << " ----------------------------- " << std::endl;
  std::cout << " Sampling time " << SamplingTime << std::endl;
  std::cout << " Electronics noise " << noise << std::endl;
  std::cout << " Trigger threshold " << threshold << std::endl;
  std::cout << " Waveform total samples " << TotSamples << std::endl;
  std::cout << " Waveform pre-samples " << PreSamples << std::endl;
  std::cout << " ----------------------------- " << std::endl;
  std::cout << " ------ WF parameters  --------- " << std::endl;
  std::cout << " ----------------------------- " << std::endl;
  std::cout << "  if ( x < " << WF1 << " ) " << std::endl;
  printf("      %3.2f*TMath::Exp(-(t-%3.2f)^2)/(2*(%3.2f*(t-%3.2f)+%3.2f)^2)\n",WF0,WF1,WF2,WF1,WF3) ;
  std::cout << "  elseif ( x >= " << WF1 << " ) " << std::endl;
  printf("      %3.2f*TMath::Exp(-(t-%3.2f)^2)/(2*(%3.2f*(t-%3.2f)+%3.2f)^2)\n",WF0,WF1,WF4,WF1,WF5) ;
  std::cout << " ----------------------------- " << std::endl;
}


void GenWF::InitiateSupportFunctions(void){
  //
  // This is the electronics response function for a fast imput light pulse. It includes all the electronics chain.
  //
  wfrf = new TF1("wfrf",WFResponsefunction,0.,300.,6);
  wfrf->SetParameters(WF0,WF1,WF2,WF3,WF4,WF5);

  wfnorm = 0;
  for(int j = 0; j < (int) wfrf->GetXmax();j++ ) {
    // The WF histogram contains the time with respect the particle t0 
    wfnorm += wfrf->Eval((double)j); 
  }
 
  //
  // This is the scintillator response of the EJ-200
  //
  fsc  = new TF1("fsc",ScintillationResponse,0.,30.,2);
  fsc->SetParameters(RaiseTime,FallTime); 

  Rayf = new TF1("Rayf","1+x*x",-1.,1.); 
  
}

bool GenWF::GeneratefromdEdx(double L, double E){

  double nphotons = photonsperMeV*E;
  
  return Generate(L,nphotons); 
}

// Transport a single photon and fill the time in the histogram twf0 and twfL for the two ends. 
// Inputs are the position and the common time for all photons.  
int GenWF::TransportPhoton(double position,double timesamplingphase) {

  // Simulate the photons.
  double d = 0.;
  double cosine = r->Uniform(-1.,1.);
  double sine = sqrt(1.-cosine*cosine); 
  double phi = r->Uniform(0.,2.*3.141592);
  double x = 0; double y = 0; double z = position; 
  double dirx,diry,dirz;
  dirx = sine*cos(phi);
  diry = sine*sin(phi); 
  dirz = cosine; 
  
  double accdist = 0; 
  
  bool absorbed = false;
  bool detected = false; 
  
  double photontime = 0.0;
  
  double Lambda = 1./(1./LambdaMie+1./LambdaRayleigh+1./Attenuation); // Average interaction length 

  double BarEnd = -9999; // This will flag the end at which the photon was absorved by the SiPM.  

  int nphotonsinend0 = 0; 
  
  // Transport the photon step by step.
  do{
    
    double distf = r->Exp(Lambda); 
    accdist += distf;
    
    x += dirx*distf;
    y += diry*distf; 
    z += dirz*distf;
    photontime += distf/speed;     
    
    if( z >= BarLength || z <= 0.0  ) { // The photon goes beyond the bar ends. 
      // Reflect the photons
      double cosint;
      // First reflection in Scint-Air transition 
      double Rsa = ReflFactorScintAir(abs(dirz),cosint);
      
      if( r->Uniform(0.,1.) > Rsa ) { // photon not reflected in Scint-Air
	if(  r->Uniform(0.,1.) < Sicoverage ) { // Is in the SiPM acceptance ?
	  double cosint2; 
	  double RaSPM = ReflFactorAirSiPM(cosint,cosint2); // Is reflected by the SIPM cover ?
	  
	  if( r->Uniform(0.,1.) > RaSPM ) { // Photon absorbed by SiPM.
	    double correctedistf;	    
	    detected = true;
	    if( z <= 0 )
	      BarEnd = 0;
	    else
	      BarEnd = 1;
	    // Correct the time.
	    if( z <= 0 ) {
	      correctedistf = TMath::Abs(z/dirz);
	    }
	    else if( z >= BarLength ) {
	      correctedistf = TMath::Abs((z-BarLength)/dirz);
	    }
	    else {
	      std::cout << " Not Possible " << std::endl; exit(0);
	    }
	    photontime -= correctedistf/speed;
	    break;
	  } 
	  else { // Photon is reflected by silicon 
	    double correctedistf;	    
	    // Move the photon back to the scintillator
	    if( z <= 0 ) {
	      correctedistf = TMath::Abs(z/dirz);
	      z = 0;
	    }
	    else if( z >= BarLength ) {
	      correctedistf = TMath::Abs((z-BarLength)/dirz);
	      z=BarLength;
	    }
	    else {
	      std::cout << " Not Possible " << std::endl; exit(0);
	    }
	    // Asume pure specular reflexion, no change in cos. 
	    dirz        = -dirz;
	    // Remove the time from expectation to real. 
	    photontime -= correctedistf/speed; 
	  }
	}
	else { // Outside SiPM acceptance

	  // Check the absorption in the 3D printed material
	  
	  if( r->Uniform(0.,1.) > Reflectivity3DMat ) { // Absorbed by plastic around SiPM
	    absorbed = true;
	    break;
	  }
	  else { 
	    double correctedistf;
	    
	    // Asume pure diffused reflexion. 
	    cosine = r->Uniform(0.,1.);  // randomised the direction, positive or negative direction decided later. 
	    sine = sqrt(1.-cosine*cosine);
	    
	    // The refraction entering from the air into the scintillator. 
	    sine = sine/n; 
	    cosine = sqrt(1.-sine*sine);
	    
	    phi = r->Uniform(0.,2.*3.141592);
	    dirx = sine*cos(phi);
	    diry = sine*sin(phi); 
	    // Position the photon back to the scintillator
	    if( z <= 0 ) {
	      correctedistf = TMath::Abs(z/dirz);
	      z    =   0;
	      dirz = cosine;  // it goes forwards. Done after time correction. 
	    }
	    else if( z >= BarLength ) {
	      correctedistf = TMath::Abs((z-BarLength)/dirz);
	      z    =   BarLength;
	      dirz = -cosine;  // it goes backwards, negative cosine. Done after time correction. 
	    }
	    else {
	      std::cout << " Not Possible " << std::endl; exit(0);
	    }
	    photontime -= correctedistf/speed;
	  }
	}
      }
      else { // Photon reflected in the Scint-Air gap
	double correctedistf; 
	if( z <= 0 ) {
	  correctedistf = TMath::Abs(z/dirz);
	  z = 0;
	}
	else if( z >= BarLength ) {
	  correctedistf = TMath::Abs((z-BarLength)/dirz);
	  z=BarLength;
	}
	else {
	  std::cout << " Not Possible " << std::endl; exit(0);
	}
	photontime -= correctedistf/speed;
	dirz = -dirz;
      }
    }
    else {  // Photon was not absorbed or detected or suffered scattering. 
    
      double selectprocess = r->Uniform(0.,1.);
      // Mie Scattering 
      if( selectprocess < Lambda/LambdaMie )
	computeMieScatteredDirection(dirx,diry,dirz);
      else if( (1.-selectprocess) < Lambda/LambdaRayleigh ) 
	computeRayleighScatteredDirection(dirx,diry,dirz);
      else {
	absorbed = true;
	break;
      }
    }
  } while( !absorbed && !detected );
  
  if( detected ) {    

    //Add the sampling phase bewteen the electronics clock and the event arrival. 
    photontime += timesamplingphase;
    
    // Add the scintillator de-excitation. 
    photontime += fsc->GetRandom();

    if( BarEnd == 0 ) 
      hpt0->Fill(photontime);
    else if( BarEnd == 1 ) 
      hpt1->Fill(photontime);

    // Take into account the time difference from the 4 blocks of SiPM 
    int ngroup = r->Uniform(0.,4.); // From 0 to 4 !
    switch(ngroup) {
    case 0:  // The readout is centered bewteen these two 
    case 1: 
      photontime += 0;
      break;
    case 2:
      photontime += delaybtwSiPM; 
      break;
    case 3:
      photontime += 2.* delaybtwSiPM; 
      break;
    }
    
    // Sample arrival time bin. 
    int tj0 = (int)(photontime/SamplingTime);
    
    double offset = (photontime-tj0*SamplingTime)/SamplingTime; // offset inside the sampling for the signal integration
    if( TMath::Abs(offset) > 1 ) std::cout << " This is strange " << offset << "  " << SamplingTime << std::endl;
    
    // Add the electronics shaper
    for(int j = 0; j < (int) wfrf->GetXmax();j++ ) {
      // The WF histogram contains the time with respect the particle t0 
      double tj = tj0+j;	
      double ftj = wfrf->Eval(((double)j+offset)*SamplingTime)/wfnorm;
      if( BarEnd == 0 ) 
	twf0->Fill(tj,ftj);
      else if ( BarEnd == 1 )
	twfL->Fill(tj,ftj);
      else
	std::cout << " Error, photon can be accepted only in one of the two ends " << BarEnd << std::endl;
    }

    if( BarEnd == 0 )
      return 0;
    else if (BarEnd == 1 )
      return 1;
  }

  return -1;  
}


bool GenWF::Generate(double L, double nphotons) {
  NPhotons = nphotons; 
  NPhcnt0 = NPhcnt1 = 0;
    
  // Photon detection efficiency 
  nphotons *= PhotonEff;

  twf0->Reset();
  twfL->Reset();
  
  double mincos = 0.; // this is 0 for specular reflection. 

  if( !Mirror ) 
    mincos = TMath::Sqrt(1.-1/n*1/n); // there is no photons below this value if full reflection in the bar.  

  double nphotonseff = nphotons*(1.-mincos);   // Effective number of photons. 
  
  int nphevt = r->Poisson(nphotonseff); // Fluctuate the number of photons. 

  // Constant for the event. This takes into consideration the fact that signals and clocks are not in synch normally. 
  double timesamplingphase = r->Uniform(0.,SamplingTime); 

  // Some counters 
  int reflected = 0; 
  int nrejected = 0;  

  NPhevt = nphevt;
  
  for(int i = 0;i < nphevt; i++ ) { 
    double indx =  TransportPhoton(L, timesamplingphase);
    if(indx < 0 )
      nrejected++;
    else if ( indx == 0 )
      NPhcnt0++;
    else if ( indx == 1 )
      NPhcnt1++;
  }
 
  // end loop one event
  // Generate the wave forms for both ends  
  
  bool goodsignal = ElectronicsSimulation(0); // End close to entry point. 

  goodsignal &= ElectronicsSimulation(1); // End farther from entry point. 

  return goodsignal;
}
  
bool GenWF::ElectronicsSimulation(int endcode ) {
  
  std::vector<double> wf;
  
  int ith = 0;
  bool found = false;

  TH1F *twf;

  if( endcode == 0 )
    twf = twf0;
  else 
    twf = twfL;
  
  // Add noise and get the threshold 
  for(int i = 0; i < twf->GetNbinsX();i++) {
    double val = Amplitude*twf->GetBinContent(i+1)+r->Gaus(0.,noise);
    // First time we cross the threshold
    if( val > threshold && !found ) {
      ith = i;
      found = true;
    }
    wf.push_back(val); 
  }
  
  if( !found ) return false; 
  
  // Copy the valid section of the WF 
  std::vector<double> wfpub(&wf[ith-PreSamples],&wf[ith+PostSamples]); 
  

  // Recompute the Cell0Time to the beginning of the WF using the time clock of the Cell0Time.
  if( endcode == 0 ) { 
    Cell0Time_0 = SamplingTime*(double)(ith-PreSamples);
    wfpub0.clear();
    wfpub0 = wfpub; 
  }
  else {
    Cell0Time_1 = SamplingTime*(double)(ith-PreSamples);
    wfpub1.clear();
    wfpub1 = wfpub;
  }
  
  return true; 
}


#endif 
 
