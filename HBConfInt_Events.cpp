#include <iostream>
#include "TInterpreter.h"
#include "TCanvas.h"
#include <TROOT.h>
#include "TSystem.h"
#include "TFile.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TFrame.h"
#include <TBenchmark.h>
#include "TCint.h"
#include "TStyle.h"
#include "TRandom.h"
#include <include/TBoundedGaussian.hh>
#include <include/Statistics.hh>


//Double_t Nexcess;


int HBConfInt_Events(int Non,int Noff,Double_t Alpha, double conf_level=95.){

  // Calculate flux error - from: MRSpectrum.C !!
  // For N < 20 use (hopefully ;-) correct poisson statistics
  
  //   double fluxErr(0.);
  //   if(flux > 20) {
  //     fluxErr = std::sqrt(flux);
  //   } else if(flux > 0) {
  //     double down = flux - Utilities::Statistics::ErrorDown((int)flux);
  //     double up = Utilities::Statistics::ErrorUp((int)flux) - flux;
  //     fluxErr = std::sqrt(down > up ? down : up);
  //   }
  
  //   sigma_flux = fluxErr*1.;
  

  Float_t Nexcess=1.;

  Utilities::TBoundedGaussian fc(conf_level/100.); // (99%) confidence level
  double upper(0),lower(0);

  std::cout << "Measured Non is: " << Non << std::endl;
  std::cout << "Measured Noff is: " << Noff << std::endl;
  std::cout << "Alpha is: " << Alpha << std::endl;

  Nexcess = (float)Non - (Alpha*(float)Noff);
  //  Nexcess = 3.14;

  std::cout << "Nexcess is Non-Alpha*Noff = " << Nexcess << endl << endl;


 // Try derivation directly from # events:
 // int Utilities::TBoundedGaussian::GetConfidenceInterval(int Non,int Noff,Double_t norm,Double_t& min,Double_t& max,int flag){
 //

 // Some values:
 // fc.GetConfidenceInterval(292.,3229.,0.0909,lower,upper);
 // fc.GetConfidenceInterval(292.,2788.,0.1004,lower,upper);
 // fc.GetConfidenceInterval(26658.,222453.,0.1192,lower,upper);


 fc.GetConfidenceInterval(Non,Noff,Alpha,lower,upper);


 if( lower<=0 )
   std::cout << conf_level << "% CL upper limit on Nevents: " << upper << std::endl;
 // else if (upper < Non)
 else {
   std::cout << "Confidence interval for true Nexcess: " << lower << " .. " << upper << std::endl;
   std::cout << "Nexcess is (at " << conf_level << "% confidence level): " << std::endl;
   std::cout << Nexcess << " -" << Nexcess-lower << " +" << upper-Nexcess << std::endl;
 }


 return ;

}
