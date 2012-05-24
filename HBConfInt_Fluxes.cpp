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
#include </nfs/astrop/d1/hess/software/hap-10-06/utilities/include/TBoundedGaussian.hh>
#include </nfs/astrop/d1/hess/software/hap-10-06/utilities/include/Statistics.hh>


int HBConfInt_Fluxes(double flux, double sigma_flux, double conf_level=99.){

  // Calculate flux error - from: MRSpectrum.C !!
  // For N < 20 use (hopefully ;-) correct poisson statistics

//   double fluxErr(0.);
//   if(flux > 20) {
//     fluxErr = std::sqrt(nON);
//   } else if(flux > 0) {
//     double down = flux - Utilities::Statistics::ErrorDown((int)flux);
//     double up = Utilities::Statistics::ErrorUp((int)flux) - flux;
//     fluxErr = std::sqrt(down > up ? down : up);
//   }
  


 Utilities::TBoundedGaussian fc(conf_level/100.); // (99%) confidence level
 double upper(0),lower(0);
 // double flux = 100.;
 // double sigma_flux = 10.;

 std::cout << "Measured flux is: " << flux << std::endl;
 std::cout << "Error of flux is: " << sigma_flux << std::endl;
 fc.GetConfidenceInterval(flux,sigma_flux,lower,upper);
 if( lower<=0 )
   std::cout << conf_level << "% CL upper limit on flux: " << upper << std::endl;
 else {
   std::cout << "Confidence interval for true flux: " << lower << " .. " << upper << std::endl;
   //   std::cout << "Flux is (at " << fc.GetCL()*100.0 << " confidence level): " << std::endl;
   std::cout << "Flux is (at " << conf_level << "% confidence level): " << std::endl;
   std::cout << flux << " -" << flux-lower << " +" << upper-flux << std::endl;
 }

 return;

}
