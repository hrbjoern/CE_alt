/*-------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------*/
// ...
/* ... hbcalceffarea.cpp ... 2010-09-21 ff.
// ... 2011-02-09: Start einer umfassenden Erweiterung, inklusive Analyse?
// ...             HBCalcLimits.cpp !!!   
// ...
// ...
// ... Note: Everything's based on MRAnalyseLightcurve2 
// ...       and mrspectraltools/src/MRLightcurve2::CalculateResults                    
// ...
// ...
/* So, let's see ... calculate eff. areas from MRAnalysis output ROOT file.      */
/*-------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------*/

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>

#include <TROOT.h>
#include <TKey.h>
#include <TGaxis.h>
#include <TMath.h>
#include <TDatime.h>
#include "TSystem.h"
#include "TH1.h"
#include <TStyle.h>
#include <TApplication.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>

// #include <crashdata/UTC.hh>
// #include <utilities/TBoundedGaussian.hh>

//#include </nfs/astrop/d1/hess/software/hap-10-06/mrdata/include/MRDataFile.hh>
#include <mrdata/MRDataFile.hh>                                         // mit ++ laeuft's ...
#include </nfs/astrop/d1/hess/software/hap-10-06/mrdata/include/MRDataRun.hh>
#include </nfs/astrop/d1/hess/software/hap-10-06/mrdata/include/MREvent.hh>

// #include <mrtools/MRStatistics.hh>
// #include <mrtools/MRDBReader.hh>

#include <mrtools/MRConfig.hh>
#include <mrtools/MRStyle.hh>
#include <mrtools/MRSkymapTools.hh>

#include <mrplotters/MRSkymapPlotter4.hh>

// #include <mrdata/MRDataSelector.hh>

// #include <mrrecotools/MRMuonCorrection.hh>

#include </nfs/astrop/d1/hess/software/hap-10-06/mrspectraltools/include/MREffAreaLookup.hh>
//#include <mrspectraltools/MRLightcurve2.hh>
// #include <map>



int HBCalcLimits(std::string fileName){

  // mrscripts/skymap/PlotSkymap4.C :
  double overSRadius = 0.1225;
  double exRadiusFact = 1.5;
  std::string exSkyRegFileName = "";
  std::string outputFile = "";
  bool doOverwrite = 0;


  MRTools::SetStandardStyle();
  //  TStopwatch timer;

  std::cout << "INFO: Plotting skymap from file : " << fileName << std::endl;
  std::cout << "INFO: Oversampling radius is set to " << overSRadius << " deg" <<std::endl;
  std::cout << "INFO: Exclusion radius around the source is " << overSRadius*exRadiusFact << " deg" << std::endl;

  TFile* f = new TFile(fileName.c_str());

  if((f == 0) || f->IsZombie()) {
    std::cout << "ERROR: Could not open file : " << fileName << std::endl;
    return 0; 
  }

  TDirectory *skymapDir = dynamic_cast<TDirectory*>(f->Get("Skymap"));
  if((skymapDir == 0) || skymapDir->IsZombie()) {
    std::cout << "ERROR: Could not get directory Skymap from file " << fileName << std::endl;
    return 0;
  }

  gROOT->cd();

  MRTools::MRExcludeSkyRegions* exSkyRegions = new MRTools::MRExcludeSkyRegions();
  TMarker* marker(0);
  skymapDir->GetObject("TMarker", marker);
  if(marker != 0) {
    exSkyRegions->AddSkyRegionCircle(marker->GetX(),
				     marker->GetY(),
				     overSRadius*exRadiusFact);
  }
  if(exSkyRegFileName != "") exSkyRegions->ReadSkyRegionsFromFile(exSkyRegFileName);

  exSkyRegions->Print();

  MRPlotters::MRSkymapPlotter4* skymap = new MRPlotters::MRSkymapPlotter4(0, 0);
  skymap->SetExcludeSkyRegions(exSkyRegions);
  //skymap->SetUseGalCoord(1);

  double skymapSizeDeg(10.);// 5.
  int nBins(100);// 50
  //int nBinsOvers(200);// 200
  int nBinsOvers(400);// 200
  double ringMin(0.2);// 0.2
  double runMax(0.3);// 0.5
  //double ringMin(0.2);// 0.2
  //double runMax(0.3);// 0.5
  double borderDeg(.6);// 1.

  skymap->SetSkymapParameters(skymapSizeDeg, nBins, nBinsOvers, overSRadius, ringMin, runMax, borderDeg);
  skymap->SetSmoothingRadiusDeg(0.1);
  //skymap->SetSmoothingRadiusDeg(0.1414);
  skymap->SetDoDZAngleCorrection(1);
  skymap->SetUseLinDZAngleCor(1);
  skymap->SetDoGausOversample(0);
  //skymap->SetDoCamAccCorrection(false);

  std::cout << "Reading trees ..." << std::endl;
  if(!skymap->ReadTreesFromDir((TDirectory*)f->Get("Skymap"), false, true, true, true)) {
    std::cout << "ERROR: Failed reading trees" << std::cout;
    return 0;
  }

  std::cout << "Calculating results ..." << std::endl;
  skymap->DrawResults();

//   if(outputFile != "") {
//     std::cout << "INFO: Writing data to file : " << outputFile << std::endl;
//     if(!gSystem->AccessPathName(outputFile.c_str(), 0) && !doOverwrite) {
//       std::cout << "ERROR: File already exists : " << outputFile << std::endl;
//     } else {
//       TFile *outF = new TFile(outputFile.c_str(), "UPDATE");
//       outF->cd();
//       outF->Delete("Skymap;*");
//       outF->Write();
//       TDirectory *dir = outF->mkdir("Skymap", "Skymap");
//       skymap->WriteResults(dir);
//       dir->Write();
//       outF->Write();
//       outF->Close();
//       delete outF;
//     }
//     std::cout << "INFO: Writing done" << std::endl;
//   }

  f->Close();

  //.........................................................//

  //int HBCalcEffArea_Test(){

  // Data file:
  //  TFile* f = new TFile(fileName.c_str());

  // Variables and flags:
  double energyThreshold = 0.264;
  double spectralIndex = -2.7;
  double theta2 = 0.01;
  double maxImpactPar = 800.;
  // rate correction factor: eventuell noch einbauen?!
  double rateCorFactor(1.);
  bool exRunMethod3(0);
  std::ostringstream exRunNMethod3;
  int effAreaCalcMode = 0;
  int effAreaNBins = 30;
  bool beVerbose = 0;
  double cor;
  double cor1 = 1.;
  bool TEST(1);
  TEST=0;

  // Eff. area normalization factor:
  double fluxFactor = 1./(std::pow(maxImpactPar, 2)*TMath::Pi()*1.e4);

  // ROOT and Histos:
  gROOT->Reset();

  TH1D *EffAreasPerEvent = new TH1D("EffAreasPerEvent","Eff. Areas Per Event",35,-1.,1.5); 
  TH1D *Events = new TH1D("Events","Events", 35,-1.,1.5); 
  TH1D *EffAreaNormalized = new TH1D("EffAreaNormalized","Eff. Area, total, normalized;log(E);Aeff(E)",35,-1.,1.5); 


  // From MRAnalyseLightcurve2: 
  /*-------------------------------------------------------------------------------*/
  // Handle Effective Areas
  /*-------------------------------------------------------------------------------*/
  // read in effective areas for 4 tel
  MRSpectralTools::MREffAreaLookup* effArea = new MRSpectralTools::MREffAreaLookup();
  effArea->ReadEffAreasFromFile("/d1/hess/config/mr/lookup/effarea/config/effArea_south_hdLowEffMC_0510_s20_newELookup_4tel.txt",0,30,0);
  if(TEST==0)
    effArea->CalculateEffectiveArea(spectralIndex, theta2, maxImpactPar, 0, beVerbose);


  // Read in data file:
  //  MRData::MRDataFile* dataFile = new MRData::MRDataFile("2010-09-07_v3.root");
  MRData::MRDataFile* dataFile = new MRData::MRDataFile((fileName.c_str()));
  MRData::MRDataRun* dataRun = 0;

  
  // Run loop:
  while((dataRun = dataFile->NextDataRun())){
    int runN = dataRun->GetInt("RunNumber");

    //    std::cout << "bin noch da \n" << std::endl;


    // Event loop: (from MRLightcurve)
    TTree* tree = dataRun->GetEventTree();
    MRData::MREvent* event = dataRun->GetEvent();
    int nEntries = (int)dataRun->GetEventTree()->GetEntries();
    bool firstEvent(true);
     
    // CHB TEST 2011-02-07:
    Int_t Esubth=0, Eplusth=0;

    for(int i = 0; i < nEntries; i++) {
      
      tree->GetEntry(i);

      //	  if(!fDataSelector || fDataSelector->EventPassCuts(event)) {
      double impactPar = event->GetImpactParM();
      double energy = event->GetEnergyTeV()*rateCorFactor;
      //      if(fUseMuonCorrection) energy = event->GetMCEnergyTeV()*rateCorFactor;

      if((energy > 0.) && (impactPar < maxImpactPar)) {
	//	std::cout << "bin noch da \n" << std::endl;
	int eventType = event->GetEventType();

	// 	//      if((eventType == 0) || (eventType == -10)) {  // CHB: Fill Histos! TESST!
	// 	EffAreasPerEvent->Fill(TMath::Log10(energy),cor1/fluxFactor);
	// 	//	EffAreasPerEvent->Fill(TMath::Log10(energy),cor1);
	// 	Events->Fill(TMath::Log10(energy));
	// 	EffAreaNormalized->Fill(TMath::Log10(energy),cor1/fluxFactor);
	
	// 	std::cout << "energy,cor1/fluxFactor = " << energy << "    " << cor1/fluxFactor << std::endl;
	// 	std::cout << "energy,fluxFactor = " << energy << "    " << fluxFactor << std::endl;
	//      }

	// CHB TEST 2011-02-07:
	if(energy < energyThreshold) {
	  std::cout << "HALLO! E < Eth" << endl;
	  Esubth++;
	  std::cout << "Esubth = " << Esubth << endl;
	}

	if(energy > energyThreshold) {
	  // CHB TEST 2011-02-07:
	  std::cout << "HALLO! E > Eth" << endl;
	  Eplusth++;
	  std::cout << "Eplusth = " << Eplusth << endl;
	  
	  if(!exRunMethod3) {
	    double camOffset = event->GetShowerCamDistDeg();
	    double cosZAngle = std::cos(event->GetZAngleDeg()*TMath::Pi()/180.);
	    cor = 
	      //	      fEffAreaLookup3Tel
	      //	      ? fEffAreaLookup3Tel->GetCorrection(event->GetEnergyTeV(), cosZAngle, camOffset)
	      //	      : 
	      effArea->GetCorrection(event->GetEnergyTeV(), cosZAngle, camOffset); 
	    // CHB: Wichtig! Bisher nur 4 Tel, und: Handling?!
	    if(cor != 0.) {
	      //	      std::cout << "cor = " << cor << endl;
 	      if((eventType == 0) || (eventType == -10)) {  // CHB: Fill Histos!
		EffAreasPerEvent->Fill(TMath::Log10(energy),cor/fluxFactor);
		//	EffAreasPerEvent->Fill(TMath::Log10(energy),cor1);
		Events->Fill(TMath::Log10(energy));
		EffAreaNormalized->Fill(TMath::Log10(energy),cor/fluxFactor);
	      }
// 		fFluxCounterMap["runM3"].AddOnCor(fluxFactor/cor);
// 		if(fDoMinLightcurve) fFluxCounterMap["minM3"].AddOnCor(fluxFactor/cor);
// 	      } else if(eventType == 2) {
// 		fFluxCounterMap["runM3"].AddOffCor(alphaFactor*fluxFactor/cor, alphaFactor);
// 		if(fDoMinLightcurve) fFluxCounterMap["minM3"].AddOffCor(alphaFactor*fluxFactor/cor, alphaFactor);
// 	      }
	    } else {
	      std::cout << "WARNING <MRLightcurve2::CalculateResults>" << std::endl
			<< "    -> cor = 0. for method 3" << std::endl
			<< "    -> NOT REALLY Excluding run #" << runN << " for method 3" << std::endl;
	      //	      exRunMethod3 = true;
	      exRunNMethod3 << " " << runN;
	    }
	  }
	}// END if(energy > fEnergyThreshold)

      } // anderes if?!

    } // end of event loop
    
  } // end of run loop
    
  std::cout << "So, jetzt haben wir den Event, seine Energie und sowas wie eine eff. Fläche ... Also brauchen wir die Histos!\n" << std::endl;
  
  
  //EffAreasPerEvent->DrawCopy();
  EffAreaNormalized->Divide(Events);
  //  EffAreaNormalized->GetYAxis->SetDraw();
  TCanvas *c1 = new TCanvas();
  gPad->SetLogy();
  EffAreaNormalized->Draw();
  //      Events->Draw();
  
  
  // Write eff. area data to output file: 
  ofstream outfile("HBCalcEffArea.out");
  for(int i=1;i<=(EffAreaNormalized->GetNbinsX());i++) {
    outfile << "log(E),Aeff(E)= " << EffAreaNormalized->GetXaxis()->GetBinCenter(i) 
	    << " " << EffAreaNormalized->GetBinContent(i) << endl; 
  }

  outfile.close();

  


  // THE END:
  return 0 ;



}

