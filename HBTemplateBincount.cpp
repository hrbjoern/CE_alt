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

// using std::cout;
// using std::cin;
// using std::endl;

// HB's variables:
Double_t xDist, yDist, rDist;
Double_t nHadrons=0.;
Double_t nGammas=0.;
//void draw();   

void HBTemplateBincount(){

  //Preliminaria:
  gROOT->Reset();
  //  gStyle->SetPalette(10,0);

  // Declare Histos:
  //  static  TH2D *h1new = new TH2D();
  TH2D *h1corrnew = new TH2D();
  TH2D *p1new = new TH2D();
  TH2D *s1new = new TH2D();

  //chb test:
  cout << "class name 1 = " << h1corrnew->ClassName() << endl;

  // Open input ROOT file:
  //  TFile *f = TFile::Open("FornaxClusterWobble_2010-05-27_template1.root");
  //  TFile *fi = new TFile("FornaxClusterWobble_2010-05-27_template1.root");
  
  //  TFile *fi = new TFile("2010-07-19_Test.root");
  //  TFile *fi = new TFile("2010-08-24_Test2.root");
  //  TFile *fi = new TFile("2010-08-24_Test3.root");
  //  TFile *fi = new TFile("2010-09-07_PlotSkymap.root");
  //  double alpha = 0.1192;
  //  TFile *fi = new TFile("2010-09-15_PlotTemplateSkymap.root");
  //  TFile *fi = new TFile("2010-09-27_SkymapOutput.root");
  //  double alpha = 0.117536;
  //  TFile *fi = new TFile("2010-10-04_point9degrees.root"); // 0.9 deg int radius; for subtraction!?
  //  double alpha = 0.117571.;
  //  TFile *fi = new TFile("2010-10-04_onedegree.root"); // 1.0 deg int radius; both from 0.5 deg Wobble root file!
  //  double alpha = 0.117615;
  //  TFile *fi = new TFile("2010-10-04_point1degrees.root"); // 0.1 deg int radius; 
  //  double alpha = 0.117688.;
  //  TFile *fi = new TFile("2010-10-04_point9degrees_v2.root"); // 0.9 deg int radius; 
  //  double alpha = 0.117571.;
  //  TFile *fi = new TFile("2010-10-04_onedegree_v2.root"); // 1.0 deg int radius; both from 0.5 deg Template root file
  //  double alpha = 0.117615;
  //  TFile *fi = new TFile("2010-10-05_point99degrees.root"); // 0.99 deg int radius; 
  //  double alpha =   0.117718;
  //  TFile *fi = new TFile("2010-10-06_onedegree_hardcuts.root"); // 1.0 deg int radius; from 0.1 deg hard cuts file
  //  double alpha = 0.0595231;
  TFile *fi = new TFile("2010-10-06_point9degrees_hardcuts.root"); // 1.0 deg int radius; from 0.1 deg hard cuts file
  double alpha = 0.0591913;


  fi->ls();
  //  f.cd("Skymap");
  //  gDirectory->ls();
//   //  MRSkymapPlotter4_skymap_h1->Draw();

  TDirectory* dir;
  //  dir = dynamic_cast<TDirectory*>(fi.Get("Skymap"));
  dir = dynamic_cast<TDirectory*>(fi->Get("Skymap"));
  //  dir->ls();


  // Go back to basic ROOT directory:
//  gDirectory->pwd();
  gROOT->cd();
  gDirectory->pwd();

  // Copy template skymap histos:
  // Hadron candidates:
  h1corrnew = dynamic_cast<TH2D*> ( (dynamic_cast<TH2D*>(dir->Get("MRSkymapPlotter4_skymap_h1_overs_cor")))->Clone("h1corrnew"));
  // Photon candidates:
  p1new = dynamic_cast<TH2D*> ( (dynamic_cast<TH2D*>(dir->Get("MRSkymapPlotter4_skymap_p1_overs")))->Clone("p1new"));
  // Template Excess High Res. Skymap:
  //  s1new = static_cast<TH2D*> ( (static_cast<TH2D*>(dir->Get("MRSkymapPlotter4_skymap_excess_template_highRes_smooth")))->Clone("s1new"));
  // Template Excess Oversampled:
  s1new = static_cast<TH2D*> ( (static_cast<TH2D*>(dir->Get("MRSkymapPlotter4_skymap_excess_overs")))->Clone("s1new"));
  // Template Excess Skymap:
  //  s1new = static_cast<TH2D*> ( (static_cast<TH2D*>(dir->Get("MRSkymapPlotter4_skymap_excess")))->Clone("s1new"));

  //  s1new = static_cast<TH2D*> ( (static_cast<TH2D*>(dir->Get("MRSkymapPlotter4_skymap_excess_overs")))->Clone("s1new"));


  //  h1new->SetTitle("h1new");
  //  h1new = dynamic_cast<TH2D*> ( (dynamic_cast<TH2D*>(dir->Get("MRSkymapPlotter4_skymap_h1")))->Clone("h1new")); 

  //  h1new->Print();

  //  cout << "class name MR = " << static_cast<TH2D*>(dir->Get("MRSkymapPlotter4_skymap_h1"))->ClassName() << endl;
  //  cout << "class name 2 = " << h1new->ClassName() << endl;
  //  cout << "Test 0" << endl << endl;

  // Close input ROOT file:
  fi->Close();

  gDirectory->pwd();
  s1new->Print();
  // cout << "Test 1" << endl << endl;

  // chb test:
  //  cout << "class name 3 = " << s1new->ClassName() << endl;
  //  h1new->DrawCopy();





  // Open output ROOT file:
  TFile *fo = new TFile("HBTemplateBincount_Out.root","RECREATE");
   
   //   cout << "Test 1" << endl << endl;
   //   fo->Write();
   //   cout << "Test 4" << endl << endl;
   //   fo->ls();

   s1new->Write();

   //   cout << "Test 5" << endl << endl;

   //   fo->ls();
   //   fo->Close();

   //   cout << "Test 6" << endl << endl;


   //*... Starting Fornax Analysis..................................................................................................................*//

   // Signal map: 

   cout  << "s1new->GetBinContent(24,24) = " <<    s1new->GetBinContent(24,24)  << endl;
   cout  << "s1new->GetBinContent(24,25) = " <<    s1new->GetBinContent(24,25)  << endl;
   cout  << "s1new->GetBinContent(25,24) = " <<    s1new->GetBinContent(25,24)  << endl;
   cout  << "s1new->GetBinContent(25,25) = " <<    s1new->GetBinContent(25,25)  << endl;

   cout  << "s1new->GetNbinsX() = " <<    s1new->GetNbinsX()  << endl;
   cout  << "s1new->GetNbinsY() = " <<    s1new->GetNbinsY()  << endl << endl;

   cout  << "s1new->GetBinContent(0,0) = " <<    s1new->GetBinContent(0,0)  << endl << endl;
   cout  << "s1new->GetBinContent(1,1) = " <<    s1new->GetBinContent(1,1)  << endl << endl;

   cout  << "s1new->GetBinContent(100,100) = " <<    s1new->GetBinContent(100,100)  << endl;
   cout  << "s1new->GetBinContent(100,101) = " <<    s1new->GetBinContent(100,101)  << endl;
   cout  << "s1new->GetBinContent(101,100) = " <<    s1new->GetBinContent(101,100)  << endl;
   cout  << "s1new->GetBinContent(101,101) = " <<    s1new->GetBinContent(101,101) << endl  << endl;

   Double_t Summe = s1new->GetBinContent(100,100)
     + s1new->GetBinContent(100,101)
     + s1new->GetBinContent(101,100)
     + s1new->GetBinContent(101,101);
   
   Double_t Summe2 =   s1new->GetBinContent(24,24) 
     +   s1new->GetBinContent(24,25) 
     +  s1new->GetBinContent(25,24) 
     +  s1new->GetBinContent(25,25) ;
   

 cout << "Summe/4. = " << Summe/4. << " events" << endl << endl;
 // cout << "Summe2/4. = " << Summe2/4. << " events" << endl << endl;



 // Draw 2D Histo:
 // s1new->DrawCopy("COLZ");


 //
 // 'corrected' hadrons and gamma cand/s: 
 //

 cout  << "h1corrnew->GetBinContent(100,100) = " <<    h1corrnew->GetBinContent(100,100)  << endl;
 cout  << "h1corrnew->GetBinContent(100,101) = " <<    h1corrnew->GetBinContent(100,101)  << endl;
 cout  << "h1corrnew->GetBinContent(101,100) = " <<    h1corrnew->GetBinContent(101,100)  << endl;
 cout  << "h1corrnew->GetBinContent(101,101) = " <<    h1corrnew->GetBinContent(101,101) << endl  << endl;

 Double_t Summe_hcorr = h1corrnew->GetBinContent(100,100)
   + h1corrnew->GetBinContent(100,101)
   + h1corrnew->GetBinContent(101,100)
   + h1corrnew->GetBinContent(101,101);

 cout << "Summe_hcorr = " << Summe_hcorr <<" \n" <<endl;
 cout << "Summe_hcorr/4. = " << Summe_hcorr/4. <<" \n" <<endl;

 cout  << "p1new->GetBinContent(100,100) = " <<    p1new->GetBinContent(100,100)  << endl;
 cout  << "p1new->GetBinContent(100,101) = " <<    p1new->GetBinContent(100,101)  << endl;
 cout  << "p1new->GetBinContent(101,100) = " <<    p1new->GetBinContent(101,100)  << endl;
 cout  << "p1new->GetBinContent(101,101) = " <<    p1new->GetBinContent(101,101) << endl  << endl;

 Double_t Summe_p = p1new->GetBinContent(100,100)
   + p1new->GetBinContent(100,101)
   + p1new->GetBinContent(101,100)
   + p1new->GetBinContent(101,101);

 cout << "Summe_p = " << Summe_p << "\n" << endl;
 cout << "Summe_p/4. = " << Summe_p/4. << "\n" << endl;

 
 cout << "Summe_p - alpha*Summe_hcorr = " << Summe_p - alpha*Summe_hcorr << endl;
 cout << "same/4. = " << (Summe_p - alpha*Summe_hcorr)/4. << endl;
 cout << "alpha = " << alpha << endl << endl;


 // Upper limit calculation:

 double conf_level = 95.;
 Utilities::TBoundedGaussian fc((conf_level/100.)); // (95%) confidence level
  double upper(0),lower(0);
  
  double Non = Summe_p/4.;
  double Noff = Summe_hcorr/4.;
  double Nexcess;

  std::cout << "Measured Non is: " << Non << std::endl;
  std::cout << "Measured Noff is: " << Noff << std::endl;
  std::cout << "alpha is: " << alpha << std::endl;

  Nexcess = (float)Non - (alpha*(float)Noff);
  //  Nexcess = 3.14;

  std::cout << "Nexcess is Non-alpha*Noff = " << Nexcess << endl << endl;

  fc.GetConfidenceInterval(Non,Noff,alpha,lower,upper);
  if( lower<=0 )
    std::cout << conf_level << "% CL upper limit on Nevents: " << upper << std::endl;
  // else if (upper < Non)
  else {
    std::cout << "Confidence interval for true Nexcess: " << lower << " .. " << upper << std::endl;
    std::cout << "Nexcess is (at " << conf_level << "% confidence level): " << std::endl;
    std::cout << Nexcess << " -" << Nexcess-lower << " +" << upper-Nexcess << std::endl;
  }

   

//   cout << "nHadrons = " << nHadrons << endl;
//   cout << "nGammas = " << nGammas << endl;

  

  // TCanvas *surf = new TCanvas("surfopt","Hallo Welt",100,200,1000,400);
  // surf->Divide(3,1);
  // surf->cd(1);
  // //   h1new->DrawCopy("surf3");
  // //  gPad->SetLogz();
  // //  h1new->Draw("surf3"); 
  // h1new->Draw("COLZ"); 
  // // h1new->Draw(); 
  //  // h1new2->Draw(); 
  //  // h1new2->DrawCopy("surf3"); 
  //  // h1new->DrawCopy("surf3"); 

  // surf->cd(2);
  // p1new->Draw("COLZ"); 
  //  p1new->Draw("surf3"); 
  //  p1new->Draw(); 

  // surf->cd(3);
  // h1->Draw(""); 


//    fo->ls();
   fo->Write();
   //   fo->ls();
   fo->Close();


   delete h1corrnew;
   delete p1new;

   delete s1new;

   //  cout << "Und hilfts denn was, wenn ich das Histo in einer Datei speichere?! Milton fragen?!" << endl << endl;



// chb: This is the END!
}


//   cout << "class name = " << h1new->ClassName() << endl;

//   TCanvas *surf = new TCanvas("surfopt","Hallo Welt",200,200,800,600);
//   surf->Divide(2,1);
// //    surf.SetFillColor(cancolor);
//   surf->cd(1);
//   gStyle->SetPalette(60,0);
//   //  h1new->DrawCopy("surf3");
//   //  h1new->Draw("surf3"); 
//   //  h1new->Draw(); 
//   //  h1new2->Draw(); 
//   //  h1new2->DrawCopy("surf3"); 
//   h1new2->Draw("surf3"); 

//   surf->cd(2);
//   //  (static_cast<TH2D*>(dir->Get("MRSkymapPlotter4_skymap_h1")))->DrawCopy("COLZ");

//   //  h1new->DefineColorLevels(99);
  
//   surf->Update();

//   fo.Write();

  //  surf->Save();
  //  surf->Modified();

  //  cout << "Und hilfts denn was, wenn ich das Histo in einer Datei speichere?! Milton fragen?!" << endl << endl;

//    TPaveText *pt = new TPaveText(0.6,0.85,0.98,0.98,"brNDC");
//    pt->SetFillColor(18);
//    pt->SetTextAlign(12);
//    pt->AddText("Use the axis Context Menu LabelsOption");
//    pt->AddText(" \"a\"   to sort by alphabetic order");
//    pt->AddText(" \">\"   to sort by decreasing values");
//    pt->AddText(" \"<\"   to sort by increasing values");
//    pt->Draw();
    
//  }

//  h1new->DrawCopy();


//  TH2D *hnew = (TH2D*)f.Get(MRSkymapPlotter4_skymap_h1);

//   TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(h1->GetName());
//   if (key2) {
//     TH1 *h2 = (TH1*)key2->ReadObj();
//     h1->Add( h2 );
//     delete h2;
//   }

//  h1new->delete();


//  delete h1new;

//  return;


//  fo.Draw();

//}
    

// void Reserve()
// {
//    TString dir = gSystem->UnixPathName(TCint::GetCurrentMacroName());
//    dir.ReplaceAll("h1draw.C","../hsimple.C");
//    dir.ReplaceAll("/./","/");
//    if (gBenchmark->GetBench("hsimple") < 0) gInterpreter->LoadMacro(dir.Data());
//    TFile *example = (TFile*)gROOT->ProcessLineFast("hsimple(1)");
//    if (!example) return;
   
//    example->ls();
//    TH1 *hpx = (TH1*)example->Get("hpx");
   
//    TCanvas *c1 = new TCanvas("c1","Histogram Drawing Options",200,10,700,900);
//    TPad *pad1 = new TPad("pad1","The pad with the function",0.03,0.62,0.50,0.92,21);
//    TPad *pad2 = new TPad("pad2","The pad with the histogram",0.51,0.62,0.98,0.92,21);
//    TPad *pad3 = new TPad("pad3","The pad with the histogram",0.03,0.02,0.97,0.57,21);
//    pad1->Draw();
//    pad2->Draw();
//    pad3->Draw();

//    // Draw a global picture title
//    TPaveLabel *title = new TPaveLabel(0.1,0.94,0.9,0.98,
//                     "Drawing options for one dimensional histograms");
//    title->SetFillColor(16);
//    title->SetTextFont(52);
//    title->Draw();

//    // Draw histogram hpx in first pad with the default option.
//    pad1->cd();
//    pad1->GetFrame()->SetFillColor(18);
//    hpx->SetFillColor(45);
//    hpx->DrawCopy();
//    TPaveLabel *label1 = new TPaveLabel(-3.5,700,-1,800,"Default option");
//    label1->SetFillColor(42);
//    label1->Draw();

//    // Draw hpx as a lego. Clicking on the lego area will show
//    // a "transparent cube" to guide you rotating the lego in real time.
//    pad2->cd();
//    hpx->DrawCopy("lego1");
//    TPaveLabel *label2 = new TPaveLabel(-0.72,0.74,-0.22,0.88,"option Lego1");
//    label2->SetFillColor(42);
//    label2->Draw();
//    TPaveLabel *label2a = new TPaveLabel(-0.93,-1.08,0.25,-0.92,"Click on lego to rotate");
//    label2a->SetFillColor(42);
//    label2a->Draw();

//    // Draw hpx with its errors and a marker.
//    pad3->cd();
//    pad3->SetGridx();
//    pad3->SetGridy();
//    pad3->GetFrame()->SetFillColor(18);
//    hpx->SetMarkerStyle(21);
//    hpx->Draw("e1p");
//    TPaveLabel *label3 = new TPaveLabel(2,600,3.5,650,"option e1p");
//    label3->SetFillColor(42);
//    label3->Draw();

//    // The following illustrates how to add comments using a PaveText.
//    // Attributes of text/lines/boxes added to a PaveText can be modified.
//    // The AddText function returns a pointer to the added object.
//    TPaveText *pave = new TPaveText(-3.78,500,-1.2,750);
//    pave->SetFillColor(42);
//    TText *t1=pave->AddText("You can move");
//    t1->SetTextColor(4);
//    t1->SetTextSize(0.05);
//    pave->AddText("Title and Stats pads");
//    pave->AddText("X and Y axis");
//    pave->AddText("You can modify bin contents");
//    pave->Draw();
//    c1->Update();
// }
