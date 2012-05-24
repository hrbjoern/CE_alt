#include <iostream>
#include "TInterpreter.h"
#include "TCanvas.h"
#include <TROOT.h>
#include <TMath.h>
#include <cmath>
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

// using std::cout;
// using std::cin;
// using std::endl;

// HB's variables:
Double_t xDist, yDist, rDist;
Double_t nHadrons=0.;
Double_t nGammas=0.;
//void draw();   

void HBTemplateAnalysis(){

  //Preliminaria:
  gROOT->Reset();
  //  gStyle->SetPalette(10,0);

  // Declare Histos:
  //  static  TH2D *h1new = new TH2D();
  TH2D *h1new = new TH2D();
  TH2D *p1new = new TH2D();
  TH2D *h1 = new TH2D();

  //chb test:
  cout << "class name 1 = " << h1new->ClassName() << endl;

  // Open input ROOT file:
  //  TFile *f = TFile::Open("FornaxClusterWobble_2010-05-27_template1.root");
  TFile *fi = new TFile("FornaxClusterWobble_2010-05-27_template1.root");
  //  fi->ls();
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
  //  h1new = static_cast<TH2D*> ( (static_cast<TH2D*>(dir->Get("MRSkymapPlotter4_skymap_h1")))->Clone("h1new"));  
  h1new = static_cast<TH2D*> ( (static_cast<TH2D*>(dir->Get("MRSkymapPlotter4_skymap_h1_cor")))->Clone("h1new"));  
  cout << "Warning! We're now using the 'corrected' hadron map." << endl << endl;
  // Photon candidates:
  p1new = static_cast<TH2D*> ( (static_cast<TH2D*>(dir->Get("MRSkymapPlotter4_skymap_p1")))->Clone("p1new"));  // 
  // Testing:
  h1 = static_cast<TH2D*> ( (static_cast<TH2D*>(h1new))->Clone("h1"));  // 
 

  //  h1new->SetTitle("h1new");
  //  h1new = dynamic_cast<TH2D*> ( (dynamic_cast<TH2D*>(dir->Get("MRSkymapPlotter4_skymap_h1")))->Clone("h1new")); 

  //  h1new->Print();

  //  cout << "class name MR = " << static_cast<TH2D*>(dir->Get("MRSkymapPlotter4_skymap_h1"))->ClassName() << endl;
  //  cout << "class name 2 = " << h1new->ClassName() << endl;
  //  cout << "Test 0" << endl << endl;

  // Close input ROOT file:
  fi->Close();
  gDirectory->pwd();
  h1new->Print();
  // cout << "Test 1" << endl << endl;

  // chb test:
  //  cout << "class name 3 = " << h1new->ClassName() << endl;
  //  h1new->DrawCopy();

  TH2D *h1new2 = new TH2D("h1new2","h1new2",10,0,10,10,0,10);
//   h1new2 = dynamic_cast<TH2D*> ( (dynamic_cast<TH2D*>(dir->Get("MRSkymapPlotter4_skymap_h1")))->Clone("h1new2"));
  //h1new2->Draw();
  

//   //  if(h1new) {
//   //  h1new->Save();
//   // Filenamen geben?  h1new->Write();



   // TH1F *h1 = new TH1F("h1","h1",100,-3,3);
   // //   TH1F *h2 = new TH1F("h2","h2",100,-3,3);
   // TRandom r;
   // for (Int_t i=0;i<100000;i++) {
   //    Double_t x1 = r.Gaus(-1,0.5);
   //    Double_t x2 = r.Gaus(1,1.5);
   //    if (i <1000){
   // 	h1->Fill(x1);
   // 	h1new2->Fill(x1,x1);
   //    }
   //    //      h2->Fill(x2);
   // }
  //    h1->Draw();


  // Open output ROOT file:
  TFile *fo = new TFile("HBTemplateAnalysis_Out.root","RECREATE");
   
   //   cout << "Test 1" << endl << endl;
   //   fo->Write();
   //   cout << "Test 4" << endl << endl;
   //   fo->ls();

   h1new->Write();
   p1new->Write();
   h1new2->Write();

   //   cout << "Test 5" << endl << endl;

   //   fo->ls();
   //   fo->Close();

   //   cout << "Test 6" << endl << endl;


   //*... Starting Fornax Analysis...................................................................*//
   //
   // Fornax Cluster:
   // SetTargetName:         Fornax Galaxy Cluster
   // SetObjectRA:           3h38m30s  ... x-axis -> 3.64
   // SetObjectDec:          -35d27'00'' ... y-axis -> -35.45



   //   cout  << "h1new->GetBinContent(25,25) = " <<    h1new->GetBinContent(25,25)  << endl << endl;

   cout  << "h1new->GetNbinsX() = " <<    h1new->GetNbinsX()  << endl;
   cout  << "h1new->GetNbinsY() = " <<    h1new->GetNbinsY()  << endl;

   cout  << "h1new->fNcells = " <<    h1new->fNcells << endl;
   cout  << "p1new->fNcells = " <<    p1new->fNcells << endl << endl;


   // for (int i=0; i<(h1new->GetNbinsX());i++) {
   //   cout <<  "h1new->GetBinContent("<<i<<"),25) = " <<    h1new->GetBinContent(i,25)  << endl;
   //   cout <<  "h1new->GetXaxis()->GetBinCenter("<<i<<") = " <<    h1new->GetXaxis()->GetBinCenter(i)  << endl;
   //   cout <<  "h1new->GetYaxis()->GetBinCenter("<<i<<") = " <<    h1new->GetYaxis()->GetBinCenter(i)  << endl
   // 	  << endl;
   // }


   // Calculate distance from nominal center, if smaller than 1 deg, fill event counter:
   for (int i=0; i<(h1new->GetNbinsX());i++) {
     xDist = (h1new->GetXaxis()->GetBinCenter(i))*15. - 54.625; // distance: note coord's!
     xDist = xDist*cos((h1new->GetYaxis()->GetBinCenter(i))*M_PI/180.);
     //     xDist = xDist*cos(35.45*M_PI/180.);

     for (int j=0; j<h1new->GetNbinsY();j++) {
       yDist =  (h1new->GetYaxis()->GetBinCenter(j)) + 35.45;
       rDist =sqrt(xDist*xDist + yDist*yDist);
       
       if (rDist < 1.){
	 //	 cout << "i, j = " << i << ", " <<  j << ", xDist, yDist, rDist = " << xDist << " " << yDist << " " << rDist << endl;
	 // chb test: 
	 h1->Fill((h1new->GetXaxis()->GetBinCenter(i)),(h1new->GetYaxis()->GetBinCenter(j)),2000);

	 nHadrons += h1new->GetBinContent(i,j) ; 
	 nGammas += p1new->GetBinContent(i,j)  ;
	 
       }
       
     }
   }
   

   cout << "nHadrons (corrected) = " << nHadrons << endl;
   cout << "nGammas = " << nGammas  << endl << endl;

   cout << "Using alpha(Template) = 0.122097 (from MR output), this corresponds to" << endl
	<< "nsignal = nGammas - alpha* nHadronsCorr = " << nGammas- 0.122097*nHadrons << endl << endl;

  
   /* Drawing: ..............................................*/

  TCanvas *surf = new TCanvas("surfopt","Hallo Welt",100,200,1000,400);
  surf->Divide(3,1);
  surf->cd(1);
  //   h1new->DrawCopy("surf3");
  //  gPad->SetLogz();
  //  h1new->Draw("surf3"); 
  h1new->Draw("COLZ"); 
  // h1new->Draw(); 
   // h1new2->Draw(); 
   // h1new2->DrawCopy("surf3"); 
   // h1new->DrawCopy("surf3"); 

  surf->cd(2);
  p1new->Draw("COLZ"); 
   p1new->Draw("surf3"); 
   p1new->Draw(); 

  surf->cd(3);
  h1->Draw("COLZ"); 


//    fo->ls();
   fo->Write();
   //   fo->ls();
   fo->Close();


   delete h1new;
   delete p1new;
   delete h1new2;
   delete h1;

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
