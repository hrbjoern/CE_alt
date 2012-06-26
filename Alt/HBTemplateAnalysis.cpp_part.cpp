  gDirectory->pwd();
  gROOT->cd();
  gDirectory->pwd();
  h1new = static_cast<TH2D*> ( (static_cast<TH2D*>(dir->Get("MRSkymapPlotter4_skymap_h1")))->Clone()); 
  h1new->SetName("h1new");
//  h1new = dynamic_cast<TH2D*> ( (dynamic_cast<TH2D*>(dir->Get("MRSkymapPlotter4_skymap_h1")))->Clone("h1new")); 

  h1new->Print();

  cout << "class name MR = " << static_cast<TH2D*>(dir->Get("MRSkymapPlotter4_skymap_h1"))->ClassName() << endl;

  cout << "class name 2 = " << h1new->ClassName() << endl;

  cout << "Test 0" << endl << endl;

  fi->Close();
  gDirectory->pwd();
  h1new->Print();
  cout << "Test 1" << endl << endl;

  cout << "class name 3 = " << h1new->ClassName() << endl;