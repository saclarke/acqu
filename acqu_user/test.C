{
gStyle->SetOptStat(1); 
gStyle->SetPalette(1) ;
TChain *h1 = new TChain("h1") ;

// TFile *f1=new TFile("Theta_generated_Alpha_12_beta_1_9_E_110_sigma_fX_130_300000.root", "recreate", "w");
// h1-> Add("/home/sokhoyan/acqu/Alpha_12_beta_1_9_E_130_sigma_fY_130.root") ;

//File *f1=new TFile("Theta_generated_Alpha_10_8_beta_4_E_130_sigma_fY_130_300000.root", "recreate", "w");

//   h1-> Add("/home/sokhoyan/acqu/Compton_145_150_MeV.root") ;
 h1-> Add("/home/scott/acqu/acqu_user/Dilepton_firstStart_200_300_MeV.root") ;

float M_outgoing_Proton ;

float Px_bm, Py_bm, Pz_bm ;
float test1, Pt_bm, En_bm, En_l0101,  En_l0101_unit, Px_l0101,  Py_l0101, Pz_l0101 ;
float En_l0203,  En_l0203_unit, Px_l0203,  Py_l0203, Pz_l0203 ;
float En_l0302,  En_l0302_unit, Px_l0302,  Py_l0302, Pz_l0302 ;
float En_l0302,  En_l0302_unit, Px_l0302,  Py_l0302, Pz_l0302 ;

float En_l0114,  En_l0114, Px_l0114,  Py_l0114, Pz_l0114 ;

 h1-> SetBranchAddress("Px_bm", &Px_bm) ;
 h1-> SetBranchStatus("Px_bm", 1) ;
 h1-> SetBranchAddress("Py_bm", &Py_bm) ;
 h1-> SetBranchStatus("Py_bm", 1) ;
 h1-> SetBranchAddress("Pz_bm", &Pz_bm) ;
 h1-> SetBranchStatus("Pz_bm", 1) ;
 h1-> SetBranchAddress("En_bm", &En_bm) ;
 h1-> SetBranchStatus("En_bm", 1) ;

 h1-> SetBranchAddress("En_l0203", &En_l0203) ;
 h1-> SetBranchStatus("En_l0203", 1) ;
 h1-> SetBranchAddress("Px_l0203", &Px_l0203) ;
 h1-> SetBranchStatus("Px_l0203", 1) ;
 h1-> SetBranchAddress("Py_l0203", &Py_l0203) ;
 h1-> SetBranchStatus("Py_l0203", 1) ;
 h1-> SetBranchAddress("Pz_l0203", &Pz_l0203) ;
 h1-> SetBranchStatus("Pz_l0203", 1) ;
 
 h1-> SetBranchAddress("En_l0302", &En_l0302) ;
 h1-> SetBranchStatus("En_l0302", 1) ;
 h1-> SetBranchAddress("Px_l0302", &Px_l0302) ;
 h1-> SetBranchStatus("Px_l0302", 1) ;
 h1-> SetBranchAddress("Py_l0302", &Py_l0302) ;
 h1-> SetBranchStatus("Py_l0302", 1) ;
 h1-> SetBranchAddress("Pz_l0302", &Pz_l0302) ;
 h1-> SetBranchStatus("Pz_l0302", 1) ;
 
 h1-> SetBranchAddress("En_l0114", &En_l0114) ;
 h1-> SetBranchStatus("En_l0114", 1) ;
 h1-> SetBranchAddress("Px_l0114", &Px_l0114) ;
 h1-> SetBranchStatus("Px_l0114", 1) ;
 h1-> SetBranchAddress("Py_l0114", &Py_l0114) ;
 h1-> SetBranchStatus("Py_l0114", 1) ;
 h1-> SetBranchAddress("Pz_l0114", &Pz_l0114) ;
 h1-> SetBranchStatus("Pz_l0114", 1) ;

int nentries = h1->GetEntries() ;
TH1F *h_theta_gen = new TH1F("h_theta_gen", "h_theta_gen", 6, 40, 160);
TH1F *h_E_gen = new TH1F("h_E_gen", "h_E_gen", 6, 40., 160);
TH2F *h_E_Theta_gen = new TH2F("h_E_Theta_gen", "h_E_Theta_gen", 1000, 0., 1000, 180, 0, 180);

 TH1F *hincoming_PhotonE = new TH1F("hincoming_PhotonE", "hincoming_PhotonE", 200, 0., 200.);

 //h1->Draw("acos(Pz_l0101/sqrt(Px_l0101**2+ Py_l0101**2 + Pz_l0101**2))*180/TMath::Pi() >> h_theta_gen") ;
 //h_theta_gen->Rebin(20) ;
 //h_theta_gen->SetLineColor(4) ;
 
 //h1->Draw("En_l0101*1000 >> h_E_gen") ;
 
 //h_E_gen->Write() ;
 //f1->Close() ;
 TH1F* hM_invar = new TH1F("hM_invar", "hM_invar", 400, 0., 400.) ;
TH1F* hM_miss = new TH1F("hM_miss", "hM_miss", 400, 0., 2000.) ; 
TH1F* hM_miss_nocut = new TH1F("hM_miss_nocut", "hM_miss_nocut", 400, 0., 2000.) ; 
TH1F* hM_miss_nocut_direct = new TH1F("hM_miss_nocut_direct", "hM_miss_nocut_direct", 400, 0., 2000.) ; 
TH1F* hE_conservation = new TH1F("hE_conservation", "hE_conservation", 1500, 0., 1500.) ; 

TH2F* hTheta_M_miss_nocut = new TH2F("hTheta_M_miss_nocut", "hTheta_M_miss_nocut", 180, 0, 180, 400, 0., 2000.) ; 
TH2F* hEnergy_M_miss_nocut = new TH2F("hEnergy_M_miss_nocut", "hEnergy_M_miss_nocut", 200, 0, 200, 400, 0., 2000.) ; 
TH2F* hEnergy_Theta_nocut_gamma = new TH2F("hEnergy_Theta_nocut_gamma", "hEnergy_Theta_nocut_gamma", 200, 0, 500, 180, 0., 180.) ; 

TH2F* hEnergy_Theta_nocut_elpos = new TH2F("hEnergy_Theta_nocut_elpos", "hEnergy_Theta_nocut_elpos", 200, 0, 500, 180, 0., 180.) ; 
TH2F* hEnergy_Theta_elpos = new TH2F("hEnergy_Theta_elpos", "hEnergy_Theta_elpos", 200, 0, 500, 180, 0., 180.); 
TH2F* hEnergy_Theta_elpos_proton_thetaCut = new TH2F("hEnergy_Theta_elpos_proton_thetaCut", "hEnergy_Theta_elpos_proton_thetaCut", 200, 0, 500, 180, 0., 180.); 

TH2F* hTheta1_Theta2_elpos = new TH2F("hTheta1_Theta2_elpos", "hTheta1_Theta2_elpos", 180, 0., 180., 180, 0., 180.); 

TH2F* hEnergy_Theta_nocut_proton = new TH2F("hEnergy_Theta_nocut_proton", "hEnergy_Theta_nocut_proton", 200, 0, 2000, 180, 0., 180.) ; 
TH2F* hEnergy_Theta_nocut_proton = new TH2F("hEnergy_Theta_nocut_proton", "hEnergy_Theta_nocut_proton", 200, 0, 200, 180, 0., 180.) ; 

TH1F* hMinv_elpos = new TH1F("hMinv_elpos", "hMinv_elpos", 300, 0., 1000.) ; 

hM_miss_nocut->SetLineColor(4) ;
hM_miss->SetLineColor(2) ;

 TLorentzVector outgoing_Proton_calc ;
 nentries = h1->GetEntries();
  
for (Int_t  i = 0; i < nentries; i ++)
{
if(i%10000 == 0)
cout << i << endl ;

TFile* currentFile = h1->GetFile() ;
// strcpy(sname, currentFile->GetName()) ;

              h1->GetEntry(i) ;
   
   float outgoing_Proton_totP = sqrt(En_l0114**2 - 0.938272**2) ;
   float outgoing_Electron_totP = sqrt(En_l0203**2 - 0.000510998**2) ;
   float outgoing_Positron_totP = sqrt(En_l0302**2 - 0.000510998**2) ;
      
   TLorentzVector incoming_Photon(Px_bm*En_bm*1000., Py_bm*En_bm*1000., Pz_bm*En_bm*1000., En_bm*1000.) ;
   TLorentzVector outgoing_Electron(Px_l0203*outgoing_Electron_totP*1000.,  Py_l0203*outgoing_Electron_totP*1000., Pz_l0203*outgoing_Electron_totP*1000., En_l0203*1000.) ;
   TLorentzVector outgoing_Positron(Px_l0302*outgoing_Positron_totP*1000.,  Py_l0302*outgoing_Positron_totP*1000., Pz_l0302*outgoing_Positron_totP*1000., En_l0302*1000.) ;
   TLorentzVector Target_4vec(0., 0., 0., 938.272) ;  
   TLorentzVector outgoing_Proton(Px_l0114*outgoing_Proton_totP*1000., Py_l0114*outgoing_Proton_totP*1000., Pz_l0114*outgoing_Proton_totP*1000.,  En_l0114*1000.) ;

   hM_miss_nocut_direct->Fill( sqrt(outgoing_Proton.E()**2 - (Px_l0114**2 + Py_l0114**2 + Pz_l0114**2)) ) ;

//  cout << incoming_Photon.Px() << endl ;
// cout << sqrt(outgoing_Electron.E()**2 - (Px_l0114**2 + Py_l0114**2 + Pz_l0114**2)) << endl ;
 //  cout << sqrt(outgoing_Proton.E()**2 - (outgoing_Proton.Px()**2 + outgoing_Proton.Py()**2 + outgoing_Proton.Pz()**2)) << endl ;
// cout << outgoing_Proton.Pz()**2<< endl ;
// cout << outgoing_Proton.Py() + outgoing_Electron.Py() + outgoing_Positron.Py() << endl;
 
   outgoing_Proton_calc = incoming_Photon + Target_4vec - outgoing_Electron - outgoing_Positron  ;
   M_outgoing_Proton = outgoing_Proton_calc.M() ;
 
 //   cout << outgoing_Proton_calc.Phi() << " " << outgoing_Proton.Phi() << endl ;
  // cout << outgoing_Proton.M() << endl ;
  // cout << outgoing_Positron.M() << endl ;
  
   hTheta_M_miss_nocut->Fill(outgoing_Electron.Theta()*180/TMath::Pi(), M_outgoing_Proton, 1) ;
   hTheta_M_miss_nocut->Fill(outgoing_Positron.Theta()*180/TMath::Pi(), M_outgoing_Proton, 1) ;
   
   hEnergy_M_miss_nocut->Fill(outgoing_Electron.E(), M_outgoing_Proton, 1) ;
   hEnergy_M_miss_nocut->Fill(outgoing_Positron.E(), M_outgoing_Proton, 1) ;
    
   hEnergy_Theta_nocut_elpos->Fill(outgoing_Electron.E(), outgoing_Electron.Theta()*180/TMath::Pi(), 1) ;
   hEnergy_Theta_nocut_elpos->Fill(outgoing_Positron.E(), outgoing_Positron.Theta()*180/TMath::Pi(), 1) ;
 //     if(outgoing_Electron.Theta()*180/TMath::Pi() > 20. > outgoing_Positron.Theta() > 20. && outgoing_Proton.Theta() > 20.)
   hTheta1_Theta2_elpos->Fill(outgoing_Electron.Theta()*180/TMath::Pi(), outgoing_Positron.Theta()*180/TMath::Pi(), 1); 
 
  if(outgoing_Electron.Theta()*180/TMath::Pi() > 20. && outgoing_Positron.Theta()*180/TMath::Pi() > 20.) 
    {
   hEnergy_Theta_elpos->Fill(outgoing_Electron.E(), outgoing_Electron.Theta()*180/TMath::Pi(), 1) ;
   hEnergy_Theta_elpos->Fill(outgoing_Positron.E(), outgoing_Positron.Theta()*180/TMath::Pi(), 1) ;
   if(outgoing_Proton.Theta()*180/TMath::Pi() > 20.)
   {
   hEnergy_Theta_elpos_proton_thetaCut->Fill(outgoing_Positron.E(), outgoing_Positron.Theta()*180/TMath::Pi(), 1) ;
   }   
   }
  
  hEnergy_Theta_nocut_elpos->Fill(outgoing_Electron.E(), outgoing_Positron.Theta()*180/TMath::Pi(), 1) ;
  hEnergy_Theta_nocut_elpos->Fill(outgoing_Positron.E(), outgoing_Positron.Theta()*180/TMath::Pi(), 1) ;
  
//  hEnergy_Theta_nocut_proton->Fill(outgoing_Proton.E() - 938.272, outgoing_Proton.Theta()*180/TMath::Pi(), 1) ;
   
//   cout << outgoing_Proton.Theta()*180/TMath::Pi() << " " << outgoing_Proton.E() << endl ;
  
   hEnergy_Theta_nocut_proton->Fill(outgoing_Proton.E() - 938.272, outgoing_Proton.Theta()*180/TMath::Pi(), 1) ;

   hM_miss_nocut->Fill(M_outgoing_Proton) ;
  
   //////////////////////////////
   
   hMinv_elpos->Fill((outgoing_Electron + outgoing_Positron).M()) ;
     
   hincoming_PhotonE->Fill(En_bm*1000) ;
   
 //   h_theta_gen->Fill(En_bm*1000) ;

////  hincoming_PhotonE->Fill(incoming_Photon.E()) ;
// cout << (outgoing_Electron + outgoing_Positron).M() << endl ;
  
}
//hM_miss_nocut_direct->Draw() ;
//hincoming_PhotonE->Draw() ;

TCanvas *check_mass = new TCanvas("check_mass","check_mass",100,100, 400, 350) ;
check_mass->Divide(1, 1) ;
check_mass->cd(1);
hMinv_elpos->Draw("") ;

	#include <TSpectrum.h>
	TSpectrum *s = new TSpectrum();
	//Search function goes into the histogram and counts the number of peaks
	double nfound = s->Search(hMinv_elpos, 1, "nobackground", 0.9);
	//printf("Found %d candidate peaks\n", nfound);

	//The Search function saves the (x,y) co-ords of the peaks, this finds these values and prints it to the screen
	float *xpeaks = s->GetPositionX();
	for(int p=0; p<nfound; p++){
		float xp = xpeaks[p];
		int bin = hMinv_elpos->GetXaxis()->FindBin(xp);
		float yp = hMinv_elpos->GetBinContent(bin);
		cout << "X Position: " <<  xp << "  Y Position: " << yp << endl;
	}


TCanvas *check_theta = new TCanvas("check_theta","check_theta",100,100, 400, 550) ;
check_theta->Divide(2, 3) ;

check_theta->cd(1) ;
hEnergy_Theta_nocut_elpos->Draw("colz") ;
check_theta->cd(2) ;
hEnergy_Theta_elpos->Draw("colz") ;
check_theta->cd(3) ;
hTheta1_Theta2_elpos->Draw("colz") ;
check_theta->cd(4) ;
hEnergy_Theta_elpos_proton_thetaCut->Draw("colz") ;
check_theta->cd(5) ;
hEnergy_Theta_nocut_proton->Draw("colz") ;

TCanvas *check = new TCanvas("check","check",100,100, 400, 350) ;
check->Divide(2, 2) ;

check->cd(1) ;
hTheta_M_miss_nocut->Draw("colz") ;
check->cd(2) ;
hEnergy_M_miss_nocut->Draw("colz") ;
check->cd(3) ;
hEnergy_Theta_nocut_gamma->Draw("colz") ;
TLine *line = new TLine(100., 0., 100., 180.); 
line->SetLineColor(2) ;
line->Draw() ;
// check->cd(5) ;
// hM_miss_nocut->Draw() ;
check->cd(4) ;
hM_miss_nocut->Draw() ;
// shEnergy_Theta_nocut_proton->Draw("colz") ;
//h_E_gen->Write() ;
//f1->Close() ;

}
