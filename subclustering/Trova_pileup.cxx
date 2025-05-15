//CALL: ./progname.exe <path to file>

#include <iostream>
#include <unistd.h>
#include <limits.h>
#include <sstream>
#include <fstream>
#include <chrono>
#include <string>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <stdio.h>
#include <cstdlib>
#include <filesystem>
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "Analyzer.h"

using namespace std;

void ScIndicesElem(int nSc, UInt_t npix, float* sc_redpixID, int &nSc_red, vector<int>& B, vector<int>& E);

int main(int argc, char** argv)
{
    if(argc<2) {cerr<<"No file path!\nSuggested use: ./progname.exe <path to file"; exit(EXIT_FAILURE);}
    string nome=argv[1];

    TFile* f = TFile::Open(Form("%s",nome.c_str()),"READ");
    TTree* tree = (TTree*)f->Get("Events");

    //define variables required
    unsigned int nSc=0;
    int nSc_red=0;
    UInt_t Nredpix=0;
    int run;
    int event;

    int nmax=500000;
    int nscmax=100;
    vector<float> sc_redpixID;
    sc_redpixID.reserve(nmax);
    vector<int> XPix;
    XPix.reserve(nmax);
    vector<int> YPix;
    YPix.reserve(nmax);
    vector<float> ZPix;
    ZPix.reserve(nmax);
    vector<float> xmean;
    xmean.reserve(nscmax);
    vector<float> ymean;
    ymean.reserve(nscmax);
    vector<float> width;
    width.reserve(nscmax);
    vector<float> length;
    length.reserve(nscmax);
    vector<float> integral;
    integral.reserve(nscmax);
    vector<float> size;
    size.reserve(nscmax);
    vector<float> nhits;
    nhits.reserve(nscmax);

    //Link the variables to the tree branches
    tree->SetBranchAddress("run",&run);
    tree->SetBranchAddress("event",&event); 
    /////////////Reco variables//////////////     
    tree->SetBranchAddress("nSc",&nSc);
    tree->SetBranchAddress("sc_redpixIdx",sc_redpixID.data());
    tree->SetBranchAddress("nRedpix",&Nredpix);
    tree->SetBranchAddress("redpix_ix",XPix.data());
    tree->SetBranchAddress("redpix_iy",YPix.data());
    tree->SetBranchAddress("redpix_iz",ZPix.data());
    tree->SetBranchAddress("sc_width",width.data());
    tree->SetBranchAddress("sc_length",length.data());
    tree->SetBranchAddress("sc_integral",integral.data());
    tree->SetBranchAddress("sc_size",size.data());
    tree->SetBranchAddress("sc_nhits",nhits.data());
    tree->SetBranchAddress("sc_xmean",xmean.data());
    tree->SetBranchAddress("sc_ymean",ymean.data());

    /////////////////////////////////Analysis Variables ////////////////////////////////////////////////
    vector<int> BeginScPix;
    vector<int> EndScPix;

    /////////////////////////////////Histograms/////////////////////////////////////////////////////////
    int nbinsx=2304;
    int nbinsy=2304;
    TH1D* skewl = new TH1D("skewl","skewl",1000,-2,2);
    TH1D* skewt = new TH1D("skewt","skewt",1000,-2,2);
    TH1D* profl = nullptr;
    TH1D* proft = nullptr;
    //TH2F* hmap_intensity = new TH2F("hmap_intensity","hmap_intensity",nbinsx,0,nbinsx,nbinsy,0,nbinsy);
    //TH2F* hmap_5_9keV = new TH2F("hmap_5_9keV","hmap_5_9keV",nbinsx,0,nbinsx,nbinsy,0,nbinsy);

    
    TFile *fout = TFile::Open(Form("map.root"),"RECREATE");
    int counter=0;

    //Cycle on events (k is the image index)
    cout<<"this run has "<<tree->GetEntries()<<" entries"<<endl;
    for(int k=0;k<tree->GetEntries();k++)
    {
        
        tree->GetEntry(k);
        if(k%500==0) {cout<<"getting entries..."<<endl; cout << "Nev: "<< k << "\nnSc:  " << nSc <<endl;}
        //for reduced pixels:
        ScIndicesElem(nSc,Nredpix,sc_redpixID.data(),nSc_red,BeginScPix,EndScPix);

        for(int clu=0;clu<nSc;clu++) //(this is the cluster index in the kth image)
        {
            if(k%500==0 && clu%20==0) cout<<"Cluster "<< clu << " integral "<<integral[clu]<<endl;         

            if(sc_redpixID[clu]!=-1)
            {
                bool cut1= (integral[clu]>20000);       //remove below 0.25 eV (clear noise) and above 300 keV (no need for alphas). In GIN a muon travelling 150 cm leaves about 300 keV
                bool cut2= (width[clu]/length[clu] < 0.7 && width[clu]<70 && length[clu]>1200);  //Remove noisy hotspot side of GEM
                //bool cut3 = (run>22800);
                if(cut1 && cut2)
                {

                    Analyzer Traccia(Form("Track%i_event%i_run%i",counter,k,run),XPix.data(),YPix.data(),ZPix.data(),BeginScPix[clu],EndScPix[clu]);
                    counter++;
                    proft=Traccia.FillProfile(0);
                    profl=Traccia.FillProfile(1);


                    TH2F* salva = Traccia.GetHistoTrack();
                    salva->SetName(Form("Track_%d_%d_%d_delta%d_size%d",k,clu,static_cast<int>(abs(profl->GetSkewness())*100),static_cast<int>(integral[clu]/nhits[clu]*100),static_cast<int>(nhits[clu]/size[clu]*100)));
                    proft->SetName(Form("ProftTrack_%d_%d",k,clu));
                    profl->SetName(Form("ProflTrack_%d_%d",k,clu));
                    skewl->Fill(abs(profl->GetSkewness()));
                    skewt->Fill(proft->GetSkewness());
                    cout<<"Saved\n";
                    cout<< "SkewT "<<proft->GetSkewness()<<endl;
                    cout<< "SkewL "<<profl->GetSkewness()<<endl;
                    cout<< "nhits/size "<<nhits[clu]/size[clu]<<endl;
                    cout<< "delta "<<integral[clu]/nhits[clu]<<endl;
                    salva->Rebin2D(4,4);
                    proft->Write();
                    profl->Write();
                    salva->Write();
                    //delete salva;
                    delete profl;
                    delete proft;
                    
                }
            }
            
            
        }
    
    }
    skewl->Write();
    skewt->Write();


//    hmap_occu->Write();
    fout->Flush();
    fout->Close();
    f->Close();

    return 0;
}

//Functions
void ScIndicesElem(int nSc, UInt_t npix, float* sc_redpixID, int &nSc_red, vector<int>& B, vector<int>& E)
{
  B.clear();
  E.clear();

  vector<float> sc_redpix_start;

  int parcount=0;

  for(int i=0; i<nSc; i++){
    if(sc_redpixID[i]>=0)  sc_redpix_start.push_back(sc_redpixID[i]);
  }

  nSc_red = sc_redpix_start.size();

  sc_redpix_start.push_back(npix);

  for(int i=0;i<sc_redpix_start.size()-1;i++){
    B.push_back(sc_redpix_start[i]);
    E.push_back(sc_redpix_start[i+1]);
    //std::cout<<B[i]<<" "<<E[i]<<endl;
  }

  sc_redpix_start.clear();

  return;

}