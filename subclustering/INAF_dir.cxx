//This in particular compile using  g++ Analyzer.cxx Base_script.cxx -o nameprog `root-config --libs --cflags` -lSpectrum
//Then use as ./nameprog path_to_rootfile

#include <iostream>
#include <string>
#include <vector>
#include "Analyzer.h"
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <chrono>
#include <cmath>  // for M_PI

using namespace std;


// Function to calculate indices for clusters based on reduced pixels
// nSc: number of superclusters
// npix: total number of pixels
// sc_redpixID: array of indices for the first reduced pixel in each cluster
// nSc_red: number of processed superclusters (output)
// B: begin indices for reduced pixels (output)
// E: end indices for reduced pixels (output)
void ScIndicesElem(int nSc, UInt_t npix, float* sc_redpixID, int &nSc_red, vector<int>& B, vector<int>& E){
  nSc_red=0;
  B.clear();
  E.clear();

  vector<float> sc_redpix_start;
  sc_redpix_start.push_back(0); // Initialize the start indices

  int parcount=0; // Counter for iterating through pixels

  // Loop through each supercluster to process its reduced pixels
  for(int i=0; i<nSc; i++){
    if(sc_redpixID[i]>0){
      sc_redpix_start.push_back(sc_redpixID[i]); // Add the start index of the reduced pixels
    }
  }

  nSc_red = sc_redpix_start.size(); // Number of processed superclusters

  sc_redpix_start.push_back(npix); // Ensure the last pixel is included

  // Calculate the begin and end indices for each cluster's reduced pixels
  for(int i=0;i<sc_redpix_start.size()-1;i++){
    B.push_back(sc_redpix_start[i]);
    E.push_back(sc_redpix_start[i+1]);
  }

  sc_redpix_start.clear(); // Clear the temporary storage
}

// Function to find the divisor of A closest to B
int closest_divisor(int A, float B){
  vector<float>  div_dist;
  vector<int> div_candidate;

  // Find all divisors of A and their distance to B
  for(int i=1; i<=A; i++){
    if(A%i==0){ // Check if i is a divisor of A
      div_dist.push_back(abs((float)i-B));
      div_candidate.push_back(i);
    }
  }

  // Find the divisor closest to B
  auto result = min_element(div_dist.begin(), div_dist.end());
  return div_candidate[result-div_dist.begin()];
}

int main(int argc, char** argv){

// Open the input ROOT file and get the "Events" TTree
  TFile* f = TFile::Open(Form("%s",argv[1]));
  TTree* tree = (TTree*)f->Get("Events");

  system("mkdir -p Tracks/");

  // Variable declarations for event data
  int nmax=2000000000;
  int nscmax=5000;
  int npixel=2304;
  int npixelsmall=250;
  float slimnesslimit=0.6;

  unsigned int nSc=0;
  int nSc_red=0;
  UInt_t Nredpix=0;
  Int_t sc_npix=0;
  int run;
  int event;

  // Reserve space for vectors to store pixel and cluster data
  vector<float> sc_redpixID;
  sc_redpixID.reserve(nmax);
  vector<UInt_t> ScNpixels;
  ScNpixels.reserve(nscmax);
  vector<int> XPix;
  XPix.reserve(nmax);
  vector<int> YPix;
  YPix.reserve(nmax);
  vector<float> ZPix;
  ZPix.reserve(nmax);
  // Additional vectors for analysis
  vector<float> xmean;
  xmean.reserve(nscmax);
  vector<float> ymean;
  ymean.reserve(nscmax);
  vector<float> ymin;
  ymin.reserve(nscmax);
  vector<float> ymax;
  ymax.reserve(nscmax);
  vector<float> xmin;
  xmin.reserve(nscmax);
  vector<float> xmax;
  xmax.reserve(nscmax);
  vector<float> scsize;
  scsize.reserve(nscmax);
  vector<float> scnhits;
  scnhits.reserve(nscmax);
  vector<float> v_sc_rms;
  v_sc_rms.reserve(nscmax);
  vector<float> v_sc_tgausssigma;
  v_sc_tgausssigma.reserve(nscmax);
  vector<float> v_sc_theta;
  v_sc_theta.reserve(nscmax);
  vector<float> width;
  width.reserve(nscmax);
  vector<float> length;
  length.reserve(nscmax);
  vector<float> integral;
  integral.reserve(nscmax);

// Prepare for the analysis by setting branch addresses for tree
  // This connects the variables defined above with the corresponding data in the ROOT file
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("event",&event);
  tree->SetBranchAddress("nSc",&nSc);
  tree->SetBranchAddress("sc_redpixIdx",sc_redpixID.data());
  //tree->SetBranchAddress("sc_size",&sc_npix);
  tree->SetBranchAddress("nRedpix",&Nredpix);
  tree->SetBranchAddress("redpix_ix",XPix.data());
  tree->SetBranchAddress("redpix_iy",YPix.data());
  tree->SetBranchAddress("redpix_iz",ZPix.data());
  tree->SetBranchAddress("sc_width",width.data());
  tree->SetBranchAddress("sc_length",length.data());
  tree->SetBranchAddress("sc_integral",integral.data());
  tree->SetBranchAddress("sc_xmean",xmean.data());
  tree->SetBranchAddress("sc_ymean",ymean.data());
  tree->SetBranchAddress("sc_xmin",xmin.data());
  tree->SetBranchAddress("sc_ymin",ymin.data());
  tree->SetBranchAddress("sc_xmax",xmax.data());
  tree->SetBranchAddress("sc_ymax",ymax.data());
  tree->SetBranchAddress("sc_size",scsize.data());
  tree->SetBranchAddress("sc_nhits",scnhits.data());
  tree->SetBranchAddress("sc_theta",v_sc_theta.data());
  tree->SetBranchAddress("sc_rms",v_sc_rms.data());
  tree->SetBranchAddress("sc_tgausssigma",v_sc_tgausssigma.data());

  // Prepare vectors for storing analysis results
  vector<int> BeginScPix;
  vector<int> EndScPix;
  vector<int> BeginScallPix;
  vector<int> EndScallPix;

  // Variables for storing intermediate analysis results
  double phi_DIR=0, phi_DIR_deg=0, phi_PCA=0, phi=0;
  double phi_HT_skew=0, phi_HT_maxpx=0, phi_HT_maxpxRebin=0, phi_HT_peakprof=0, phi_HT_maxprof=0, phi_HT_integral=0;
  double scint = 0;
  double skewness_L=0, skewness_T=0, kurtosis_L=0, kurtosis_T=0;
  double xbar=0, ybar=0;
  double xl=0, yl=0, xr=0, yr=0;
  double reco_theta=0, x_mean=0, y_mean=0, x_min=0, x_max=0, y_min=0, y_max=0;
  int procTime=0;
  double profileMean=0, profileSkew=0; 
  bool puCandidate=false;

  double time=0.;

  //impact point and directionality
  Int_t NPIP=905;     //105 LIME    300 da fiorina
  Float_t wFac=2.5;     //1 da fiorina
  double xIP=0, yIP=0;

  vector<std::pair<double,double>> peakslong;
  vector<std::pair<double,double>> peakstrans;
  int npeaks=0;
  double peak_density=0, tracklength=0, track_width=0, recolength=0, recowidth=0, reco_sc_rms=0, reco_sc_tgausssigma=0;
  vector<double> peak_distance, peak_width;
  

  //ANALYSIS
  int counter=0;
  int ired=0;
  //int counterall=0;

  // Prepare output file path and open output file
  // Assuming 'argv[1]' is the input file and 'argv[2]' is the output directory
  string inputFilePath = argv[1];
  string outputDir = argc > 2 ? argv[2] : "."; // Use current directory if not provided
  string filename = inputFilePath.substr(inputFilePath.find_last_of("/\\") + 1);
  string outputFilePath = outputDir + "/Analysis_" + filename;

  //string filename(argv[1]);
  //filename = filename.substr(filename.find_last_of("/\\")+1);
  cout<<filename<<endl;
  TFile* fout = new TFile(outputFilePath.c_str(), "recreate");
  cout<<Form("%s/Analysis_%s",argv[2],filename.c_str())<<endl;
  fout->cd();
  fout->mkdir("Tracks");

  TTree* out_tree = new TTree("track_info","track_info");
  out_tree->Branch("run",&run);
  out_tree->Branch("event",&event);
  out_tree->Branch("nSc",&nSc);
  out_tree->Branch("nSc_red",&nSc_red);
  out_tree->Branch("Integral",&scint);
  out_tree->Branch("ScSize",&sc_npix);
  out_tree->Branch("AngleDIR",&phi_DIR);
  out_tree->Branch("DegreeDIR",&phi_DIR_deg);
  out_tree->Branch("AnglePCA",&phi_PCA);
  out_tree->Branch("RecoTheta",&reco_theta);
  out_tree->Branch("RecoScRMS",&reco_sc_rms);
  out_tree->Branch("RecoScTGaussSigma",&reco_sc_tgausssigma);
  out_tree->Branch("X_ImpactPoint",&xIP);
  out_tree->Branch("Y_ImpactPoint",&yIP);
  out_tree->Branch("Xmean",&x_mean);
  out_tree->Branch("Ymean",&y_mean);
  out_tree->Branch("Xmin",&x_min);
  out_tree->Branch("Ymin",&y_min);
  out_tree->Branch("Xmax",&x_max);
  out_tree->Branch("Ymax",&y_max);
  out_tree->Branch("XBar",&xbar);
  out_tree->Branch("YBar",&ybar);
  out_tree->Branch("procTime",&procTime);
  out_tree->Branch("ProfileMean",&profileMean);
  out_tree->Branch("ProfileSkew",&profileSkew);
  out_tree->Branch("PileUpCandidate",&puCandidate);
  /*
  out_tree->Branch("Npeaks",&npeaks);
  out_tree->Branch("PeakDensity",&peak_density);
  out_tree->Branch("TrackLength",&tracklength);
  out_tree->Branch("TrackWidth",&track_width);
  out_tree->Branch("RecoLength",&recolength);
  out_tree->Branch("RecoWidth",&recowidth);
  out_tree->Branch("PeakDistance",&peak_distance);
  out_tree->Branch("Skewness",&skewness_L);
  out_tree->Branch("Kurtosis",&kurtosis_T);
  out_tree->Branch("phi_HT_skew",&phi_HT_skew);
  out_tree->Branch("phi_HT_maxpx",&phi_HT_maxpx);
  out_tree->Branch("phi_HT_maxpxRebin",&phi_HT_maxpxRebin);
  out_tree->Branch("phi_HT_peakprof",&phi_HT_peakprof);
  out_tree->Branch("phi_HT_maxprof",&phi_HT_maxprof);
  out_tree->Branch("phi_HT_integral",&phi_HT_integral);
  */

  int pileUpCounter=0;
  int totalentries=tree->GetEntries();
  cout<<"this run has "<<totalentries<<" entries"<<endl;
  auto t_tot0 = std::chrono::steady_clock::now();
  for(int k=0;k<totalentries;k++)
  //for(int k=0;k<1;k++)
  {
    if (k%100==0) cout<<"Entry "<<k<<" Percentage "<<static_cast<float>(k)/totalentries*100<<"%"<<endl;
    sc_redpixID.clear();
    tree->GetEntry(k);
    //cout << "Nev: "<< k << "\nnSc:  " << nSc << " event "<< event <<endl;
    //for reduced pixels:
    ScIndicesElem(nSc,Nredpix,sc_redpixID.data(),nSc_red,BeginScPix,EndScPix);

    //cout<<"nSc "<<nSc<<" nSc_red "<<nSc_red<<" Nredpix "<<Nredpix<<endl;

    //Start the cycle on the supercluster of the event
    int pixcounter =0;
    for(int i=0;i<nSc_red;i++)
    {
      scint = integral[i];
      recolength=length[i];
      recowidth=width[i];
      sc_npix = scnhits[i];
      reco_theta=v_sc_theta[i];
      x_mean=xmean[i];
      y_mean=ymean[i];
      x_min=xmin[i];
      x_max=xmax[i];
      y_min=ymin[i];
      y_max=ymax[i];
      reco_sc_rms=v_sc_rms[i];
      reco_sc_tgausssigma=v_sc_tgausssigma[i];
      //cout<<"SC "<<i<<endl;

      pixcounter += EndScPix[i] - BeginScPix[i];
      //cout<<"counted pix: "<<pixcounter<<endl;

      //! Condition to filter out certain events based on physical properties
      //For Polarized 8Kev photon in MANGO
      //if(scint>25000 && scint<50000 && recowidth/recolength>0.7 && recowidth/recolength<1 && x_mean>900 && x_mean<1350 && y_mean<1350 && y_mean>900 && run>22700)
      if(scint>25000 && scint<50000 && recowidth/recolength>0.7 && recowidth/recolength<1 && x_mean>1000 && x_mean<1100 && y_mean<1100 && y_mean>1000 && run>22700)
      //For Polarized 17Kev photon in MANGO
      //if( x_mean>900 && x_mean<1350 && y_mean<1350 && y_mean>900 && scint<90000 && scint>60000 && sc_npix<6000  && run>22700)
      // For LIME 55Fe
      //if (y_max < 1250 && y_min > 1050 && x_max < 1250 && x_min > 1050 && scint>2000 && reco_sc_rms>5 && reco_sc_tgausssigma>2.63 && reco_sc_tgausssigma<4.5 && recowidth/recolength>0.6 )
      {
        // Start timers for the entire track processing
        

      bool cut1 = (scint>26000 && scint<46200);      //For 8.7 keV Fusion 26000 to 46200 (2 sigma region)
      bool cut2 = (x_mean>800 && x_mean<1400);
      bool cut3 = (y_mean>800 && y_mean<1400);
      bool cut4 = (run>22751 && run<22800);
      bool cut5 = (recolength<250);   //Per i 8.7keV

      //if(y_max>1600 && y_max<2000 &&scint>0 && recowidth/recolength<0.4)//for HeCF4 60-40 Fusion
      if(cut1 && cut2 && cut3 && cut4 && cut5)
      {
        auto t0 = std::chrono::steady_clock::now();

        Analyzer Traccia(Form("Track%i_event%i_run%i",counter,k,run),XPix.data(),YPix.data(),ZPix.data(),BeginScPix[i],EndScPix[i]);
        if (Traccia.Getbinmax()<=0) continue; //skip if the track is empty
        //Traccia.SavePic(Form("Track%i_event%i_%i.png",counter,event,run));
        Traccia.SetWScal(wFac);
        Traccia.SetNPIP(NPIP);//Threshold in number of selected points during IP search
        Traccia.ApplyThr();   //10 da fiorina
        Traccia.RemoveNoise(0);   //30 da fiorina
        //puCandidate = Traccia.PileUpCandidate(false, counter, true, 0.2, 2.0); //for 17keV
        puCandidate = Traccia.PileUpCandidate(false, counter, false, 0.,0.); //for 8keV
        Traccia.GetTrackProfileStats(profileMean, profileSkew);
        Traccia.ImpactPoint(Form("TrackIPRegion%i_run%i_evt%i",k,run,counter));
        Traccia.ScaledTrack(Form("TrackScaled%i_run%i_evt%i",k,run,counter));
        Traccia.Direction();
        Traccia.ImprCorrectAngle();
        Traccia.BuildLineDirection();



        // Grab results
        xIP         = Traccia.GetXIP();
        yIP         = Traccia.GetYIP();
        phi_DIR     = Traccia.GetDir();
        phi_DIR_deg = phi_DIR * (180.0 / M_PI);
        phi_PCA     = Traccia.AngleLineMaxRMS();
        xbar        = Traccia.GetXbar();
        ybar        = Traccia.GetYbar();

        // End timing for the entire block
        auto tEnd = std::chrono::high_resolution_clock::now();

        if (puCandidate) {
          pileUpCounter++;
          //Traccia.TrackProfilePlotSave(Form("Track%i_event%i_run%i_entry%i.png", counter, event, run,k));
          continue;
        }

        //! Print only every N events
        if (k % 1000 == 0)
        {
            // Save a diagnostic image
            if(k%500==0 && i%5==0) Traccia.SavePicDir(Form("Track%i_event%i_run%i.png",counter,event,run));
            if (fabs(phi_DIR_deg) < 25. || fabs(phi_DIR_deg) > 150.){
            Traccia.SavePicDir(Form("Track%i_event%i_run%i_entry%i.png", counter, event, run,k));
            }

            std::cout << "Processing entry " << k << " / " << tree->GetEntries() << std::endl;
            std::cout << "counter: " << counter << std::endl;
            std::cout << "XIP: " << xIP << "  YIP: " << yIP << std::endl;
            std::cout << "XIPPrev: " << Traccia.GetXIPPrev() << "  YIPPrev: " << Traccia.GetYIPPrev() << std::endl;
            std::cout << "Degree: " << phi_DIR_deg << "  tan(angle): " << std::tan(phi_DIR) << std::endl;

            // Compute durations in microseconds

            // Print each as "us" and as "% of total"
            auto printTime = [&](const char* label, long long dt) {
                double pct = 100.0 * (double)dt / (double)dt_total;
                std::cout << label << dt << " us (" << pct << " %)\n";
            };

            printTime("Time Analyzer ctor:         ", dt_ctor);
            printTime("Time ApplyThr:              ", dt_applyThr);
            printTime("Time RemoveNoise:           ", dt_removeNoise);
            printTime("Time PileUpCandidate:       ", dt_pileUpCandidate);
            printTime("Time ImpactPoint:           ", dt_impactPoint);
            printTime("Time ScaledTrack:           ", dt_scaledTrack);
            printTime("Time Direction:             ", dt_direction);
            printTime("Time BuildLineDirection:    ", dt_buildLineDirection);

            // Lastly, the total is always 100%
            std::cout << "TOTAL time for all steps: " << dt_total << " us (100%)\n\n";
        }

        // Fill the tree with the total time for this event
        procTime = std::chrono::duration_cast<std::chrono::microseconds>(tEnd - t0).count();
        out_tree->Fill();
        

      }
    counter++;
    }//superclusters

  sc_redpixID.resize(nSc);//Madonna fai il resize altriemnti lui si ricorda la dimensione precedente

  }//ttree entries
  auto t_tot1 = std::chrono::steady_clock::now();
  std::chrono::duration<double> durtmp2=t_tot1-t_tot0;
  double timetot=durtmp2.count();
  cout << "Total time taken in seconds is: " << timetot << "\n"<<endl;

  out_tree->Write();
  fout->Close();
  cout<<"pile up percentage: "<<(double)pileUpCounter/(double)counter<<endl;

  return 0;
}