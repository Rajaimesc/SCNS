#pragma once 

#include "SCTripleMatch.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"


#include <math.h>
#include <stdio.h>
//using namespace ana;

void SCTripleMatch( ){
 
 std::ifstream file("DataList.txt"); //We take a number of CAF files  from which we take the data of interest 
 std::string filename;
 int files=6;




TH2F* AllBarycentersCRT = new TH2F(" AllBarycentersCRT", " AllBarycentersCRT", 25, 0., 100.,50,0.,200.); //This histogram will fill all the individual histograms for this event
for(int filenum = 0; filenum < files; filenum++) //Loop over all the files that are into the file
{
 file>>filename;
 //filename= getline(file,text[0]);
 const char * c = filename.c_str();
 SpectrumLoader loader(c);


 const Binning kTest = Binning:: Simple(50,0.,200.);
 const Binning kTestBary = Binning:: Simple(25,0.,100.);

Spectrum sSliceCRTDistances("SliceCRTDistances",kTest,loader,kSliceCRTReconstuctorTree,kNoSpillCut); //The spectrum runs the code developed in the SCTripleMatch.h and fills the value of each entry with the value/values of interest  
Spectrum sSliceCRTBarycenters("SliceCRTBarycenters",kTestBary,loader,kSliceCRTReconstuctorBary,kNoSpillCut);  //Spectrum without cuts on the barycenters
Spectrum sSlice2D("Barycenters_CRTDISt",loader,kTestBary,kSliceCRTReconstuctorBary,kTest,kSliceCRTReconstuctorTree,kNoSpillCut); //2D spectrum to see how barycenters and distance to the CRT relate





 
 loader.Go(); //The loader process all
 double factor = 1e20; //Each spectrum comes with a given weight given by CAFANA
 
 TH1* hSliceCRTDistances = sSliceCRTDistances.ToTH1(factor); //The cafana histogram is weighted to the factor 
 //TCanvas *C1 = new TCanvas("C1", "C1",0,65,1200,900);
 //hSliceCRTDistances->Draw();

  
 TH1* hSliceCRTBarycenters = sSliceCRTBarycenters.ToTH1(factor); //The cafana histogram is weighted to the factor 
 //TCanvas *C2 = new TCanvas("C2", "C2",0,65,1200,900);
 //hSliceCRTBarycenters->Draw();
 

 TH2* hSlice2D = sSlice2D.ToTH2(factor);
 //TCanvas *C3= new TCanvas("C3","C3",0,65,1200,900);
 AllBarycentersCRT->Add(hSlice2D);
 //hSlice2D->Draw("colz");
  

 delete hSlice2D;
 std::cout<<"Here Goes file number "<<filenum<<" :D"<<std::endl;
}
TCanvas *C4 = new TCanvas("C4", "C4",0,65,1200,900);
AllBarycentersCRT->Draw("colz"); //here the sum of the spectrums for all the loops gets saved

/*
Here we fill the root tree with the variables of interest defined on the CAF-file 
*/
int evt,npfp,cryostat,Plane,mult,CRTRegion,runnumber=0;
float crtdist,crtdistpandora,flashtime,crtpmtdiff,PCAx,PCAz,PANDORAx,PANDORAz,barycenterdist,DELTA_X_PCA,DELTA_Z_PCA,DELTA_X_PANDORA,DELTA_Z_PANDORA,DIR_X_PCA,DIR_Y_PCA,DIR_Z_PCA,DIR_X_PANDORA,DIR_Y_PANDORA,DIR_Z_PANDORA,trklen,CRTx,CRTz,CRTtime,PCA_MEAN_X,PCA_MEAN_Y,PCA_MEAN_Z,PANDORA_MEAN_X,PANDORA_MEAN_Y,PANDORA_MEAN_Z=0;
double a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36=0;
FILE *fpm = fopen("TripleMatchCoordinates.txt", "r");


TFile *outfile=new TFile("SCTripleMatchCoords.root","RECREATE");
TTree *T = new TTree("T","TripleMatchtree");
T->Branch("evt",&evt,"evt/I");
T->Branch("npfp",&npfp,"npfp/I");
T->Branch("mult",&mult,"mult/I");
T->Branch("trklen",&trklen,"trklen/F");
T->Branch("crtdist",&crtdist,"crtdist/F"); //Check
T->Branch("crtdistpandora",&crtdistpandora,"crtdistpandora/F"); //Check
T->Branch("flashtime",&flashtime,"flashtime/F"); //Check
T->Branch("crtpmtdiff",&crtpmtdiff,"crtpmtdiff/F");
T->Branch("PCAx",&PCAx,"PCAx/F"); //check
T->Branch("PCAz",&PCAz,"PCAz/F"); //check
T->Branch("PANDORAx",&PANDORAx,"PANDORAx/F"); //check
T->Branch("PANDORAz",&PANDORAz,"PANDORAz/F"); //check
T->Branch("Plane",&Plane,"Plane/I"); //check
T->Branch("cryostat",&cryostat,"Cryostat/I"); //check
T->Branch("barycenterdist",&barycenterdist,"barycenterdist/F"); //check
T->Branch("DELTA_X_PCA",&DELTA_X_PCA,"DELTA_X_PCA/F"); //check
T->Branch("DELTA_Z_PCA",&DELTA_Z_PCA,"DELTA_Z_PCA/F"); //check
T->Branch("DELTA_X_PANDORA",&DELTA_X_PANDORA,"DELTA_X_PANDORA/F"); //check
T->Branch("DELTA_Z_PANDORA",&DELTA_Z_PANDORA,"DELTA_Z_PANDORA/F"); //check
T->Branch("DIR_X_PCA",&DIR_X_PCA,"DIR_X_PCA/F"); //check
T->Branch("DIR_Y_PCA",&DIR_Y_PCA,"DIR_Y_PCA/F"); //check
T->Branch("DIR_Z_PCA",&DIR_Z_PCA,"DIR_Z_PCA/F"); //check
T->Branch("DIR_X_PANDORA",&DIR_X_PANDORA,"DIR_X_PANDORA/F"); //check
T->Branch("DIR_Y_PANDORA",&DIR_Y_PANDORA,"DIR_Y_PANDORA/F"); //check
T->Branch("DIR_Z_PANDORA",&DIR_Z_PANDORA,"DIR_Z_PANDORA/F"); //check
T->Branch("CRTx",&CRTx,"CRTx/F"); //check
T->Branch("CRTz",&CRTz,"CRTz/F"); //check
T->Branch("CRTtime",&CRTtime,"CRTtime/F"); //check
T->Branch("CRTRegion",&CRTRegion,"CRTRegion/I"); //check
T->Branch("runnumber",&runnumber,"runnumber/I"); //check
T->Branch("PCA_MEAN_X",&PCA_MEAN_X,"PCA_MEAN_X/F"); //check
T->Branch("PCA_MEAN_Y",&PCA_MEAN_Y,"PCA_MEAN_Y/F"); //check
T->Branch("PCA_MEAN_Z",&PCA_MEAN_Z,"PCA_MEAN_Z/F"); //check
T->Branch("PANDORA_MEAN_X",&PANDORA_MEAN_X,"PANDORA_MEAN_X/F"); //check
T->Branch("PANDORA_MEAN_Y",&PANDORA_MEAN_Y,"PANDORA_MEAN_Y/F"); //check
T->Branch("PANDORA_MEAN_Z",&PANDORA_MEAN_Z,"PANDORA_MEAN_Z/F"); //check

//Adding PCAy,Pandoray,CRT_time,CRT region



while (fscanf(fpm, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&a1,&a2,&a3,&a4,&a5,&a6,&a7,&a8,&a9,&a10,&a11,&a12,&a13,&a14,&a15,&a16,&a17,&a18,&a19,&a20,&a21,&a22,&a23,&a24,&a25,&a26,&a27,&a28,&a29,&a30,&a31,&a32,&a33,&a34,&a35,&a36)!=EOF)
{
 evt= a1;
 npfp= a2;
 flashtime=a3;
 crtpmtdiff= a4;
 cryostat= a5;
 mult = a6;
 Plane= a7;
 trklen= a8;
 barycenterdist= a9;
 crtdist= a10;
 crtdistpandora =a11;
 DELTA_X_PCA= a12;
 PCAx=a13;
 DELTA_X_PANDORA= a14;
 PANDORAx=a15;
 DELTA_Z_PCA=a16;
 PCAz=a17;
 DELTA_Z_PANDORA=a18;
 PANDORAz=a19;
 DIR_X_PCA=a20;
 DIR_Y_PCA=a21;
 DIR_Z_PCA=a22;
 DIR_X_PANDORA=a23;
 DIR_Y_PANDORA=a24;
 DIR_Z_PANDORA=a25;
 CRTx=a26;
 CRTz=a27;
 CRTtime=a28;
 CRTRegion= a29;
 PCA_MEAN_X = a30;
 PCA_MEAN_Y = a31;
 PCA_MEAN_Z = a32;
 PANDORA_MEAN_X = a33;
 PANDORA_MEAN_Y = a34;
 PANDORA_MEAN_Z = a35;
 runnumber= a36;


T->Fill();
}
    
    
 T->Write();
AllBarycentersCRT->Write("Barycenters-Distances"); //And save the plot of the coincidences
   




 outfile->Close();


}
