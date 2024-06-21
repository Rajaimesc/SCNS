


#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/EnsembleRatio.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
//#include "sbnana/CAFAna/Cuts/TruthCuts.h"
//#include "LocalCuts.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Analysis/ExpInfo.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h" //after v09_44 release
#include "sbnana/CAFAna/Systs/NuMIFluxSysts.h"
#include "sbnana/SBNAna/Cuts/NumuCuts.h"
//#include "FlashesHelper.h"

#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/SBNAna/Vars/Binnings.h"
#include "sbnana/SBNAna/Vars/NueVars.h"
#include "sbnana/SBNAna/Vars/NumuVars.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"
#include "sbnana/SBNAna/Cuts/TruthCuts.h"
#include "sbnana/SBNAna/Cuts/NumuCuts.h"
#include "TVector3.h"
#include <algorithm>

#include <fstream>
#include <iostream>
#include <vector> 
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "stdio.h"
#include "TProfile.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
using namespace ana;

/*We want to save the output of our variables in a text file since CAFana works iterating over each event and saves only the variables we want to fill
  into an histogram called spectrum
 */ 
ofstream MyFile("TripleMatchCoordinates.txt");
ofstream MyFile2("PROJECTEDVALUESPCA.txt");

/////////////////////////////////////////////////////////////////////////////////////////////////////////////// CONSTANTS //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//ICARUS is designed as a Liquid Argon TPC from which the ionized particles reach the cable that collect the charges inside the volume. 
//This proccess of collection and reconstruction is set up to a certain constants as drift velocity and the tics to clock conversion of the electronics
double const vdrift=0.157; // The drift velocity in the projection chamber is set to 0.157 cm/us Due to the electric field in it // 0.173 with clock 850 tics in MC
double const tics=0.4; // each clock tick in the electronics corresponds to 0.4 us(400 ns)

//ICARUS is divided into two cryostats(East and West) that have coordinates -358.63 to -61.7 cm in X for the first one and 61.7 to 358.73 cm for the second
//This values are kept as constants to ensure that the reconstruction done by the triple Match is still in the fiducial volume.

double const minLimitW=61.7; //cm Anode West-East position 
double const maxLimitW=358.73; //cm Anode West-West position
double const minLimitE=-61.7; //cm Anode East-West position
double const maxLimitE=-358.73; //cm Anode East-East position
double const cathW=210; // cm Cathode W position. Each cryostate has a cathode for which the the electric field is generated located at the middle of the x coordinate
double const cathE=-210; // cm Cathode E position. Each cryostate has a cathode for which the the electric field is generated located at the middle of the x coordinate
double const exc=2; // cm max displacement out of the cathode

bool isInFV (double x, double y, double z) //Fiducial volume of the detector
{
  if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false; //A simple boolean designed to check wheter a coordinate is within the fiducial volume of icarus wihthin an error.

  return (( ( x < -61.94 - 25 && x > -358.49 + 25 ) ||
			( x >  61.94 + 25 && x <  358.49 - 25 )) &&
		  ( ( y > -181.86 + 25 && y < 134.96 - 25 ) &&
		  ( z > -894.95 + 30 && z < 894.95 - 50 ) ));
}


//SpillCutsare cuts put into the cut for which the spectrum iterator skips if the condition is nor fullfiled (same as the ones in root histograms).

//This cut is used for testing and reduces the iteration to one just even in the CAF-tree.
const SpillCut kJustOneEventCut([](const caf::SRSpillProxy* sr){
  int event= sr->hdr.evt;
  if(event == 30||event == 8 ||event == 21 || event == 15 ||event == 17 || event == 3 ) return true; //18908 first and original  24552,24924,25608,23964 MC= 30,8,21,15,32,17,3
  return false;
});



////////////////////////////////////////////////////////////////////////////////////////////////////////// LETS PLAT WITH STRUCTS/CLASSES///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
Since the process of the TripleMatch implies looking  the information that each subsystem takes separately and work with them together a way to manipulate this info
is to store them in given structures for the variables of interest given each detector/object.
*/

//CrossPoint storages the projections of the selected and drifted track onto the Cosmic Ray Tagger Plane
struct CrossPoint {
    double X;
    double Y;
    double Z;
};


//DriftedTrack is a collection of the spatial coordinates  of the track after it is drifted on the X coordinates along its initial and end coordinates.
// Outbound defines the number of points that fall outside the fiducial volume.
struct DriftedTrack {
    std::vector<double> spx; // Drifted Track Hit Points X coordinate
    std::vector<double> spy; // Drifted Track Hit Points Y coordinate
    std::vector<double> spz; // Drifted Track Hit Points Z coordinate
    int outbound; // Number of hit points out of the logical volume of the TPC
    double startx; // Drifted Track StartX coordinate
    double endx; // Drifted Track EndX coordinate
};

//Direction saves the Director cosines and mean values for which the line that defines the track are constructed from fitting the points.
struct Direction {
    double dirx; // Direction of the Track: X 
    double diry; // Direction of the Track: Y 
    double dirz; // Direction of the Track: Z 
    double meanx; // Mean Point of the Track: X 
    double meany; // Mean Point of the Track: Y 
    double meanz; // Mean Point of the Track: Z 
};

// CRT hit saves the coordinates that we need to take from the the Cosmic ray tagger system as: the time of the signal, position, the region of the CRT
// Also the CRT system and the photoelectron quantity that triggers the signal.
struct CRThit {
	double t1; // CRT Hit time
	double x; // CRT Hit X Coordinate
	double y; // CRT Hit Y Coordinate
	double z; // CRT Hit Z Coordinate
	int region; // CRT Hit region
	int sys; // CRT Hit subsystem (0 is Top, 1 is Side)
  double pe; // CRT Hit amplitude [pe]
};

/* 
   The optical flash structure saves the variables  from the optical system that serves to perform the TripleMatch: the type of the flash, the time
   the total apmplitude, the cryo in which the flash is registered, the flash barycenters, and wheter the flash was registered in beam time or in Collection time(Beam + an enclosing interval)
   The flash outputs are already connected to a CRT-hit from the DAQ since flash and CRT-times can be assigned to the same particle it they occur within +/- 100 ns
   All of the variables are storaged into vectors as the loop that CAFAna does is over the TPC objects referred as slices. The optical flashes and CRT-hits live in a different structre of the CAF file
   for this reason, first the CAF loops over the Opflash-CRT and storages those values in vectors for the iterations in the TPC.
*/ 
struct OpFlash {
    std::vector<int> flash_types; // Classification of the Flashes according to CRT-Match: 0 is no CRT-Match, 1 is Top-CRT signal before OpFlash, 2 is side-CRT signal before OpFlash, 3 is matched with a top CRT first, then a flash then a side CRT, 4 is matched with a top CRT after an optical flash, 5 is matched with a side CRT after an optical flash, 6 is matched with multiple top CRT-hits before the flash, 7 is matched with multiple top CRT-hits and 1 side CRT-hits before the flash, 9 is all other cases.

    std::vector<double> flash_times; // Flash Times measured by the photomultiplier systems.
    std::vector<double> pes; // Photo-electron amplitudes measured by the PMT systems
    std::vector<int> cryos; //  0 East, 1 West
    std::vector<double> Barycenters_X; // Flash barycenters X coordinate
    std::vector<double> Barycenters_Y; // Flash barycenters Y coordinate
    std::vector<double> Barycenters_Z; // Flash barycenters Z coordinate
    vector<bool> inGate_times; // If flash times in gate window
    vector<bool> inBeam_times; // If Flash time in beam spill window
    std::vector<int> counts_perFlash; // number of tracks for which this flash was the best match-candidate
    
    std::vector<std::vector<std::vector<double>>> CRTPos; // Eventually, position of the CRT Hits associated with the Flash-Track match //There can be more than one CRT-Candidate per flash
    std::vector<std::vector<double>> CRTTime; // Eventually, time of the CRT Hits associated with the Flash-Track match.
    //////////////////////////////////
    
    //The CRT variables matched to the flash.
    std::vector<std::vector<CRThit>> CRT; // Eventually, vector of CRT struct matched with this Flash
    std::vector<int> flash_id; // Id of the Optical Flashes for matching with the OpFlash_rawMap
};



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////// HERE WE PUT THE USEFUL FUNCTIONS FOR THE ANALYSIS ///////////////////////////////////////////////////////////////////////////////////////////



int DriftDirection(int TPC) //Dirft direction gives the direction of the drift  for each TPC inside each cryostat, regions 0/1 are left of the cathode (so they drift in x negative) and regions 3/4 are right of the cathode (so they drift in x  positive)
{
  int drift=0;
  if(TPC==0 || TPC== 1) // Drift direction for this TPCs are left of the cathode (to X negative)
  {
    drift= -1;

  }
  else if(TPC==2 || TPC==3) // Drift direction for this TPCs are left of the cathode (to X positive)
  {
    drift=1;

  }
  return drift;

}
float delta_bar_expected_from_mc(float x, int cryo){ // The distribution of electric field in the cryostats is not uniform so to correctly assign the weights at the end points of the cryostats the following empirical formula was developed following previous studies
    float value=0;
    if(cryo==0)
        {
        if(x< -500)
        {
        value=    + 9371.09 +  74.2292*x+ 0.233642*x*x+ 0.000365071*x*x*x+ 2.83091e-07*x*x*x*x+ 8.72829e-11*x*x*x*x*x;
        }
    else if(x> 500)
        {
        value= -1047.96 +12.1278*x - 0.0496555*x*x +9.45929e-05*x*x*x  -8.59464e-08*x*x*x*x+3.03292e-11*x*x*x*x*x;
        }
        }
    if(cryo==1)
        {
        if(x< -500)
        {
        value=    -7313.39 -54.5114*x-0.160099*x*x-0.00023108*x*x*x-1.63539e-07*x*x*x*x-4.51835e-11*x*x*x*x*x;
        }
    else if(x> 500)
        {
        value=5983.29 -41.805*x +0.113756*x*x -0.000149932*x*x*x+9.48963e-08*x*x*x*x-2.26312e-11*x*x*x*x*x;
        }
        }
        return value;
}


/* The Barycenters are collected by the 3 GetopFlashes_center functions
 */
std::vector <double> GetopFlashes_center(const caf::SRSpillProxy* sr) //CAFana loops over all the standard record (event by event)
{ std::vector <double> barycenters_z;
  double center_z=0;
  std::cout<<"Number of flashes in this event"<<" "<<sr->nopflashes<<std::endl;
  for(auto const& opflash : sr->opflashes){ 
   center_z=opflash.center.z; //Each flash of the event gets saved in a vector for which it will be associated with the TPC-CRT.
   barycenters_z.push_back(center_z);
  
  }
  return barycenters_z;

}
std::vector <double> GetopFlashes_center_x(const caf::SRSpillProxy* sr) //CAFana loops over all the standard record (event by event)
{ std::vector <double> barycenters_x;
  double center_x=0;
  std::cout<<"Number of flashes in this event"<<" "<<sr->nopflashes<<std::endl;
  for(auto const& opflash : sr->opflashes){
   center_x=opflash.center.x; //Each flash of the event gets saved in a vector for which it will be associated with the TPC-CRT.
   barycenters_x.push_back(center_x);
  
  }
  return barycenters_x;


}
std::vector <double> GetopFlashes_times(const caf::SRSpillProxy* sr) //Same as with the barycenters the times are saved  after scanning over the pmt signals of the event
{ std::vector <double> times;
  double time=0;
  for(auto const& opflash : sr->opflashes){
   time=opflash.time;
  times.push_back(time);
  
  }
  return times;

}
std::vector <double> GetopFlashes_cryo(const caf::SRSpillProxy* sr) //Same as with the barycenters the cryostat information is saved  after scanning over the pmt signals of the event
{ std::vector <double> times;
 std::vector <double> cryos;
  double cryo=0;
  for(auto const& opflash : sr->opflashes){
   cryo=opflash.cryo;
 cryos.push_back(cryo);
  
  }
  return cryos;

}

std::vector<double> ChargeBarycenter(const caf::Proxy<caf::SRPFP>& ipfp)  //The Barycenter of charge is computed at code level using the coordinates in the TPC wires and the charge amplitude deposited in it.
{ 
  double average_x_charge=0; 
  double average_y_charge=0;
  double average_z_charge=0;
  double  x_charge=0;
  double  y_charge=0;
  double  z_charge=0;
  double total_charge=0;
  std::vector<double> vector_active;

  
    for(std::size_t ihit(0); ihit < ipfp.trk.calo[2].points.size(); ++ihit){ //first we lopp  over all the hits of the track in the plane of collection. ICARUS has 3 planes:induction1(0),induction2(1) and collection(2). The collection plane is used for reconstruction since is the less noisy of the 3  
    double x = ipfp.trk.calo[2].points[ihit].x; //each of the coordinates is collected from the hit sored in the points vector
		double y = ipfp.trk.calo[2].points[ihit].y; //each of the coordinates is collected from the hit sored in the points vector
		double z = ipfp.trk.calo[2].points[ihit].z; //each of the coordinates is collected from the hit sored in the points vector
		double integral_charge= ipfp.trk.calo[2].points[ihit].integral; //The integral of the carge for the hit is used as the weight for the barycenter
      
        int size=  ipfp.trk.calo[2].points.size(); // the total number of hits that form the track

		if(std::isnan(z)) continue; //Sometimes the hits are saved as nan as they are badly reconstructed, we avoid those ones
		if(z> -1000 && z <1000 && isInFV(x,y,z)){ // we check that all of the hits are within the fiducial volume.
		 x_charge+= x*integral_charge; //calculate the value of the charge in the x plane
		 y_charge+= y*integral_charge;  //calculate the value of the charge in the y plane
		 z_charge+= z*integral_charge;  //calculate the value of the charge in the z plane
		total_charge+= integral_charge;  //The total charge is the sum of the integral charge on each hit
    }  


    }
    average_x_charge= x_charge/total_charge; // The average charge is the just calculated by the charge in x over the total charge
    average_y_charge= y_charge/total_charge; // The average charge is the just calculated by the charge in y over the total charge
    average_z_charge = z_charge/total_charge; // The average charge is the just calculated by the charge in z over the total charge
    vector_active= {average_x_charge,average_y_charge,average_z_charge};
    return vector_active;


 
 } 

Direction PCAfit (const caf::Proxy<caf::SRPFP>& ipfp, float time,int cryo){ // The PCA fit is the method of track reconstruction that is used along the one done by Pandora to perform the TripleMatch and further projection on the CRT-Plane
  //The PCA fit gets the information for a track, stored in the pfp object of CAFana and from whicht the TPC signals are taken, the flash time from which the track was drifted and the cryostat in which the track was reconstructed 
   std::vector <double> x; //Position vector in x
   std::vector <double> y; //Position vector in y
   std::vector <double> z; //Position vector in z
   int size=0; //The size of the vector depends on the number of hits that the track has in the collection plane
   int Nhits=0; //Number of hits that form the track
   double t0 = ipfp.t0; //time saved for the TPC track, for non cathode crossers that time is Nan
   
  if(ipfp.trk.calo[2].nhit != -999) // Pandora doesnt save the hit information of all the tracks since it rejects most of the cosmic rays, for this reason we make sure the information exists.
   {
    Nhits= ipfp.trk.calo[2].nhit; //Nhits is a variable saved in the DAQ and reconstruction of the TPC
   }
   else 
   {
    Nhits= 0; // This conditional is a sort of "sainitizer" to avoid segmentation faults due to the fact that sometimes the pandora reconstructor sets the hit counter at -99 for empty tracks.
              //This way the sizes are all set to number >=0 
   }
  

    for(std::size_t ihit(0); ihit <ipfp.trk.calo[2].points.size() ; ++ihit){ //Loop over all the hits in the collection plane that form the track in the Tpc
    // Sub x, sub y and sub z are auxiliar variables from which the information of each hit the collection plane is stored
    double suby=0;
    double subz=0;
    double recoX=0; //Reco X stores the information of the drifted X coordinate after being moved to the flash time correspondent position
    double hit_time = ipfp.trk.calo[2].points[ihit].t; // Stored time of the hit in the electronic readout system
    int tpc= ipfp.trk.calo[2].points[ihit].tpc; // The tpc plane in which the track is reconstructed, information taken by PANDORA.
    int drift= DriftDirection(tpc); //Knowing the tpc plane the direction of the drift for this track is computed.
    recoX = ipfp.trk.calo[2].points[ihit].x; //Cathode crosser tracks are correctly reconstructed in the detector, so they are not drifted.
    suby = ipfp.trk.calo[2].points[ihit].y; //The drift occurs only in x, so for this reason y and z remain the same
    subz=  ipfp.trk.calo[2].points[ihit].z; //The drift occurs only in x, so for this reason y and z remain the same
    x.push_back(recoX); //once the drift happens all of the coordinates of the Track are stored in a vector for x, y and z.
    y.push_back(suby);  //once the drift happens all of the coordinates of the Track are stored in a vector for x, y and z.
    z.push_back(subz);  //once the drift happens all of the coordinates of the Track are stored in a vector for x, y and z.

  }
  
   size=x.size(); //Since every hit has coordinates x, y and z the size of the vectors is the same
    double xavg=0, yavg=0, zavg=0; // For the next step the PCA needs to get the averages of the positions in each vector
    for(int k=0; k<size; k++) {
        xavg+=x[k]; yavg+=y[k]; zavg+=z[k]; //loop over the values of the vectors to get the total value
    }
    
    xavg=xavg/size; //average value as the total/size
    yavg=yavg/size; //average value as the total/size
    zavg=zavg/size; //average value as the total/size
    float x2=0, y2=0, z2=0, xiy=0, xiz=0, yiz=0; //Is necessary also to determine the dispersion of the points in each vector, aka how far they are from the mean
    for(int k=0; k<size; k++){
        float xadj = x[k] - xavg; //for each value in the vector the distance from the mean is calculated
        float yadj = y[k] - yavg; //for each value in the vector the distance from the mean is calculated
        float zadj = z[k] - zavg; //for each value in the vector the distance from the mean is calculated
        x2+=xadj*xadj; // The variances of the variables are then computed for the PCA fit
        y2+=yadj*yadj; // The variances of the variables are then computed for the PCA fit
        z2+=zadj*zadj; // The variances of the variables are then computed for the PCA fit
        xiy+=xadj*yadj; // The covariances of the variables are also needed for the PCA fit
        xiz+=xadj*zadj; // The covariances of the variables are also needed for the PCA fit
        yiz+= yadj*zadj; // The covariances of the variables are also needed for the PCA fit
    }//end loop calculating covariance matrix elements
    x2=x2/size; //The elements of the covariance matrix are also weighted by the mean
    y2=y2/size; //The elements of the covariance matrix are also weighted by the mean
    z2=z2/size; //The elements of the covariance matrix are also weighted by the mean
    xiy=xiy/size; //The elements of the covariance matrix are also weighted by the mean
    xiz=xiz/size; //The elements of the covariance matrix are also weighted by the mean
    yiz=yiz/size; //The elements of the covariance matrix are also weighted by the mean

    //Now All the variance and covariance elements are placed into the covariance matrix,for which the PCA analysis is performed
    TMatrixD covmat(3,3); //covariance matrix is computed place by place
    covmat[0][0] = x2; covmat[1][1]=y2; covmat[2][2]=z2;
    covmat[0][1] = xiy; covmat[1][0] = xiy;
    covmat[0][2] = xiz; covmat[2][0] = xiz;
    covmat[2][1] = yiz; covmat[1][2] = yiz;
    TMatrixDEigen thiscov(covmat); // The Matrix goes through eigen to get EigeValues and EigenVectors. 
    TMatrixD this_eval = thiscov.GetEigenValues(); //eigenValues are taken from the matrix. The principal components of the Track are obtained by means of reducing the correlations between the variables and therefore getting the most signal of the components
    TMatrixD this_evec = thiscov.GetEigenVectors(); //eigenVectors correspond to the principal directions that refer to the lineal combination of the coordinates that defines the principal components.
    if(this_eval.GetNrows()!=3 || this_eval.GetNcols()!=3 || this_evec.GetNrows()!=3 || this_evec.GetNcols()!=3){ 
        std::cout << "evals/evects wrong size! continuing....\n"; // A Security check to see if both matrices are 3X3 as expected
    }

    int max_eval = -999999; int maxevalpos = -1;
    for(int k = 0; k < 3; k++){
        if(this_eval[k][k]>max_eval){
            max_eval = this_eval[k][k];
            maxevalpos = k;  //As the EigenValues refer to the magnitude of the variances in the new principal components, the one with the larger variance is chosen as a larger variance means a bigger difference between the coordinates
        }
    }//end loop looking for best eval
    Direction ThisDirection= {this_evec[0][maxevalpos] , this_evec[1][maxevalpos], this_evec[2][maxevalpos], xavg, yavg, zavg}; // The Directions of the reconstructed  tracks are taken as the eigenvecotrs that correspond to the linear combination of bigger variance. This directions define the director cosines of the track, which are saved along the mean points in x,y and z.
    return ThisDirection; //The direction of the track gets saved in a struct called direction for further analysis and projection
}


Direction Pandorafit(const caf::Proxy<caf::SRPFP>& ipfp,float time,int cryo) //For comparison of the performance of the PCA reconstruction with respect of the Pandora reconstruction ones the director cosines that pandora generate are also taken from the event data
{
  double suby=0; //Same as with the PCA fit, we take the coordinates and director cosines for each point of the track
  double subz=0;
  double dirx=0;
  double diry=0;
  double dirz=0;
  std::vector <double> x; //This vector save the value of the coordinates on each point of the track hit
  std::vector <double> y;
  std::vector <double> z;
   double recoX=0;
   double t0 = ipfp.t0;
   
  
  for(std::size_t ihit(0); ihit < ipfp.trk.calo[2].points.size(); ++ihit){
   int tpc= ipfp.trk.calo[2].points[ihit].tpc; 
   double hit_time = ipfp.trk.calo[2].points[ihit].t;
   int drift= DriftDirection(tpc);
   if(isnan(t0)==false)
   {
    recoX = ipfp.trk.calo[2].points[ihit].x; //Cathode crossers have tehir time reconstructed correctly so they dont need to be drifted.
   }
   if(isnan(t0)== true)
   {
    recoX =(hit_time- 846 - (time/tics))*vdrift; // Non cathode crossers are reconstructed at the trigger time so they are shifted to the flash time of the hit
   }
  

    suby = ipfp.trk.calo[2].points[ihit].y;
    subz=  ipfp.trk.calo[2].points[ihit].z;
    x.push_back(recoX);
    y.push_back(suby);
    z.push_back(subz);
    

  }
  int size=x.size();
  double xavg=0, yavg=0, zavg=0;
  for(int k=0; k<size; k++) 
  {
    xavg+=x[k]; yavg+=y[k]; zavg+=z[k];
  }
  xavg=xavg/size;
  yavg=yavg/size;
  zavg=zavg/size;
/* here is a problem with the reconstruction by hand of:
 -Barycenters
 -Mean coordinates
 -DirCos(using the PCA).
  There should be then a mistake in the way the coordinates in the calohits are selected for then manipulated in both tha matrices and the averages, because it only works on the first iteration. The first version of the TripleMatch works then with  Pandora Values and trk.start and trk.end. 
  Only the first iteration works and gives a very similar reconstruction (+/-5 cm per each track)  AvgX, AvgY and AvgZ are then defined under the pandora given variables mentioned earlier.
  
*/




 double AvgX=  ipfp.trk.start.x; // The director cosines that are saved for the Pandora reconstruction at the start and end points. I take the start ones since they dont experience multiple scattering as much as the end ones
 double AvgY=  ipfp.trk.start.y;
 double AvgZ=  ipfp.trk.start.z;

  dirx= ipfp.trk.dir.x; //Director cosines of the starting points in pandora
  diry= ipfp.trk.dir.y;
  dirz= ipfp.trk.dir.z;




 Direction PandoraDir= {dirx,diry,dirz,AvgX,AvgY,AvgZ};  //Saved in the struct that is used for doing the projection after
  



 return PandoraDir;
}

DriftedTrack Drift_tracks_outbound (const caf::Proxy<caf::SRPFP>& ipfp, double time,int cryo) { //The first constraint that is made over the drifted tracks for the triplematch is that the drifted track should still be contained within the fiducial volume. Drift_tracks outbund checks that feature
   
   /*As the structure of CAFana is an iterative one, in which each functions iterates over a particle reconstruction, instead of saving a vector for each track on the position
    Is more efficient to just iterate over the positions and director cosines for each event on each function and then check each variable of interest
   */
   int outBound=0; //Outbound tells us how many points are outside of the fiducial volume of the detector
   double cath=0; //Coordinates of the cathode
   double maxLimit=0; //Maxlimit of the cryostat, changes if 0 or 1
   double minLimit=0;//Minlimit of the cryostat, changes if 0 or 1
   double startx= ipfp.trk.start.x; //Start point of the track 
   double endx= ipfp.trk.end.x; //End point of the track 
   double t0= ipfp.t0;
   std::vector<double> recX; //Vector that stores the drifted X coordinates 
   std::vector<double> recY; //Stores the Y coordinates
   std::vector<double> recZ; //Stores the Z coordinates
   
   
  /*
   To check wether the drifted track is on the fiduacial veolume or not we define the limits on each cryostat
  */
   if(cryo==0) //the east cryo
    { 
      cath= cathE; //The value of the cathode is set to -360 cm
      maxLimit= maxLimitE;//Max coordinate in x is -58.
      minLimit= minLimitE; //Minimum coordinate in x is -210 cm

    }
   else if(cryo==1) //the west cryo
   { 
    cath= cathW; //The value of the cathode is set to 210 cm
    maxLimit=maxLimitW; //Max coordinate in x is 360
    minLimit=minLimitW;  //Min coordinate in x is 58

   }
   
    
  for(std::size_t ihit(0); ihit <ipfp.trk.calo[2].points.size(); ++ihit)
   {
    double realY= ipfp.trk.calo[2].points[ihit].y; //Each of the points reconstructed in the TPC is collected for each vector
    double realZ= ipfp.trk.calo[2].points[ihit].z;
    double x= ipfp.trk.calo[2].points[ihit].x;
    double hit_time= ipfp.trk.calo[2].points[ihit].t; //the time of the hit given by the photomultiplier system too
    int tpc= ipfp.trk.calo[2].points[ihit].tpc; // the tpc type is takes (EE(0),EW(1),WE(2),WW(3))
    double recoX=0;
    double planeX=0;
    double realX=0;
    if(std::isnan(x)) continue;

    int drift= DriftDirection(tpc); //Drif direction determines wheter the track drifts to east or west part of each cryostat
    //will add the driftdirection here.
    if(isnan(t0)==false)
    {
      realX= x; //Cathode crossers have a good reconstruction of x and therefore are not needed to correct
      continue;
    }
    if(isnan(t0==true))
    {
     recoX=(hit_time- 846 - (time/tics))*vdrift; //The reconstruction takes the time trigger time(for non cathode crossers) and correct the drifted hit by means of the flash time 
    }
    
    if(drift==-1 && cryo==1) // In this case we  are on the  west cyostat, so x>0 and the point is drifting towards the east wall
    {
      planeX=minLimit; 
      realX= planeX + recoX; 

    }
    if(drift==1 && cryo==1) // In this case we  are on the  west cyostat, so x>0 and the point is drifting towards the west wall
    {
      planeX=maxLimit;
      realX=planeX - recoX;

    }
    else if(drift==-1 && cryo==0) // In this case we  are on the  east cyostat, so x<0 and the point is drifting towards the east wall
    {
      planeX= maxLimit; 
      realX=planeX+recoX;


    }
    else if(drift== 1 && cryo==0) // In this case we are on the  east cyostat, so x<0 and the point is drifting towards the east wall
    {
      planeX=minLimit; 
      realX=planeX-recoX;

     }
     // To check if the track is still fully contained is necessary to check the coordinates of each cryostat and wheter every point is contained or not
    if(cryo==1)  
    {
      if(tpc==0 || tpc==1) 
      {
        if(realX>(cath + exc) || realX<(minLimit-exc)) { //the excess quantity refers to an uncertainty of 2cms on the reconstruction . For the first 2 tpcs if the point is below the east limit and the cathode is outbound  
          outBound++;
        }

      }
      else if(tpc==2 || tpc==3) //For the last 2 tpcs if the point is outside the west limit or behind the cathode is outbound 
      {
        if(realX < (cath- exc) || realX>( maxLimit + exc))
        {
          outBound++;

        }

      }

    }
    else if(cryo== 0)
    {
      if(tpc==0 || tpc==1) // For the first 2 tpcs  
      {
        if(realX>(cath + exc) || realX<(maxLimit-exc)) //Every coordinate reconstructed from the cathode to the higher limit is outbound for the first two TPCs
        {
          outBound++;

        }

      }
      if(tpc==2 || tpc==3) //For the last 2 TPCs
      {
        if(realX < (cath- exc) || realX>( minLimit + exc)) //Every coordinate reconstructed from the lower limit to the cathode  is outbound for the first two TPCs
        {
          outBound++;
        }


      }

    }
    recX.push_back(realX); //The values in X are stored afterd drifted
    recY.push_back(realY); //Values in Y and Z remain the same
    recZ.push_back(realZ); //Values in Y and Z remain the same
    if(x == ipfp.trk.start.x && realY==ipfp.trk.start.y && realZ==ipfp.trk.start.z ) //For the further projection on the TPC plane the start and end coordinates of the track are also saved
    {
      startx=realX;
    }
    if(x == ipfp.trk.end.x && realY==ipfp.trk.end.y && realZ==ipfp.trk.end.z )
    {
      endx=realX;
    }
   }



   DriftedTrack thisDriftedTrack = {recX,recY,recZ,outBound,startx,endx}; //Finally the drifted track is saved in a struct that keeps information about the coordinates,star,end and outbound


  
    
    return thisDriftedTrack;
}
/*The information of some hits in the CRT detector is already matched with the optical flashes coordinates.
This function does an scan of the CRT hits matched with a flash for a given event and stores the information so it is stored to be handled along the track information
The information is saved in a struct called CRThit from which each of the variables of interest is called and used in relation the reconstructed tracks
*/
std::vector<CRThit> MatchedCRThits(const caf::Proxy<caf::SRCRTPMTMatch>& crtpmt) //Scan through the library of the Caf file that stores the information of the CRT 
{
  /* Thee coordinates of interest fot the projection are the time and position of all the hits on the CRT plane, the region where it happened 
  and the signal that sets the flash
  */
  double t1=0;
  double x=0;
  double y=0;
  double z=0;
  int region=0;
  int sys=0;
  double pe=0;
  std::vector<CRThit> hits_vector;
  pe = crtpmt.flashPE;
  for(std::size_t icrt(0); icrt < crtpmt.matchedCRTHits.size() ; ++icrt) // The hits are already stored by event on a vector, so a loop on each CRT hit is done
  { 
    t1= crtpmt.matchedCRTHits[icrt].time; 
    x = crtpmt.matchedCRTHits[icrt].position.x;
    y = crtpmt.matchedCRTHits[icrt].position.y;
    z = crtpmt.matchedCRTHits[icrt].position.z;
    region = crtpmt.matchedCRTHits[icrt].region;
    sys =  crtpmt.matchedCRTHits[icrt].sys;
    CRThit This_hits = {t1,x,y,z,region,sys,pe}; //All values for each hit are saved in the CRT_hit struct
    hits_vector.push_back(This_hits); //Since sometimes there are multiple CRT hits matched with the flash time the hits are stored in a vector
  }
    
     




return hits_vector;


};

OpFlash CRTPMTValues(const caf::SRSpillProxy* sr) //CRPMTValues stores not only the values of the CRT hits but also of the flash that is associated with the CRT hit like the barycenter position
{ 
  // Variables regarding the flashes intrinsic info
  std::vector<int> flash_types;  // The classification of the flash refers to the signal of the flash received before or after the CRT signal and how many hits per flash. Goes from 0 to 9
  std::vector<double> flash_times;
  std::vector<double> flash_pe;  //Photoelectrons that set the trigger 
  std::vector<int> cryos;  //Cryostate of the flash
  std::vector<double> Bary_x;  //Barycenter of light positions
  std::vector<double> Bary_y; 
  std::vector<double> Bary_z; 
  vector<bool> flashes_inGate;  //If the flashes of the event occut within the Gate of time acquistion (10 microsendonds)
  vector<bool> flashes_inBeam;  //If the flases occur within the beam time (2.2 microseconds)
  std::vector<int> flashes_id; //id of the flash 
  //Variables that are filled after the track and CRT matching
  std::vector<int> counts_perflash;
  std::vector<std::vector<std::vector<double>>> Candidate_CRTPositions; 
  std::vector<std::vector<double>> Candidate_CRTTimes;

 //Variables that regard te Op_Hits of the flash
  std::map<int,std::vector<double>> BaryPlanes1; 
  std::map<int,std::vector<double>> BaryPlanes2; 
  std::vector<double> sumAmpPlanes1; 
  std::vector<double> sumAmpPlanes2; 

  //Variables that contain the CRHits of the matched CRT-PMT  
  std::vector<std::vector<CRThit>> Flashes_CRT; 
  

//Here the things to fill (a lot of initialized maps will be empty the same)
  int type=0; 
  double time=0; 
  double pe=0; 
  int cryo=-1; 
  double X=0; 
  double Y=0; 
  double Z=0; 
  bool inGate=0; 
  bool inBeam=0; 
  int counts=0; 
  std::vector<std::vector<double>> CRTPos; 
  std::vector<double> CRTTime;  
  std::vector<double> BaryPlane1; 
  std::vector<double> BaryPlane2; 
  double sumAmpPlane1=0; 
  double sumAmpPlane2=0; 
  std::vector<CRThit> CRT; 
  int flash_id=0;

 //So we loop over the CRTPMTValues 
 for(auto const& crtpmt : sr->crtpmt_matches){ 
  int counter=0;
  //All the information about the flashes is collected from this iteration
  type = crtpmt.flashClassification;
  inBeam= crtpmt.flashInBeam;
  inGate= crtpmt.flashInGate;
  pe =  crtpmt.flashPE;
  X =   crtpmt.flashPosition.x;
  Y =   crtpmt.flashPosition.y;
  Z =   crtpmt.flashPosition.z;
  time = crtpmt.flashTime_us;
  flash_id= crtpmt.flashID;

  flash_types.push_back(type);
  flash_times.push_back(time);
  flash_pe.push_back(pe);
  flashes_inGate.push_back(inGate);
  flashes_inBeam.push_back(inBeam);
  Bary_x.push_back(X);
  Bary_y.push_back(Y);
  Bary_z.push_back(Z);
  if(X<0)
  {
    cryo=0;
    cryos.push_back(cryo);

  }
  if(X>1)
  {
    cryo=1;
    cryos.push_back(cryo);

  }

  CRT = MatchedCRThits(crtpmt); //The CRT  hits associated to the given flash are collected by calling the MatchedCRThits function defined before
  Flashes_CRT.push_back(CRT); //Each CRT hit is stored in a vector for each flash
  counter++;
}
  
  OpFlash Matches= {flash_types,flash_times,flash_pe,cryos,Bary_x,Bary_y,Bary_z,flashes_inGate,flashes_inBeam,counts_perflash,Candidate_CRTPositions,Candidate_CRTTimes,Flashes_CRT,flashes_id}; //All the information from the light and cosmic systems is then avaliable to be worked at code level.

  return Matches;




}

/* Once a track gets assigned a flash time by means of the barycenter clossenes for the TCP and PMT systems the CRT hits that correspond to the flash
are also transposed to the track, where the projection is then calculated and compared with the signal in the detector.
Determine plane assigns the type of the CRT (top,side or bottom) and the coordinate that remains constant (y for the top, x for the side, z for the bottom)
*/
std::pair<int, double> DeterminePlane(CRThit CRThit){ 
    int CRTPlane; //Plane gets the type of CRT
    double CRTPos; //Position the coordinate of the constant coordinate 

    if(CRThit.region==30) { // 30 defines the top CRT 
        CRTPlane=0;
        CRTPos=CRThit.y;
    }
    /*
    From 31 to 45 defines the side CRT 
    */
    else if(CRThit.region==31 || CRThit.region== 32 || CRThit.region==40 || CRThit.region==41 || CRThit.region==42 || CRThit.region==43 || CRThit.region==44 || CRThit.region==45) {
        CRTPlane=1;
        CRTPos=CRThit.x;
    }
    else {
        CRTPlane=2;
        CRTPos=CRThit.z;
    }
    return std::make_pair(CRTPlane, CRTPos);
}

/*
The cross point struct saves the coodinates projected on the track to the CRT plane matched with the flash
CalculateProjection uses the director cosines and the mean points of the track to elongate the track to the CRT plane 
*/
CrossPoint CalculateProjection(double dirx, double diry, double dirz, double x0, double y0, double z0, double position) 
{
  double Lambda = (position-x0)/dirx; //The lambda defines the step of the track
  double PosAtY = (Lambda*diry+y0); //Projection in y
  double PosAtZ = (Lambda*dirz+z0); //Projection in z

  CrossPoint CrossingPoint = {position, PosAtY, PosAtZ};
  return CrossingPoint;
}


/*
 The crossing point gets determined for the stored track and the resultant crossing point depends of the CRT plane.
*/

CrossPoint DetermineProjection(Direction thisDir, int plane, double position) 
{
 CrossPoint CrossingPoint;
  if(plane==0) // Plane at Y e.g. Top CRT Horizontal Plane
  {
    CrossPoint thisCase=CalculateProjection(thisDir.diry, thisDir.dirx, thisDir.dirz, thisDir.meany, thisDir.meanx, thisDir.meanz, position);
    CrossingPoint={thisCase.Y, thisCase.X, thisCase.Z};
  } 
  else if (plane==1) // Plane at X e.g. Side CRT West, East, Top CRT Vertical Rim West and East
  {
    CrossPoint thisCase=CalculateProjection(thisDir.dirx, thisDir.diry, thisDir.dirz, thisDir.meanx, thisDir.meany, thisDir.meanz, position);
    CrossingPoint={thisCase.X, thisCase.Y, thisCase.Z};
  } 
  else if (plane==2) // Plane at Z e.g. Side CRT South, North, Top CRT Vertical Rime South and North 
  {
    CrossPoint thisCase=CalculateProjection(thisDir.dirz, thisDir.diry, thisDir.dirx, thisDir.meanz, thisDir.meany, thisDir.meanx, position);
    CrossingPoint={thisCase.Z, thisCase.Y, thisCase.X};
  }
  return CrossingPoint;


}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
In this part the loop over the even is done with all the structures and functions declared before.
Cafana defines and object called SpillMultivar as an iteration along each event on a file from which an histogram gets filled with the variable of interest to be filled
The variables of interest for the tripleMatch are the CRTdistance and the Barycenter distance, since doing an histogram for each variable would be time consuming using CAFANA
what is done is a root file from which the variables get saved for further analysis 
*/

const SpillMultiVar kSliceCRTReconstuctorTree([](const caf::SRSpillProxy* sr)-> std::vector<double>  
{ 

   int evtnumber = sr->hdr.evt; //Event number stored in the flatCaf
   int runnumber = sr->hdr.run; //run number stored in the flatcaf

  /*
 As a first step all of the variables related to the charge and flash barycenters, along the TPC track and CRT hit coordinates are declared and
 filled with their corresponed value for an event
  */
  // The Charge variables
  double average_x_charge=0;
  double average_y_charge=0;
  double average_z_charge=0;
  // The flashes variables
  double average_opflash_x=0;
  double average_opflash_y=0;
  double barydistance=0;
  float average_opflash_z=0;
  float opflash_time=0;
  int opflash_cryo=0;
  int opflash_type=0;


  std::vector<double> flash_barycenters_z;
  std::vector<double> flash_barycenters_x;
  std::vector<double> flash_times;
  std::vector<int> flash_cryo; 
  std::vector<int> flash_type;

  //The vectors of vectors of the unmatched and matched CRThits (The matched come inside the code)
  std::vector<std::vector<CRThit>> AllCRThits; 
  std::vector<std::vector<CRThit>> MatchedCRThits; 

  //The vector of vectors that contains the candidates:
  std::vector<std::vector<std::vector<double>>> CRT_candidates;
  
  //The values to be returned for the CAFANA histogram 
  std::vector<double> vector_active;
  


  

  //Lets get all of the variables of interest for the CRTPMT-matches

  OpFlash ThisEventFlashes = CRTPMTValues(sr);

  //Lets get the opflash variables of interest:
  flash_barycenters_z= ThisEventFlashes.Barycenters_Z;
  flash_barycenters_x= ThisEventFlashes.Barycenters_X;
  flash_times = ThisEventFlashes.flash_times;
  flash_cryo = ThisEventFlashes.cryos;
  flash_type =ThisEventFlashes.flash_types;

  // And The CRT Hits for each flash
  MatchedCRThits =ThisEventFlashes.CRT;
  CRT_candidates = ThisEventFlashes.CRTPos;

  // Bools for the Matched Hits
  bool triggFlash= false;




 for(auto const& islc : sr->slc){  //Loop over the slices of the event to reconstruct all of the tracks present in the event. A slice correspond to a piece of event that was set by a trigger
   //lets get the barycenters
  average_z_charge= islc.charge_center.z;  //The barycenters  of charge of each slice are tendentially the same barycenters of the tracks. BaryZ gets the distance with the flash match. So all the tracks are tested with each charge center
  

  average_x_charge= islc.charge_center.x;  // The barycenter of charge in X is useful to determine if the reconstructed track and flash are in the same region of the TPC.

  
  
  int hitssizecounter=0;
  int totpfp=0;
  int npfp = islc.reco.npfp;
  for(std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp) //Loop over the members of the slice reconstruct all of the tracks present in the slice
   
  { 
     
    int tracklen = islc.reco.pfp[ipfp].trk.len;
    int plane=islc.reco.pfp[ipfp].trk.bestplane;
    int pfpSize=islc.reco.pfp[ipfp].trk.calo[2].points.size();
    int pfpHits=islc.reco.pfp[ipfp].trk.calo[2].nhit;
    
    if(islc.reco.pfp[ipfp].trk.calo[2].nhit== -999)  //Some hits are initialized at -99 so they mess with the loop over the side of the CaloHits, there fore set to 0
    {
      pfpHits=0;
    }

    if(abs(pfpSize -pfpHits)>10)
    {
      hitssizecounter++;
    }


    std::map<int,double> flash_goodCands; //Good cands takes the tracks for which the constraints are satisfied and stores the value of the distances
    
     if(plane !=0 || plane != 1 || plane != 2) //Same discourse that with nhits, sometimes the bestplane is not initialized correctly, so in this case is ponited to collection
     {
      plane =2;
     }
    if(islc.reco.pfp[ipfp].trk.calo[2].points.size()==0) continue; //Empty tracks are discarded. Most of the cosmic have no 

   
    std::vector<double> Track_charge_barycenter=ChargeBarycenter(islc.reco.pfp[ipfp]); //The barycenters are stored for each pfp

    // Aux Variables.
    double deltaXPCA=0;
    double deltaZPCA=0;
    double deltaXPANDORA=0;
    double deltaZPANDORA=0;

      
    

    for(size_t counter=0; counter < flash_barycenters_z.size();counter++) // Each member of the slice(PFP) that reconstructs a track is "matched" with each one of the  opFlashes contained in the event to check if they are part of the same signal.
    {
     average_opflash_z= flash_barycenters_z[counter]; //Barycenter of light at the current flash [counter]
     opflash_time = flash_times[counter]; // Time at the current flash [counter]
     opflash_cryo = flash_cryo[counter]; //In which cryostat was the flash measured( the info of the cryostat is then taken from the flash and not from the reconstructed track a the TPC)
     barydistance= abs(Track_charge_barycenter[2] - average_opflash_z- delta_bar_expected_from_mc(Track_charge_barycenter[2],opflash_cryo)); //Distance between the two mentioned barycenters, delta_bar_mc_expected corresponds to an empirical weight to consider also tracks far from the center of the TPC
     if(std::isnan(islc.reco.pfp[ipfp].trk.len)) continue; //Tracks without lenght are fake ones
     if(islc.reco.pfp[ipfp].trk.len <= 40) continue; //Tracks that are short are difficult to reconstruct, so a quality cut is set on 40cm
    
     if(islc.reco.pfp[ipfp].trk.start.x < 0 && opflash_cryo !=0) continue; //To make sure opflashes and TPC-reco are in the same region
     if(islc.reco.pfp[ipfp].trk.start.x > 0 && opflash_cryo !=1) continue; //To make sure opflashes and TPC-reco are in the same region
     if(barydistance >= 150) continue; //The conditional on the barycenter distance, if they are within 1 meter and a half they are considered as a candidate
     if(std::isnan(islc.reco.pfp[ipfp].t0) == false) 
     { 
       if(abs(islc.reco.pfp[ipfp].t0 - opflash_time )>10) continue; // Cathode crossers have a definite time, and from them is also useful to check if they occur within a reasonable window of time with respect to the flash       
     }

     //Here the function that compares the start and end of the trk with the flash time to see if its inside the FV
  
    DriftedTrack testTrack = Drift_tracks_outbound(islc.reco.pfp[ipfp], opflash_time, opflash_cryo);
    if(testTrack.outbound !=0) continue;  // if the outbound is greater than zero the track is clearly not matchable
    flash_goodCands.insert({counter,barydistance}); //the map of good candidates saves the number of the flash and the distance flash-charge for the tracks that are not outbound for the flash

       
    }
    //Now we go here to the CRT anaylsis, which can be done all inside the pfp in question to get the slice.

    std::vector<double> CRTsx, CRTsy, CRTsz;
    int ismatched=-1; 
    int ccmatched=0;
    if(flash_goodCands.size()>=1)
    { 
      std::map<int, double>::iterator itCands; //The track candidates are mapped with the distance from the flash barycenter for the crt projection
      std::map<int, std::pair<double,double>> bestTripleMatch; //Besttriplematch saves the distance of the projection and the distance of the barycenter distances
      for(itCands = flash_goodCands.begin();itCands != flash_goodCands.end();itCands++) //Loop over the good candidates to see their projection on the CRT 
      {
        int k= itCands->first; //the "number" of the flash in the event that passed the constraints
        double distancesBar= itCands->second; // This one takes de distance that passed the filter
        /*
        There are two directions one for the PCA reconstruction and the other for pandora
        */
        Direction thisDir;
        Direction thisDir2;
  
        DriftedTrack thisCand = Drift_tracks_outbound(islc.reco.pfp[ipfp], flash_times[k], flash_cryo[k]); //The outbound is checked for the good track candidates 
        
        thisDir = PCAfit(islc.reco.pfp[ipfp],flash_times[k], flash_cryo[k]); // The projections of the tracks are done or with the PCA reconstructed director cosines or with the Pandora ones; these are the PCA ones
        thisDir2= Pandorafit(islc.reco.pfp[ipfp],flash_times[k],flash_cryo[k]);//The projections of the tracks are done or with the PCA reconstructed director cosines or with the Pandora ones; these are the PANDORA ones



        double crtdistance= -99999;  //Arbitrary value for the crtdistance
        double crtdistancePandora= -999999;
        std::vector<double> CRTDistances; // To store the multiple values of the CRT Distances and save in the map just the shortest one.
        std::vector<double> CRTDistancesPandora; //FOR THE CRT COMPARISON PANDORA/PCA


        ////////// Here we loop over the CRTHITS found to be matched with the flash identified by k in CRTPMTMATCH, so we  acces that zone of the CAF and get the size. (crtpmt_matches.matchedCRTHits.size())
        for(size_t pp=0; pp<MatchedCRThits[k].size();pp++)
        {   
          double  projected_x_PCA=0;
          double  projected_y_PCA=0;
          double  projected_z_PCA=0;
           double  projected_x_PANDORA=0;
          double  projected_y_PANDORA=0;
          double  projected_z_PANDORA=0;
          double CRT_x=0;
          double CRT_z=0;
          //Here  the functions for the CRT projections:
          std::pair<int,double> CRTPlane = DeterminePlane(MatchedCRThits[k][pp]);
          CrossPoint crossPoint = DetermineProjection(thisDir, CRTPlane.first, CRTPlane.second);
          CrossPoint PandoraCross = DetermineProjection(thisDir2, CRTPlane.first, CRTPlane.second);

          if(CRTPlane.first==0)
          {
           projected_x_PANDORA = PandoraCross.X;
           projected_x_PCA= crossPoint.X;
           projected_z_PANDORA= PandoraCross.Z;
           projected_z_PCA =crossPoint.Z;
           CRT_x=MatchedCRThits[k][pp].x;
           CRT_z=MatchedCRThits[k][pp].z;
           deltaXPCA=MatchedCRThits[k][pp].x-crossPoint.X;
           deltaZPCA=MatchedCRThits[k][pp].z -crossPoint.Z;
           deltaXPANDORA=MatchedCRThits[k][pp].x-PandoraCross.X;
           deltaZPANDORA=MatchedCRThits[k][pp].z -PandoraCross.Z;
           crtdistance=sqrt(pow(MatchedCRThits[k][pp].x-crossPoint.X,2)+pow(MatchedCRThits[k][pp].z -crossPoint.Z,2));
           crtdistancePandora=sqrt(pow(MatchedCRThits[k][pp].x-PandoraCross.X,2)+pow(MatchedCRThits[k][pp].z -PandoraCross.Z,2));
        
           CRTDistances.push_back(crtdistance);
           CRTDistancesPandora.push_back(crtdistancePandora);

          

          }else if(CRTPlane.first==1)
          {
           projected_x_PANDORA = PandoraCross.Y;
           projected_x_PCA= crossPoint.Y;
           projected_z_PANDORA= PandoraCross.Z;
           projected_z_PCA =crossPoint.Z;
           CRT_x=MatchedCRThits[k][pp].y;
           CRT_z=MatchedCRThits[k][pp].z;
           deltaXPCA=MatchedCRThits[k][pp].y-crossPoint.Y;
           deltaZPCA=MatchedCRThits[k][pp].z -crossPoint.Z;
           deltaXPANDORA=MatchedCRThits[k][pp].y-PandoraCross.Y;
           deltaZPANDORA=MatchedCRThits[k][pp].z -PandoraCross.Z;
           crtdistance=sqrt(pow(MatchedCRThits[k][pp].y-crossPoint.Y,2)+pow(MatchedCRThits[k][pp].z -crossPoint.Z,2));
           crtdistancePandora=sqrt(pow(MatchedCRThits[k][pp].y-PandoraCross.Y,2)+pow(MatchedCRThits[k][pp].z -PandoraCross.Z,2));
           CRTDistances.push_back(crtdistance);
           CRTDistancesPandora.push_back(crtdistancePandora);
        
          
           
        

          }else if(CRTPlane.first==2)
          {
            projected_x_PANDORA = PandoraCross.X;
            projected_x_PCA= crossPoint.X;
            projected_z_PANDORA= PandoraCross.Y;
            projected_z_PCA =crossPoint.Y;
            CRT_x=MatchedCRThits[k][pp].x;
            CRT_z=MatchedCRThits[k][pp].y;
           deltaXPCA=MatchedCRThits[k][pp].x-crossPoint.X;
           deltaZPCA=MatchedCRThits[k][pp].y -crossPoint.Y;
           deltaXPANDORA=MatchedCRThits[k][pp].x-PandoraCross.X;
           deltaZPANDORA=MatchedCRThits[k][pp].x -PandoraCross.Y;
            crtdistance=sqrt(pow(MatchedCRThits[k][pp].x-crossPoint.X,2)+pow(MatchedCRThits[k][pp].y -crossPoint.Y,2));
            crtdistancePandora=sqrt(pow(MatchedCRThits[k][pp].x-PandoraCross.X,2)+pow(MatchedCRThits[k][pp].y -PandoraCross.Y,2));
            CRTDistances.push_back(crtdistance);
            CRTDistancesPandora.push_back(crtdistancePandora);
        

           
          }
            
            int Multiplicty= CRTDistances.size();
            //Here all of the variables of interest are written into a text file that after gets converted into a root file
            MyFile<<evtnumber<<" "<<npfp<<" "<<flash_times[k]<<" "<<MatchedCRThits[k][pp].t1-flash_times[k]<<" "<<flash_cryo[k]<<" "<<Multiplicty<<" "<<CRTPlane.first<<" "<<tracklen<<" "<<distancesBar<<" "<<crtdistance<<" "<<crtdistancePandora<<"  "<< deltaXPCA<<" "<<projected_x_PCA<<" "<<deltaXPANDORA<<" "<<projected_x_PANDORA<<" "<<deltaZPCA<<" "<<projected_z_PCA<<" "<<deltaZPANDORA<<" "<<projected_z_PANDORA<<" "<<thisDir.dirx<<" "<<thisDir.diry<<" "<<thisDir.dirz<<" "<<thisDir2.dirx<<" "<<thisDir2.diry<<" "<<thisDir2.dirz<<" "<<CRT_x<<" "<<CRT_z<<" "<<MatchedCRThits[k][pp].t1<<" "<<MatchedCRThits[k][pp].region<<" "<<thisDir.meanx<<" "<<thisDir.meany<<" "<<thisDir.meanz<<" "<<thisDir2.meanx<<" "<<thisDir2.meany<<" "<<thisDir2.meanz<<" "<<runnumber<<std::endl;
            /////<<<1<<" "  <<      2   <<" "<<      3       <<" "<<                   4                    <<" "<<     5      <<" "<<      6     <<" "<<      7      <<" "<<     8  <<" "<<     9      <<" "<<    10     <<" "<<        11        <<"  "<<     12   <<" "<<       13      <<" "<<      14     <<" "<<       15          <<" "<<    16   <<" "<<      17       <<" "<<     18      <<" "<<           19      <<" "<<      20    <<" "<<     21     <<" "<<      22    <<" "<<     23      <<" "<<       24    <<" "<<     25      <<" "<<  26 <<" "<< 27  <<" "<<          28            <<" "<<               29           <<std::endl;

          
           if(CRTDistances.size()>1) //Getting the minimal CRT-Distance for the matched flash in case it has two(or more) CRT Vectors in it (a map only saves a value for a key, indicatively the last one)
           { 
             int iMin=0;
            for(size_t aa=1; aa<CRTDistances.size();aa++)
             {
              if(CRTDistances[iMin] > CRTDistances[aa])
              {
                iMin=aa;
              }
                    

             }
             crtdistance = CRTDistances[iMin];
           }
            
          
           

          
        }
        std::pair<double,double> thisCandDist = {itCands->second,crtdistance}; //The pair that adds Barydistances and the projection-crtdistances 
        bestTripleMatch.insert({k,thisCandDist}); //Here i have the set of distances from barycenter and crt-dist for each flash number of the event in this slice
             

       
      }
     auto BestCand = std::min_element(bestTripleMatch.begin(), bestTripleMatch.end(), []( const std::pair<int, std::pair<double,double>> &x, const std::pair<int,std::pair<double,double>> &y)   
     { 
      return (x.second.first+x.second.second) < (y.second.first+y.second.second);
     });
     //For the first histogram the spectrum saves the "best CRT distance in the detector"
     vector_active.push_back(BestCand->second.second); // We find the distance of the "best-candidate" and save it into the spectrum after the selection

      flash_goodCands.clear(); 
      bestTripleMatch.clear();
    
  }
  }
  }

return vector_active;
MyFile.close();
});


///////////////////////////////////////////////////////////////////////////////////////////// THE SECOND LOOP//////////////////////////////////////
//CAFANA IS NOT VERY FRIENDLY WHEN IT COMES TO  SET MULTIVARIATE SPECTRUMS, SO BY ITSELF IT ONLY STORES ONE-DIMENSIONAL SPECTRUMS THAT CAN BE COMBINED INTO TWO
//SO FOR EXAMPLE IF I WANTED AN SPECTRUM FOR EACH VARIABLE I USED IN THE ROOT I WOULD NEED TO DO AT LEAST 10 TIMES THE SAME ALGORITHM (THERE'S A REASON WHY I'M USING A TREE)
// IN ANY CASE I WILL USE ANOTHER SPECTRUM TO GET THE BEST BARYCENTER DISTANCE AND TWO A 2 DIMENSIONAL SPECTRUM FOR THE SAKE OF DEMONSTRATION, THE CODE IS THE SAME OF BEFORE, JUST WITH A DIFFERENT PUSH BACK ON THE END VECTOR (AND NO WRITE FILE)
 const SpillMultiVar kSliceCRTReconstuctorBary([](const caf::SRSpillProxy* sr)-> std::vector<double>  
{ 

   int evtnumber = sr->hdr.evt;
   int runnumber = sr->hdr.run;


  /*
 As a first stepp all of the variables related to the charge and flash barycenters, along the TPC track and CRT hit coordinates are declared and
 filled with their corresponed value for an event
  */
  // The Charge variables
  double average_x_charge=0;
  double average_y_charge=0;
  double average_z_charge=0;
  // The flashes variables
  double average_opflash_x=0;
  double average_opflash_y=0;
  double barydistance=0;
  float average_opflash_z=0;
  float opflash_time=0;
  int opflash_cryo=0;
  int opflash_type=0;


  std::vector<double> flash_barycenters_z;
  std::vector<double> flash_barycenters_x;
  std::vector<double> flash_times;
  std::vector<int> flash_cryo; 
  std::vector<int> flash_type;

  //The vectors of vectors of the unmatched and matched CRThits (The matched come inside the code)
  std::vector<std::vector<CRThit>> AllCRThits; 
  std::vector<std::vector<CRThit>> MatchedCRThits; 

  //The vector of vectors that contains the candidates:
  std::vector<std::vector<std::vector<double>>> CRT_candidates;
  
  //The values to be returned for the CAFANA histogram 
  std::vector<double> vector_active;
  
  //Test for the plot of PCA and Pandora


  
  //AllTripleMatchAnaVars.evtnumber = sr->hdr.evt;

  //Lets get all of the variables of interest for the CRTPMT-matches

  OpFlash ThisEventFlashes = CRTPMTValues(sr);

  //Lets get the opflash variables of interest:
  flash_barycenters_z= ThisEventFlashes.Barycenters_Z;
  flash_barycenters_x= ThisEventFlashes.Barycenters_X;
  flash_times = ThisEventFlashes.flash_times;
  flash_cryo = ThisEventFlashes.cryos;
  flash_type =ThisEventFlashes.flash_types;

  // And The CRT Hits for each flash
  MatchedCRThits =ThisEventFlashes.CRT;
  CRT_candidates = ThisEventFlashes.CRTPos;

  // Bools for the Matched Hits
  bool triggFlash= false;




  
 for(auto const& islc : sr->slc){  //Loop over the slices of the event to reconstruct all of the tracks present in the event. A slice correspond to a piece of event that was set by a trigger
   //lets get the barycenters
  average_z_charge= islc.charge_center.z;  //The barycenters  of charge of each slice are tendentially the same barycenters of the tracks. BaryZ gets the distance with the flash match. So all the tracks are tested with each charge center
  

  average_x_charge= islc.charge_center.x;  // The barycenter of charge in X is useful to determine if the reconstructed track and flash are in the same region of the TPC.

  
  
  int hitssizecounter=0;
  int totpfp=0;
  int npfp = islc.reco.npfp;

  for(std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp) //Loop over the members of the slice reconstruct all of the tracks present in the slice
   
  { 
     
    int tracklen = islc.reco.pfp[ipfp].trk.len;
    int plane=islc.reco.pfp[ipfp].trk.bestplane;
    int pfpSize=islc.reco.pfp[ipfp].trk.calo[2].points.size();
    int pfpHits=islc.reco.pfp[ipfp].trk.calo[2].nhit;
    
    if(islc.reco.pfp[ipfp].trk.calo[2].nhit== -999)  //Some hits are initialized at -99 so they mess with the loop over the side of the CaloHits, there fore set to 0
    {
      pfpHits=0;
    }

    if(abs(pfpSize -pfpHits)>10)
    {
      hitssizecounter++;
    }


    std::map<int,double> flash_goodCands; //Good cands takes the tracks for which the constraints are satisfied and stores the value of the distances
    
     if(plane !=0 || plane != 1 || plane != 2) //Same discourse that with nhits, sometimes the bestplane is not initialized correctly, so in this case is ponited to collection
     {
      plane =2;
     }
    if(islc.reco.pfp[ipfp].trk.calo[2].points.size()==0) continue; //Empty tracks are discarded. Most of the cosmic have no 


   
    std::vector<double> Track_charge_barycenter=ChargeBarycenter(islc.reco.pfp[ipfp]); //The barycenters are stored for each pfp

    // Aux Variables.
    double deltaXPCA=0;
    double deltaZPCA=0;
    double deltaXPANDORA=0;
    double deltaZPANDORA=0;

      
    

    for(size_t counter=0; counter < flash_barycenters_z.size();counter++) // Each member of the slice(PFP) that reconstructs a track is "matched" with each one of the  opFlashes contained in the event to check if they are part of the same signal.
    {
     average_opflash_z= flash_barycenters_z[counter]; //Barycenter of light at the current flash [counter]
     opflash_time = flash_times[counter]; // Time at the current flash [counter]
     opflash_cryo = flash_cryo[counter]; //In which cryostat was the flash measured( the info of the cryostat is then taken from the flash and not from the reconstructed track a the TPC)
     barydistance= abs(Track_charge_barycenter[2] - average_opflash_z- delta_bar_expected_from_mc(Track_charge_barycenter[2],opflash_cryo)); //Distance between the two mentioned barycenters, delta_bar_mc_expected corresponds to an empirical weight to consider also tracks far from the center of the TPC
     if(std::isnan(islc.reco.pfp[ipfp].trk.len)) continue; //Tracks without lenght are fake ones
     if(islc.reco.pfp[ipfp].trk.len <= 40) continue; //Tracks that are short are difficult to reconstruct, so a quality cut is set on 40cm
    
     if(islc.reco.pfp[ipfp].trk.start.x < 0 && opflash_cryo !=0) continue; //To make sure opflashes and TPC-reco are in the same region
     if(islc.reco.pfp[ipfp].trk.start.x > 0 && opflash_cryo !=1) continue; //To make sure opflashes and TPC-reco are in the same region
     if(barydistance >= 150) continue; //The conditional on the barycenter distance, if they are within 1 meter and a half they are considered as a candidate
     if(std::isnan(islc.reco.pfp[ipfp].t0) == false) 
     { 
       if(abs(islc.reco.pfp[ipfp].t0 - opflash_time )>10) continue; // Cathode crossers have a definite time, and from them is also useful to check if they occur within a reasonable window of time with respect to the flash       
     }

     //Here the function that compares the start and end of the trk with the flash time to see if its inside the FV
  
    DriftedTrack testTrack = Drift_tracks_outbound(islc.reco.pfp[ipfp], opflash_time, opflash_cryo);
    if(testTrack.outbound !=0) continue;  // if the outbound is greater than zero the track is clearly not matchable
    flash_goodCands.insert({counter,barydistance}); //the map of good candidates saves the number of the flash and the distance flash-charge for the tracks that are not outbound for the flash

       
    }
    //Now we go here to the CRT anaylsis, which can be done all inside the pfp in question to get the slice.

    std::vector<double> CRTsx, CRTsy, CRTsz;
    int ismatched=-1; 
    int ccmatched=0;
    if(flash_goodCands.size()>=1)
    { 
      std::map<int, double>::iterator itCands; //The track candidates are mapped with the distance from the flash barycenter for the crt projection
      std::map<int, std::pair<double,double>> bestTripleMatch; //Besttriplematch saves the distance of the projection and the distance of the barycenter distances
      for(itCands = flash_goodCands.begin();itCands != flash_goodCands.end();itCands++) //Loop over the good candidates to see their projection on the CRT 
      {
        int k= itCands->first; //the "number" of the flash in the event that passed the constraints
        double distancesBar= itCands->second; // This one takes de distance that passed the filter
        /*
        There are two directions one for the PCA reconstruction and the other for pandora
        */
        Direction thisDir;
        Direction thisDir2;
  
        DriftedTrack thisCand = Drift_tracks_outbound(islc.reco.pfp[ipfp], flash_times[k], flash_cryo[k]); //The outbound is checked for the good track candidates 
        
        thisDir = PCAfit(islc.reco.pfp[ipfp],flash_times[k], flash_cryo[k]); // The projections of the tracks are done or with the PCA reconstructed director cosines or with the Pandora ones; these are the PCA ones
        thisDir2= Pandorafit(islc.reco.pfp[ipfp],flash_times[k],flash_cryo[k]);//The projections of the tracks are done or with the PCA reconstructed director cosines or with the Pandora ones; these are the PANDORA ones



        double crtdistance= -99999;  //Arbitrary value for the crtdistance
        double crtdistancePandora= -999999;
        std::vector<double> CRTDistances; // To store the multiple values of the CRT Distances and save in the map just the shortest one.
        std::vector<double> CRTDistancesPandora; //FOR THE CRT COMPARISON PANDORA/PCA


        ////////// Here we loop over the CRTHITS found to be matched with the flash identified by k in CRTPMTMATCH, so we  acces that zone of the CAF and get the size. (crtpmt_matches.matchedCRTHits.size())
        for(size_t pp=0; pp<MatchedCRThits[k].size();pp++)
        {   
          /*
          Here we have a set of projected coordinates onto the CRT plane, done with the angles projected by pandora and the PCA to check their efficiencies 
          Spoiler: The PCA performs better
          */
          double  projected_x_PCA=0; 
          double  projected_y_PCA=0;
          double  projected_z_PCA=0;
           double  projected_x_PANDORA=0;
          double  projected_y_PANDORA=0;
          double  projected_z_PANDORA=0;
          double CRT_x=0;
          double CRT_z=0;
          //Here  the functions for the CRT projections:
          std::pair<int,double> CRTPlane = DeterminePlane(MatchedCRThits[k][pp]); // finds the plane of the CRT 
          CrossPoint crossPoint = DetermineProjection(thisDir, CRTPlane.first, CRTPlane.second); //Projects the track onto the plane
          CrossPoint PandoraCross = DetermineProjection(thisDir2, CRTPlane.first, CRTPlane.second); //Projects the track onto the plane using pandora angles

          if(CRTPlane.first==0)
          {
           projected_x_PANDORA = PandoraCross.X;
           projected_x_PCA= crossPoint.X;
           projected_z_PANDORA= PandoraCross.Z;
           projected_z_PCA =crossPoint.Z;
           CRT_x=MatchedCRThits[k][pp].x;
           CRT_z=MatchedCRThits[k][pp].z;
           deltaXPCA=MatchedCRThits[k][pp].x-crossPoint.X;
           deltaZPCA=MatchedCRThits[k][pp].z -crossPoint.Z;
           deltaXPANDORA=MatchedCRThits[k][pp].x-PandoraCross.X;
           deltaZPANDORA=MatchedCRThits[k][pp].z -PandoraCross.Z;
           crtdistance=sqrt(pow(MatchedCRThits[k][pp].x-crossPoint.X,2)+pow(MatchedCRThits[k][pp].z -crossPoint.Z,2));
           crtdistancePandora=sqrt(pow(MatchedCRThits[k][pp].x-PandoraCross.X,2)+pow(MatchedCRThits[k][pp].z -PandoraCross.Z,2));
           CRTDistances.push_back(crtdistance);
           CRTDistancesPandora.push_back(crtdistancePandora);
      
            
          
          }else if(CRTPlane.first==1)
          {
           projected_x_PANDORA = PandoraCross.Y;
           projected_x_PCA= crossPoint.Y;
           projected_z_PANDORA= PandoraCross.Z;
           projected_z_PCA =crossPoint.Z;
           CRT_x=MatchedCRThits[k][pp].y;
           CRT_z=MatchedCRThits[k][pp].z;
           deltaXPCA=MatchedCRThits[k][pp].y-crossPoint.Y;
           deltaZPCA=MatchedCRThits[k][pp].z -crossPoint.Z;
           deltaXPANDORA=MatchedCRThits[k][pp].y-PandoraCross.Y;
           deltaZPANDORA=MatchedCRThits[k][pp].z -PandoraCross.Z;
           crtdistance=sqrt(pow(MatchedCRThits[k][pp].y-crossPoint.Y,2)+pow(MatchedCRThits[k][pp].z -crossPoint.Z,2));
           crtdistancePandora=sqrt(pow(MatchedCRThits[k][pp].y-PandoraCross.Y,2)+pow(MatchedCRThits[k][pp].z -PandoraCross.Z,2));
          
           CRTDistancesPandora.push_back(crtdistancePandora);
         

           
        

          }else if(CRTPlane.first==2)
          {
            projected_x_PANDORA = PandoraCross.X;
            projected_x_PCA= crossPoint.X;
            projected_z_PANDORA= PandoraCross.Y;
            projected_z_PCA =crossPoint.Y;
            CRT_x=MatchedCRThits[k][pp].x;
            CRT_z=MatchedCRThits[k][pp].y;
           deltaXPCA=MatchedCRThits[k][pp].x-crossPoint.X;
           deltaZPCA=MatchedCRThits[k][pp].y -crossPoint.Y;
           deltaXPANDORA=MatchedCRThits[k][pp].x-PandoraCross.X;
           deltaZPANDORA=MatchedCRThits[k][pp].x -PandoraCross.Y;
            crtdistance=sqrt(pow(MatchedCRThits[k][pp].x-crossPoint.X,2)+pow(MatchedCRThits[k][pp].y -crossPoint.Y,2));
            crtdistancePandora=sqrt(pow(MatchedCRThits[k][pp].x-PandoraCross.X,2)+pow(MatchedCRThits[k][pp].y -PandoraCross.Y,2));
           
            CRTDistancesPandora.push_back(crtdistancePandora);

           
          }
            
            int Multiplicty= CRTDistances.size();
            

          
           if(CRTDistances.size()>1) //Getting the minimal CRT-Distance for the matched flash in case it has two(or more) CRT Vectors in it (a map only saves a value for a key, indicatively the last one)
           { 
             int iMin=0;
            for(size_t aa=1; aa<CRTDistances.size();aa++)
             {
              if(CRTDistances[iMin] > CRTDistances[aa])
              {
                iMin=aa;
              }
                    

             }
             crtdistance = CRTDistances[iMin];
           }
          
           

          
        }
        std::pair<double,double> thisCandDist = {itCands->second,crtdistance}; //The pair that adds Barydistances and the projection-crtdistances 
        bestTripleMatch.insert({k,thisCandDist}); //Here i have the set of distances from barycenter and crt-dist for each flash number of the event in this slice
             

       
      }
     auto BestCand = std::min_element(bestTripleMatch.begin(), bestTripleMatch.end(), []( const std::pair<int, std::pair<double,double>> &x, const std::pair<int,std::pair<double,double>> &y)   
     { 
      return (x.second.first+x.second.second) < (y.second.first+y.second.second);
     });
     vector_active.push_back(BestCand->second.first); // We find the barycenter distance of the "best-candidate" and save it into the spectrum after the selection

      flash_goodCands.clear();  //Clear vectors to not save info not needed
      bestTripleMatch.clear();  //Clear vectors to not save info not needed
    
  }
  }
  
  }

return vector_active;


});




