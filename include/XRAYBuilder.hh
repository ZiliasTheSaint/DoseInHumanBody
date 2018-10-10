
//D. Fulea, National Institute of Public Health Bucharest, Cluj-Napoca Regional Center, Romania

#ifndef XRayBuilder_h
#define XRayBuilder_h 1

#include <iostream>
#include <fstream>
#include <vector>

#include "globals.hh"//for G4cout test!!
using namespace std;

class XRayBuilder
{
public:

  XRayBuilder(double kv, double filtration, double anodAngle, double ripple, string anodMaterial);
  ~XRayBuilder();
  //=============  
  int GetBins(){return NEFF;}; 
  double* GetEnergyArray(){
	  double* meven=new double[NEFF];
	  for (int i=0;i<NEFF;i++)
		  meven[i]=en_array[i]/1000.0;//toMeV
	  return meven;//en_array;
  };  
  
  double* GetProbArray(){
	  double* p=new double[NEFF];
	  for (int i=0;i<NEFF;i++)
		  p[i]=YF[i]/TOT_PHOTONS;
	  return p;
  };

  double Get_photons_per_mm2_per_uGy(){
	  return TOT_PHOTONS/KERMAPERMASAT750MM;
  }

  double GetAirKerma_per_mAs_at75cm(){
	  return KERMAPERMASAT750MM;
  }
  //==============
  static string filename;
  static const string dataS;
  static const string data_filterS;
  static const string data_kermaS;
  static const string data_wS;
  static const string data_mS;
  static const string data_rS;
  static const string data_constantS;
  static const string data_rippleS;
  static const string defaultSpectrumExt;
  static const string filterKermaExt;
  static const string rippleExt_R0;
  static const string rippleExt_R;
  //================
  

  private:
	  double kv;
      double filtration;
      double anodAngle;
      double ripple;
      string anodMaterial;
	  //---------
	  int ndata;
	  int NEFF;//store the effective number of energies..must be equal to ndata for corect data!
	  double* en_array;//store the energy array in keV
	  double* YF;//store attenuated photon flux [photons/(mAs*mm2) at 75 cm] array.
	  double* XVAL;//store the [photons/(mAs*mm2) at 750 mm] read from SPECTRUM (.SPC) file.
	  //---------
	  int ianod;//index for W, MO or RH
	  double DE;//delta E=>energy BIN width according to spectrum file!
	  double* KVAL;//photonsToKerma conversion factors array [uGy*mm2/photons] at 750 mm
	  //-----------
	  double* YS;//when spectrum is read=> we get: En [keV] and YS=XVAL, which is photons/(mAs*mm2) at 75 cm!!!
	  double* TMM;//filtration or thickness of each attenuator [mm].
	  int NOATT;//number of attenuators
	  double** ATT;//matrix for linear attenuation coefficients of multiple attenuators (number of attenuators=NOATT). mm-1.
	  double TOT_PHOTONS;//store total photon flux = SUM of YF; used in normalization
	  double* ATTVAL;//used in single attenuator case
	  int MAXITER;// = 10000;
	  double THI_HVL1;// = 120.0;//Not possible to have more than 120 mmAl HVL. Set the upper bound for iterations
	  double MEANSPECTRUMENERGY;
	  double KERMAPERMASAT750MM;
	  double HVL1;
	  double HVL1_TISS;
	  //---------
	  void generateSpectrum();
	  void readSpectrum(string filename);
	  void readKerma(string filename);
	  void readAttCoef(string filename, int j);
	  void buildSpectra();
	  void computeHVL1(string filename, bool tissB);
	  void resetYS();
	  void readAttCoef(string filename);
	  //-------------------------------
	  static double stringToDouble(string s){
		  double dbl=0.0;
		  istringstream buffer(s);
		  buffer>>dbl;
		  return dbl;
	  };
	  static string doubleToString(double d){
		  ostringstream strs;
		  strs << d;
          string s = strs.str();
		  return s;
	  };
	  static string intToString(int d){
		  ostringstream strs;
		  strs << d;
          string s = strs.str();
		  return s;
	  };
};
#endif
