//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//D. Fulea, National Institute of Public Health Bucharest, Cluj-Napoca Regional Center, Romania
//
#ifndef HumanPhantomPrimaryGeneratorAction_h
#define HumanPhantomPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include <vector>

class HumanPhantomPrimaryGeneratorMessenger;
class G4ParticleGun;
class G4Event;
class HumanPhantomConstruction;

class HumanPhantomPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    HumanPhantomPrimaryGeneratorAction();
    ~HumanPhantomPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);
    void GenerateBeamAlongZ();
    void GenerateBeamAlongX();
    void GenerateBeamAlongY();
    void GenerateIsotropicFlux();
    void SetBeam(G4String);
	void GenerateRectangleField();//like X-ray radiography
	void GenerateCTScan();//like a standard CT fan beam
	void GenerateMammoField();

	G4ParticleGun* GetParticleGun() const { return particleGun; }

	static void SetPrintModulo(G4int    val)  {printModulo = val;};
	static G4int     printModulo;

	static void SetFieldCenterY(G4double val){HumanPhantomPrimaryGeneratorAction::fieldCenterY=val;};
	static G4double GetFieldCenterY();//{return HumanPhantomPrimaryGeneratorAction::fieldCenterY;};
	static G4double fieldCenterY;

	static void SetFieldCenterX(G4double val){HumanPhantomPrimaryGeneratorAction::fieldCenterX=val;};
	static G4double GetFieldCenterX();//{return HumanPhantomPrimaryGeneratorAction::fieldCenterY;};
	static G4double fieldCenterX;

	static void SetFieldHeight(G4double val){HumanPhantomPrimaryGeneratorAction::fieldHeight=val;};
	static G4double GetFieldHeight();
	static G4double fieldHeight;

	static void SetFieldWidth(G4double val){HumanPhantomPrimaryGeneratorAction::fieldWidth=val;};
	static G4double GetFieldWidth();
	static G4double fieldWidth;

	static void SetAzimuthalRotationAngle(G4double val){HumanPhantomPrimaryGeneratorAction::azimuthalRotationAngle=val;};
	static G4double GetAzimuthalRotationAngle();
	static G4double azimuthalRotationAngle;

	static void SetPolarRotationAngle(G4double val){HumanPhantomPrimaryGeneratorAction::polarRotationAngle=val;};
	static G4double GetPolarRotationAngle();
	static G4double polarRotationAngle;

	//===============
	static void SetIfSpectrum(G4String val){
		if (val=="yes")
		 HumanPhantomPrimaryGeneratorAction::isSpectrum=true;
		else
		 HumanPhantomPrimaryGeneratorAction::isSpectrum=false;
	};
	static G4bool IsSpectrum();
	static G4bool isSpectrum;
	static double* ws_array;
	static int* ibin_array;
	static double* en_array;
	static int ndata;
	static void SetSpectrumAliasSamplingData(double* ws_a, int* ibin_a, double*en_a, int ndat){
		ws_array=ws_a;
		ibin_array=ibin_a;
		en_array=en_a;
		ndata=ndat;
	}
	static G4double incidentEnergy;
	static void SetDAP_uGymm2(std::string val)
	 {std::istringstream buffer(val); buffer>>HumanPhantomPrimaryGeneratorAction::DAP_uGymm2;};
	static double GetDAP_uGymm2();
	static double DAP_uGymm2;
	static void SetFocusToPhantomCentralAxisDistance(G4double val){HumanPhantomPrimaryGeneratorAction::FCA=val;};
	static G4double GetFCA();
	static G4double FCA;
	static void SetCTDI(std::string val)
	 {std::istringstream buffer(val); buffer>>HumanPhantomPrimaryGeneratorAction::CTDI_uGy;};
	static double GetCTDI();
	static double CTDI_uGy;
	//in CTCASE we must have CTDI as dose per rotation INTDz/nslices*slicethickness->e.g. [INT(-50mm;50mm)Dose]/100; 100 slices of 1 mm thick.
	//that CTDI must be multiplyed with according area to obtain a DAP version of CTDI. That ARea= fieldWidth*sliceThickness
	//===============@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	static void SetMAS(std::string val)
	 {std::istringstream buffer(val); buffer>>HumanPhantomPrimaryGeneratorAction::MAS;};
	static double GetMAS();
	static double MAS;//miliAmpereSeconds
	static void SetIfUseMASforDAP(G4String val){
		if (val=="yes")
		 HumanPhantomPrimaryGeneratorAction::useMASforDAP=true;
		else
		 HumanPhantomPrimaryGeneratorAction::useMASforDAP=false;
	};	
	static G4bool useMASforDAP;

	static double kv;
	static void SetKv(std::string val){std::istringstream buffer(val); buffer>>kv;};
    static double filtration;
	static void SetFiltration(std::string val){std::istringstream buffer(val); buffer>>filtration;};
    static double anodAngle;
	static void SetAnodAngle(std::string val){std::istringstream buffer(val); buffer>>anodAngle;};
    static double ripple;
	static void SetRipple(std::string val){std::istringstream buffer(val); buffer>>ripple;};
	static std::string anodMaterial;
	static void SetAnodMaterial(std::string val){anodMaterial=val;};
	//================@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	static void SetNumberOfEvents(G4int val){HumanPhantomPrimaryGeneratorAction::nbOfEvents=val;};	
	static G4int nbOfEvents;
	
	static void SetSliceThickness(G4double val){HumanPhantomPrimaryGeneratorAction::sliceThickness=val;};
	static G4double GetSliceThickness();
	static G4double sliceThickness;

	static void SetAngleIncrement(G4double val){HumanPhantomPrimaryGeneratorAction::angleIncrement=val;};
	static G4double GetAngleIncrement();
	static G4double angleIncrement;

	static void SetFanBeam(G4String val){
		if (val=="on")
		 HumanPhantomPrimaryGeneratorAction::isFanBeam=true;
		else
		 HumanPhantomPrimaryGeneratorAction::isFanBeam=false;
	};

	static void SetPitch(std::string val){std::istringstream buffer(val); buffer>>HumanPhantomPrimaryGeneratorAction::pitch;};
	//static double GetPitch(){return pitch;};//not used, static access!!
	static double pitch;

	static G4bool GetFanBeam();
	static G4bool isFanBeam;

	static G4bool GetHelicalScan();
	static void SetHelicalScan(G4String val){
		if (val=="on")
		 HumanPhantomPrimaryGeneratorAction::helicalScan=true;
		else
		 HumanPhantomPrimaryGeneratorAction::helicalScan=false;
	};
	static G4bool helicalScan;
	
	static G4bool GetHalfField();
	static void SetHalfField(G4String val){
		if (val=="on")
		 HumanPhantomPrimaryGeneratorAction::halfField=true;
		else
		 HumanPhantomPrimaryGeneratorAction::halfField=false;
	};
	static G4bool halfField;

	static G4bool GetDentalPanoramic();
	static void SetDentalPanoramic(G4String val){
		if (val=="on")
		 HumanPhantomPrimaryGeneratorAction::dentalPanoramic=true;
		else
		 HumanPhantomPrimaryGeneratorAction::dentalPanoramic=false;
	};
	static G4bool dentalPanoramic;


	static G4bool isCTScan;
	static G4int nbOfEventsPerScan;
	static G4int GetNumberOfEventsPerScan();
	static G4int minimumNbOfEvents;
	static G4bool isNumberOfEventsValid(){//called from run action at begining of run!
		if (!isCTScan)
			return true;
		//if here we must check if total number of events can generate at least 1 number of events per a single scan.
		//G4double DELTA = sliceThickness*(pitch-1.0);//usually 0!!
		//t+t+...=>l+l+l....;l=p x t
		
		if (dentalPanoramic){
		nbOfEventsPerScan=nbOfEvents*angleIncrement/(180.0*deg);
		minimumNbOfEvents=180.0*deg/(angleIncrement);
		//int division is always rounded down!!
		if (nbOfEventsPerScan>=1) return true;
		else return false;
		
		}else{
		nbOfEventsPerScan=nbOfEvents*pitch*sliceThickness*angleIncrement/(360.0*deg*fieldHeight);
		minimumNbOfEvents=fieldHeight*360.0*deg/(pitch*sliceThickness*angleIncrement);
		//int division is always rounded down!!
		if (nbOfEventsPerScan>=1) return true;
		else return false;
		}
	};

	static int rotations;
	static double GetRotations(){
		return rotations/360.0;
	};

	static int helical_rotations;
	static double GetHelicalRotations(){
		return helical_rotations/360.0;
	};
	//static G4double scaleXY;
	//static G4double scaleZ;

	static G4int IACTIVITY;//take number of quantas from source activity
	static G4int IKERMA;//take number of quantas from dosimetric measurements
	static G4int INUMBER;//take number of quantas from direct measurements using a particle counter
	static G4int IEVENTS;//take number of quantas directly from simulation number of runs
	static G4int IQUANTA;//this is choosen from above values; default is 3 which is IEVENTS
	static void SetIQUANTA(G4int    val)  {IQUANTA = val;};
	static G4int GetIQUANTA();

	//emmission quanta source can be isotropic or not (have a confined emmission area)
	static G4bool GetIsotropicSource();
	static void SetIsotropicSource(G4String val){
		if (val=="on")
		 HumanPhantomPrimaryGeneratorAction::isIsotropicSource=true;
		else
		 HumanPhantomPrimaryGeneratorAction::isIsotropicSource=false;
	};
	static G4bool isIsotropicSource;

	static double emissionArea_mm2;//emmission area in mm2. FCA,w,h from getter are all in mm!!!
	static void SetEmmissionArea_mm2(std::string val) 
	{std::istringstream buffer(val); buffer>>HumanPhantomPrimaryGeneratorAction::emissionArea_mm2;};
	static double GetEmmissionArea_mm2();

	static double activity_Bq;//source activity in Bq
	static void SetActivity_Bq(std::string val) 
	{std::istringstream buffer(val); buffer>>HumanPhantomPrimaryGeneratorAction::activity_Bq;};
	static double GetActivity_Bq();

	static double exposureTime_s;//exposure time in seconds
	static void SetExposureTime_s(std::string val) 
	{std::istringstream buffer(val); buffer>>HumanPhantomPrimaryGeneratorAction::exposureTime_s;};
	static double GetExposureTime_s();

	static double kerma_Gy;//kerma in Gy
	static void SetKerma_Gy(std::string val) 
	{std::istringstream buffer(val); buffer>>HumanPhantomPrimaryGeneratorAction::kerma_Gy;};
	static double GetKerma_Gy();

	static double massEnergyTransferCoefficient_cm2Perg;//mass energy transfer coefficient
	//for given energy taken from NIST tables or from RadiationHelper, g module
	//if radiation consists of photons!!!!!
	static void SetMassEnergyTransferCoefficient_cm2Perg(std::string val) 
	{std::istringstream buffer(val); buffer>>HumanPhantomPrimaryGeneratorAction::massEnergyTransferCoefficient_cm2Perg;};
	static double GetMassEnergyTransferCoefficient_cm2Perg();

	static double particleCounter;//particle counter as a number
	static void SetParticleCounter(std::string val) 
	{std::istringstream buffer(val); buffer>>HumanPhantomPrimaryGeneratorAction::particleCounter;};
	static double GetParticleCounter();

	static double radiationYield;//emmission probability
	static void SetRadiationYield(std::string val) 
	{std::istringstream buffer(val); buffer>>HumanPhantomPrimaryGeneratorAction::radiationYield;};
	static double GetRadiationYield();

  private:
    G4ParticleGun* particleGun;

    G4double x0;
    G4double y0;
    G4double z0;
    
    std::vector<G4double> probability;
   
    HumanPhantomPrimaryGeneratorMessenger* messenger;
	HumanPhantomConstruction*    phantom;     //pointer to the geometry
    
    G4String beamKind;
    G4double worldLength;
	G4int sliceIndex;
	G4int angleIndex;
	
	G4int eventIndex;	
	
	G4double aliasSample();
};

#endif


