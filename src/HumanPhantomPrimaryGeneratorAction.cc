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
// Authors: D. Fulea, National Institute of Public Health Bucharest, Cluj-Napoca Regional Center, Romania
//
#include "HumanPhantomPrimaryGeneratorAction.hh"
#include "HumanPhantomConstruction.hh"
#include "HumanPhantomPrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

#include "G4RunManager.hh"

#include "G4ios.hh"

G4double HumanPhantomPrimaryGeneratorAction::fieldCenterX=1.0;
G4double HumanPhantomPrimaryGeneratorAction::fieldCenterY=1.0;
G4double HumanPhantomPrimaryGeneratorAction::fieldHeight=1.0;
G4double HumanPhantomPrimaryGeneratorAction::fieldWidth=1.0;
G4double HumanPhantomPrimaryGeneratorAction::azimuthalRotationAngle=1.0;
G4double HumanPhantomPrimaryGeneratorAction::polarRotationAngle=0.0;
G4int HumanPhantomPrimaryGeneratorAction::printModulo=1;
G4int HumanPhantomPrimaryGeneratorAction::nbOfEvents=1;
G4bool HumanPhantomPrimaryGeneratorAction::isCTScan=false;
G4double HumanPhantomPrimaryGeneratorAction::sliceThickness=1.0*cm;
G4double HumanPhantomPrimaryGeneratorAction::angleIncrement=1.0*deg;
G4bool HumanPhantomPrimaryGeneratorAction::isFanBeam=true;
G4bool HumanPhantomPrimaryGeneratorAction::helicalScan=true;
G4bool HumanPhantomPrimaryGeneratorAction::halfField=false;
G4bool HumanPhantomPrimaryGeneratorAction::dentalPanoramic=false;
G4int HumanPhantomPrimaryGeneratorAction::nbOfEventsPerScan=1;
G4int HumanPhantomPrimaryGeneratorAction::minimumNbOfEvents=1;//for info
G4bool HumanPhantomPrimaryGeneratorAction::isSpectrum=false;
double* HumanPhantomPrimaryGeneratorAction::ws_array;//=new double[0];
int* HumanPhantomPrimaryGeneratorAction::ibin_array;//=new int[0];
double* HumanPhantomPrimaryGeneratorAction::en_array;
int HumanPhantomPrimaryGeneratorAction::ndata;
G4double HumanPhantomPrimaryGeneratorAction::incidentEnergy=0.100*MeV;
double HumanPhantomPrimaryGeneratorAction::kv=0.0;
double HumanPhantomPrimaryGeneratorAction::filtration=0.0;
double HumanPhantomPrimaryGeneratorAction::anodAngle=0.0;
double HumanPhantomPrimaryGeneratorAction::ripple=0.0;
std::string HumanPhantomPrimaryGeneratorAction::anodMaterial="";
double HumanPhantomPrimaryGeneratorAction::DAP_uGymm2=0.0;
G4double HumanPhantomPrimaryGeneratorAction::FCA=0.0;
double HumanPhantomPrimaryGeneratorAction::CTDI_uGy=0.0;
int HumanPhantomPrimaryGeneratorAction::rotations=0;
int HumanPhantomPrimaryGeneratorAction::helical_rotations=0;
double HumanPhantomPrimaryGeneratorAction::MAS=0.0;
G4bool HumanPhantomPrimaryGeneratorAction::useMASforDAP=false;
double HumanPhantomPrimaryGeneratorAction::pitch=1.0;

G4int HumanPhantomPrimaryGeneratorAction::IACTIVITY=0;//constant
G4int HumanPhantomPrimaryGeneratorAction::IKERMA=1;
G4int HumanPhantomPrimaryGeneratorAction::INUMBER=2;
G4int HumanPhantomPrimaryGeneratorAction::IEVENTS=3;
G4int HumanPhantomPrimaryGeneratorAction::IQUANTA=3;//default
G4bool HumanPhantomPrimaryGeneratorAction::isIsotropicSource=true;
double HumanPhantomPrimaryGeneratorAction::emissionArea_mm2=1000000.0;//~1m x 1m
double HumanPhantomPrimaryGeneratorAction::activity_Bq=-1.0;
double HumanPhantomPrimaryGeneratorAction::radiationYield=1.0;//100%
double HumanPhantomPrimaryGeneratorAction::exposureTime_s=1.0;
double HumanPhantomPrimaryGeneratorAction::kerma_Gy=-1.0;
double HumanPhantomPrimaryGeneratorAction::massEnergyTransferCoefficient_cm2Perg=0.0289;//photons 1 MeV in air
double HumanPhantomPrimaryGeneratorAction::particleCounter=-1.0;

HumanPhantomPrimaryGeneratorAction::HumanPhantomPrimaryGeneratorAction()
  :beamKind("beamAlongZ"),worldLength(200.)
{
  G4int n_particle = 1;

  messenger= new HumanPhantomPrimaryGeneratorMessenger(this);

  particleGun = new G4ParticleGun(n_particle);
  
  phantom = (HumanPhantomConstruction*)
             G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  sliceIndex =0;
  angleIndex =0;
  eventIndex =0;
  rotations=0;
  helical_rotations=0;

  probability.push_back(0.1667);
  probability.push_back(0.1667);
  probability.push_back(0.1667);
  probability.push_back(0.1667);
  probability.push_back(0.1666);
  probability.push_back(0.1666);

  // Default Particle
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  
  particleGun->SetParticleDefinition(particleTable->FindParticle("gamma"));
  particleGun->SetParticleEnergy(0.080*MeV);
  
  SetBeam("beamAlongX");
  SetFieldCenterY(0.0*cm);SetFieldHeight(10.0*cm);SetFieldWidth(15.0*cm);
  SetSliceThickness(1.0*cm);//1.0*mm or even 0.6*mm
  SetAngleIncrement(1.0*deg);
  SetFanBeam("off");  

  isSpectrum=true;//to be handled via run.mac

  FCA=100.0*cm;
  worldLength=2.0*FCA;//real world length
}

HumanPhantomPrimaryGeneratorAction::~HumanPhantomPrimaryGeneratorAction()
{
  delete particleGun;
  delete messenger;
}

G4int HumanPhantomPrimaryGeneratorAction::GetNumberOfEventsPerScan(){
	return nbOfEventsPerScan;
}

G4double HumanPhantomPrimaryGeneratorAction::GetAzimuthalRotationAngle()
{
	return azimuthalRotationAngle;
}

G4double HumanPhantomPrimaryGeneratorAction::GetPolarRotationAngle()
{
	return polarRotationAngle;
}

G4double HumanPhantomPrimaryGeneratorAction::GetFieldCenterY()
{
	return fieldCenterY;
}

G4double HumanPhantomPrimaryGeneratorAction::GetFieldCenterX()
{
	return fieldCenterX;
}

G4double HumanPhantomPrimaryGeneratorAction::GetFieldHeight()
{
	return fieldHeight;
}

G4double HumanPhantomPrimaryGeneratorAction::GetFieldWidth()
{
	return fieldWidth;
}

G4double HumanPhantomPrimaryGeneratorAction::GetSliceThickness()
{
	return sliceThickness;
}

G4double HumanPhantomPrimaryGeneratorAction::GetAngleIncrement()
{
	return angleIncrement;
}

G4bool HumanPhantomPrimaryGeneratorAction::GetFanBeam()
{
	return isFanBeam;
}

G4bool HumanPhantomPrimaryGeneratorAction::GetHelicalScan()
{
	return helicalScan;
}

G4bool HumanPhantomPrimaryGeneratorAction::GetHalfField()
{
	return halfField;
}

G4bool HumanPhantomPrimaryGeneratorAction::GetDentalPanoramic()
{
	return dentalPanoramic;
}

G4bool HumanPhantomPrimaryGeneratorAction::IsSpectrum()
{
	return isSpectrum;
}

double HumanPhantomPrimaryGeneratorAction::GetDAP_uGymm2()
{
	return DAP_uGymm2;
}

G4double HumanPhantomPrimaryGeneratorAction::GetFCA()
{
	return FCA;
}

double HumanPhantomPrimaryGeneratorAction::GetCTDI()
{
	return CTDI_uGy;
}

double HumanPhantomPrimaryGeneratorAction::GetMAS()
{
	return MAS;
}

G4int HumanPhantomPrimaryGeneratorAction::GetIQUANTA(){
	return IQUANTA;
}

G4bool HumanPhantomPrimaryGeneratorAction::GetIsotropicSource()
{
	return isIsotropicSource;
}

double HumanPhantomPrimaryGeneratorAction::GetEmmissionArea_mm2()
{
	return emissionArea_mm2;
}

double HumanPhantomPrimaryGeneratorAction::GetActivity_Bq()
{
	return activity_Bq;
}

double HumanPhantomPrimaryGeneratorAction::GetExposureTime_s()
{
	return exposureTime_s;
}

double HumanPhantomPrimaryGeneratorAction::GetKerma_Gy()
{
	return kerma_Gy;
}

double HumanPhantomPrimaryGeneratorAction::GetMassEnergyTransferCoefficient_cm2Perg()
{
	return massEnergyTransferCoefficient_cm2Perg;
}

double HumanPhantomPrimaryGeneratorAction::GetParticleCounter()
{
	return particleCounter;
}

double HumanPhantomPrimaryGeneratorAction::GetRadiationYield()
{
	return radiationYield;
}

void HumanPhantomPrimaryGeneratorAction::GenerateMammoField(){	 

	//G4cout<<"Generate Mammo field "<<G4endl;

  G4double height = phantom->GetWorldSizeZ();  
  G4double diam = phantom->GetWorldSizeRadius();//here, its radius
  G4double distz=phantom->GetPointSourceOrFrontalBeamToDetectorDistance();

  double pi = 3.14159265359;
  //G4cout<<" PI= "<<pi<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<G4endl;//yeah, ok both!!!
  //init-------------
  //x0=0.0;
  //y0=0.0;
  //z0=0.0;
  G4double uin=0.0;
  G4double vin=0.0;
  G4double win=1.0;
  //-------------POLAR ANGLE-----------
  G4double costmax = distz/sqrt(distz * distz + diam * diam);
  //emmision probability for angles less than theta max:
  G4double dom = (1.0 - costmax) / 2.0;// >0 and <1/2, costmax <1 and >0, thetamax<90!
  //now a random angle less than thetamax
  G4double r=G4UniformRand();
  r = r * dom;//assure costet < costmax
  G4double costet = 1.0 - 2.0 * r;//>0, positive z axis; 2*r-1;//<0 i.e. negativ z axis!!
  //-----------------------------
  G4double sintet = sqrt(1.0 - costet * costet);
  G4double tgtet = sintet / costet;
  G4double teta = abs(atan(tgtet));// -pi/2,pi/2; teta is polar angle and should be >0!
  //----------AZIMUTHAL ANGLE----------	  
  r = G4UniformRand();
  G4double phi2 = 2.0 * pi * r;
  
  //uin=sin(teta) * cos(phi2);
  //vin=sin(teta) * sin(phi2);
  //THIS IS REQUIRED BECAUSE teta = abs(atan(tgtet)) IS NOT GOOD COMPUTED BY LINUX C COMPILER!!!!!!!!!!!!!!!!!
  //I THINK ATAN IS PROBLEMATIC OR THE FACT WE GET A G4double...WHO KNOWS!!??!!!
  uin=sintet * cos(phi2);
  vin=sintet * sin(phi2);
  win = costet;

  x0 = 0.0;
  y0 = 0.0;
  z0=-0.5 * height;

  G4ThreeVector direction(uin,vin,win);
  particleGun->SetParticleMomentumDirection(direction);
}

void HumanPhantomPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{		
	if(HumanPhantomConstruction::isMammo){
		GenerateMammoField();
	}
	else{
	worldLength=2.0*FCA;//real world length

  if (beamKind == "beamAlongZ") GenerateBeamAlongZ();
  if (beamKind == "beamAlongY") GenerateBeamAlongY();
  if (beamKind == "beamAlongX") GenerateBeamAlongX();
  if (beamKind == "isotropicFlux") GenerateIsotropicFlux();
  if (beamKind == "rectangleField") GenerateRectangleField();
  if (beamKind == "CTScan"){

	  //===required number of slices for current field (Head, etc.)
	  int nSlices=fieldHeight/(pitch*sliceThickness);
	  //==============
	eventIndex++;//+1 event is begin
	GenerateCTScan();
	if (eventIndex>=nbOfEventsPerScan){
		rotations++;
		helical_rotations++;
		angleIndex++;
		eventIndex=0;
	}

	if (dentalPanoramic){
		if (angleIndex*angleIncrement>=180.0*deg){
	   angleIndex=0;
	   //sliceIndex++;
	   G4cout << "Simulation has covered the whole examination field. The next event is set to begin from starting position!"<< G4endl;
	   
	}
	}else{
	if (angleIndex*angleIncrement>=360.0*deg){
	   angleIndex=0;
	   sliceIndex++;

	   //============test if out of field and if it is, then RESET to starting position	
		if (sliceIndex>nSlices){
			//reset, start over , from starting position			
			sliceIndex=0;
			helical_rotations=0;
			G4cout << "Simulation has covered the whole examination field. The next event is set to begin from starting position!"<< G4endl;	
		}		
		//===============
	}
	}

  }

	}//is not mamo
	//Ntot required=Round (fieldHeight/sliceThikness)*360*nperscan
  particleGun->SetParticlePosition(G4ThreeVector(x0, y0,z0));

  //============
  if (isSpectrum){
	  incidentEnergy = aliasSample();
	  particleGun->SetParticleEnergy(incidentEnergy);
  } else{
	  incidentEnergy = particleGun->GetParticleEnergy();
  }
  //=============

//G4cout << " CENTERY= "<<fieldCenterY<<"; W= "<<fieldWidth<<"; H= "<<fieldHeight<< G4endl;
  particleGun->GeneratePrimaryVertex(anEvent);

}

G4double HumanPhantomPrimaryGeneratorAction::aliasSample(){
	   //(int nsbin, double[] xs_array,// )
	   // double[] ws_array, int[] ibin_array) {
		// "===============================================================
		// "
		// " samples from an alias table which must have been prepared
		// " using prepare_alias_table
		// "
		// "===============================================================

		// ;Copyright NRC;
		// implicit none;

		// $INTEGER nsbin,ibin_array(nsbin);
		// $REAL xs_array(0:nsbin),ws_array(nsbin);
	    double alias_sample = 0.0;
		double v1 = 0.0;
		double v2 = 0.0;
		double aj = 0.0;
		int j = 0;
		int nsbin=ndata;
		v1 = G4UniformRand();//random01();
		v2 = G4UniformRand();//random01();
		aj = 1.0 + v1 * nsbin;
		//Double dbl = new Double(aj);
		//j = dbl.intValue();// minim1 maxim nsbin!!
		j=aj;
		if (j > nsbin)
			j = nsbin; // " this happens only if $RANDOMSET produces
						// " numbers in (0,1]--------> is not the DEFAULT
						// case!!!
		aj = aj - j;
		if (aj > ws_array[j - 1]) {
			j = ibin_array[j - 1];
		}
		alias_sample = (1.0 - v2) * en_array[j - 1] + v2 * en_array[j];// ok, xs=0 biased

		alias_sample=alias_sample*MeV;

		return alias_sample;
}

void HumanPhantomPrimaryGeneratorAction::GenerateCTScan()
{	
	//======
	//Here, field width play the role of fan width to encompass the patient
	//The field height play the role of distance from start scan to end scan position relative to
	//fieldCenterY....middle of the scan
	//i.e. startScan=fieldCenterY+fieldHeight/2.0; endScan=fieldCenterY-fieldHeight/2.0;
	//===
   G4double phantomWorldSize=phantom->GetWorldSize();//geometry specific world size!!	
   G4double startScan=fieldCenterY+fieldHeight/2.0;

   //G4double endScan=fieldCenterY-fieldHeight/2.0;

   G4double r=G4UniformRand();
   G4double x=(2.0*r-1.0)*fieldWidth*0.5;//-W/2;+W/2
   
   //HalfField or assymmetrical exposure geometry:
   if (halfField){
	   x=(r-1.0)*fieldWidth*0.5;//-W/2;0
   }
   //////////////////////////////

   r=G4UniformRand();
   G4double y=(2.0*r-1.0)*sliceThickness;//fieldHeight;
   //G4double z=-0.5*(worldLength)*cm;//worldLength(200.)=preset~FSD+phthickness/2.0
   G4double z=-0.5*(worldLength);//already in cm DO NOT EXCEED WORLDLENGTH OR ERRORS WILL OCCUR!!
   G4double D=sqrt(x*x+y*y+z*z);

   G4double u = 0;//x/D;
   G4double v = 0;//y/D;
   G4double w = 1;//-z/D;//against zaxis to hit

   if (isFanBeam){
	u = x/D;
    v = y/D;
    w = -z/D;//against zaxis to hit
   }
   //The above are given by choose point in rectangle according to theory.

   //initial coordinates
   G4double xin=0.0;
   G4double yin=0.0;
   G4double zin=z;
   //do not want to set yin=fieldCenterY here, first we make rotations then translation for sake of simplicity.

   //ADJUST INITIAL COORDINATES IN CASE OF worldlength=2*FCA exceeds phantomWorldSize
   G4double D0=D*(worldLength-phantomWorldSize)/worldLength;
   if(D0>0.0){
	   zin=-0.5*phantomWorldSize;//DO NOT EXCEED WORLDLENGTH OR ERRORS WILL OCCUR!!
	   xin=D0*u;
	   yin=D0*v;
   }
   //-------END Adjustments

   //rotation by an incremental azimuthal angle phi of field central axis versus phantom symmetry axis.
   G4double phi = angleIncrement*angleIndex;
   ///
   if (dentalPanoramic){
	   phi = 270.0*deg + angleIncrement*angleIndex;//from LLAT to the PA to the RLAT...image receptor moves in front of patient
   }
   ////
   G4double u1=0.0;
   G4double v1=0.0;
   G4double w1=0.0;
   G4double x1=0.0;
   G4double z1=0.0;
   //Rotation about y axis=>azimuthal
   u1=u*cos(phi)+w*sin(phi);
   w1=-u*sin(phi)+w*cos(phi);
   v1=v;
   
   //adjust entry point
   x1=xin*cos(phi)+zin*sin(phi);
   z1=-xin*sin(phi)+zin*cos(phi);
   
   G4double uin=u1;
   G4double vin=v1;
   G4double win=w1;
   x0=x1;
   y0=startScan-sliceThickness/2.0 -sliceIndex*pitch*sliceThickness;

   //if helical scan we have different y0...all other variables are the same.
   if (helicalScan){
	   y0=startScan-sliceThickness/2.0 -helical_rotations*pitch*sliceThickness*angleIncrement/(360.0*deg);

	   //number of events per single rotation is n=360 degrees /angleIncrement in degrees
	   //table feed or table movement per rotation is c= p x t; p=pitch factor; t = sliceThickness
	   //Number of rotations Nr = H/c=H/(p x t); H is startScan-endScan distance.

	   //Therefore minimum number of events to cover exposure area: Nmin = Nr x n = H/(p x t) x 360/angleIncrement
	   //So, number of events per scan (how many succesive exposure before rotation) is: Nevscan = Nchoosen/Nmin
	   //which is N x p x t x angleIncerement/(H x 360)

	   //in helical scan, the focus is moved along central axis at every rotation. If not helical scan, the focus
	   //rotates 360 degrees then it is moved by c =  p x t and start again.
	   //The movement between two steps is delta = c/360 or more generally delta = c x angleIncrement/(360 degrees)
	   //which is delta = p x t x angleIncrement/(360 degrees).
   }
   //[[[[[[[[[[[[[[[[[[[[[[[[[[

   if (dentalPanoramic){
	   y0=startScan-sliceThickness/2.0;//one half rotation and that's it; Here slice thickness is the field height
   }

   z0=z1;
   G4ThreeVector direction(uin,vin,win);
   particleGun->SetParticleMomentumDirection(direction);
}

void HumanPhantomPrimaryGeneratorAction::GenerateRectangleField()
{
   G4double phantomWorldSize=phantom->GetWorldSize();//geometry specific world size!!
   G4double r=G4UniformRand();
   G4double x=(2.0*r-1.0)*fieldWidth*0.5;//-w/2,+w/2
   r=G4UniformRand();
   G4double y=(2.0*r-1.0)*fieldHeight*0.5;//-h/2,+h/2
   //G4double z=-0.5*(worldLength)*cm;//worldLength(200.)=preset~FSD+phthickness/2.0
   G4double z=-0.5*(worldLength);//already in cm DO NOT EXCEED WORLDLENGTH OR ERRORS WILL OCCUR!!
   G4double D=sqrt(x*x+y*y+z*z);
   G4double u = x/D;
   G4double v = y/D;
   G4double w = -z/D;//against zaxis to hit
   //The above are given by choose point in rectangle according to theory.

   //initial coordinates
   G4double xin=0.0;
   G4double yin=0.0;//fieldCenterY;
   G4double zin=z;
   //do not want to set yin=fieldCenterY here, first we make rotations then translation for sake of simplicity.

   //ADJUST INITIAL COORDINATES IN CASE OF worldlength=2*FCA exceeds phantomWorldSize
   G4double D0=D*(worldLength-phantomWorldSize)/worldLength;
   if(D0>0.0){
	   zin=-0.5*phantomWorldSize;//DO NOT EXCEED WORLDLENGTH OR ERRORS WILL OCCUR!!
	   xin=D0*u;
	   yin=D0*v;
   }
   //-------END Adjustments

   //now assume a general rotation with a polar angle teta and an azimuthal angle phi of field central axis versus phantom symmetry axis.
   G4double teta =polarRotationAngle;
   G4double phi = azimuthalRotationAngle;

   G4double u1=0.0;
   G4double v1=0.0;
   G4double w1=0.0;
   G4double x1=0.0;
   G4double y1=0.0;
   G4double z1=0.0;

   G4double u2=0.0;
   G4double v2=0.0;
   G4double w2=0.0;
   G4double x2=0.0;
   G4double y2=0.0;
   G4double z2=0.0;

   //BETTER: FIRST ROTATION about xAxis then ROTATION about yAxis..both axis are for main Coordonate System!!
   //rotation about xaxis=>polar   
   //what was x now is z and what was z now is y   
   w2=w*cos(teta)+v*sin(teta);
   v2=-w*sin(teta)+v*cos(teta);
   u2=u;

   z2=zin*cos(teta)+yin*sin(teta);
   y2=-zin*sin(teta)+yin*cos(teta);
   x2=xin;
   //Rotation about y axis=>azimuthal   
   u1=u2*cos(phi)+w2*sin(phi);
   w1=-u2*sin(phi)+w2*cos(phi);
   v1=v2;
   
   //adjust entry point
   x1=x2*cos(phi)+z2*sin(phi);
   z1=-x2*sin(phi)+z2*cos(phi);
   y1=y2;
   
   /* 
   //Rotation about y axis=>azimuthal
   u1=u*cos(phi)+w*sin(phi);
   w1=-u*sin(phi)+w*cos(phi);
   v1=v;
   
   //adjust entry point
   x1=xin*cos(phi)+zin*sin(phi);
   z1=-xin*sin(phi)+zin*cos(phi);
   y1=yin;

   //rotation about xaxis=>polar   
   //what was x now is z and what was z now is y   
   w2=w1*cos(teta)+v1*sin(teta);
   v2=-w1*sin(teta)+v1*cos(teta);
   u2=u1;

   z2=z1*cos(teta)+y1*sin(teta);
   y2=-z1*sin(teta)+y1*cos(teta);
   x2=x1;
   */
   //===================================================================

   /*
   G4double uin=u2;
   G4double vin=v2;
   G4double win=w2;
   x0=x2+fieldCenterX;
   y0=y2+fieldCenterY;//perform translation here!//y1;
   z0=z2;
   */
   G4double uin=u1;
   G4double vin=v1;
   G4double win=w1;
   x0=x1+fieldCenterX;
   y0=y1+fieldCenterY;//perform translation here!//y1;
   z0=z1;
  
   G4ThreeVector direction(uin,vin,win);
   particleGun->SetParticleMomentumDirection(direction);
}

void HumanPhantomPrimaryGeneratorAction::GenerateBeamAlongZ()
{
  G4double phantomWorldSize=phantom->GetWorldSize();//geometry specific world size!!
//  z0 = 0.5*(worldLength)*cm
  z0=0.5*(phantomWorldSize);//already in cm DO NOT EXCEED WORLDLENGTH OR ERRORS WILL OCCUR!!
  y0 = (phantomWorldSize)*(G4UniformRand()-0.5);//*cm;
  x0 = (phantomWorldSize)*(G4UniformRand()-0.5);//*cm;
 
   G4ThreeVector direction(0.,0.,-1.);
   
  particleGun->SetParticleMomentumDirection(direction);
}

void HumanPhantomPrimaryGeneratorAction::GenerateBeamAlongX()
{
	G4double phantomWorldSize=phantom->GetWorldSize();//geometry specific world size!!
  x0 = -0.5*(phantomWorldSize);//*cm;
  y0 = (phantomWorldSize)*(G4UniformRand()-0.5);//*cm;
  z0 = (phantomWorldSize)*(G4UniformRand()-0.5);//*cm;
 
  G4ThreeVector direction(1.,0.,0.);
  particleGun->SetParticleMomentumDirection(direction);
}
void HumanPhantomPrimaryGeneratorAction::GenerateBeamAlongY()
{
	G4double phantomWorldSize=phantom->GetWorldSize();//geometry specific world size!!
  y0 = 0.5*(phantomWorldSize);//*cm;
  x0 = (phantomWorldSize)*(G4UniformRand()-0.5);//*cm;
  z0 = (phantomWorldSize)*(G4UniformRand()-0.5);//*cm;
 
  G4ThreeVector direction(0.,-1.,0.);
  particleGun->SetParticleMomentumDirection(direction);
}

void HumanPhantomPrimaryGeneratorAction::GenerateIsotropicFlux()
{
	G4double phantomWorldSize=phantom->GetWorldSize();//geometry specific world size!!
  G4double random = G4UniformRand();
  G4double sum = 0.;
  G4int i = 0;

  while(sum<random){sum += probability[i]; i++;}
  
  if(i==1) 
    {
      z0 = -0.5*(phantomWorldSize-2.*cm);//*cm;
      y0 = (phantomWorldSize)*(G4UniformRand()-0.5);//*cm;
      x0 = (phantomWorldSize)*(G4UniformRand()-0.5);//*cm;
    }

  if(i==2) 
    {
      y0 = -0.5*(phantomWorldSize-2.*cm);//*cm;
      z0 = (phantomWorldSize)*(G4UniformRand()-0.5);//*cm;
      x0 = (phantomWorldSize)*(G4UniformRand()-0.5);//*cm;
    }

  if(i==3) 
    {
      x0 = -0.5*(phantomWorldSize-2.*cm);//*cm;
      z0 = (phantomWorldSize)*(G4UniformRand()-0.5);//*cm;
      y0 = (phantomWorldSize)*(G4UniformRand()-0.5);//*cm;
    }

  if (i==4)
    {
      z0 = 0.5*(phantomWorldSize-2.*cm);//*cm;
      y0 = (phantomWorldSize)*(G4UniformRand()-0.5);//*cm;
      x0 = (phantomWorldSize)*(G4UniformRand()-0.5);//*cm;
    }

 if(i==5) 
    {
      y0 = 0.5*(phantomWorldSize-2.*cm);//*cm;
      z0 = (phantomWorldSize)*(G4UniformRand()-0.5);//*cm;
      x0 = (phantomWorldSize)*(G4UniformRand()-0.5);//*cm;
    }

  if(i==6) 
    {
      x0 = 0.5*(phantomWorldSize-2.*cm);//*cm;
      z0 = (phantomWorldSize)*(G4UniformRand()-0.5);//*cm;
      y0 = (phantomWorldSize)*(G4UniformRand()-0.5);//*cm;
    } 

  G4double a,b,c;
  G4double n;
  do{
    a = (G4UniformRand()-0.5)/0.5;
    b = (G4UniformRand()-0.5)/0.5; 
    c = (G4UniformRand()-0.5)/0.5;
    n = a*a+b*b+c*c;
  }while(n > 1 || n == 0.0);
  n = std::sqrt(n);
  a /= n;
  b /= n;
  c /= n;

  G4ThreeVector direction(a,b,c);
  particleGun->SetParticleMomentumDirection(direction);  
}

void HumanPhantomPrimaryGeneratorAction::SetBeam(G4String beam)
{
  isCTScan=false;
  if((beam == "beamAlongZ")||(beam == "beamAlongX")||
  (beam == "beamAlongY")||(beam == "isotropicFlux")||(beam == "rectangleField")||
  (beam == "CTScan"))
		beamKind = beam;
  
  else G4cout<<"This option is not valid "<<
               "---> beamAlongZ/beamAlongY/beamAlongX/isotropicFlux/rectangleField"
             <<G4endl;

  if(beam == "CTScan")isCTScan=true;
}

