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
// D. Fulea, National Institute of Public Health Bucharest, Cluj-Napoca Regional Center, Romania
//

#include "HumanPhantomRunAction.hh"
#include "G4ios.hh"
#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "HumanPhantomConstruction.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4VVisManager.hh"
#include "G4Box.hh"
#include "G4VisAttributes.hh"
#include "HumanPhantomPrimaryGeneratorAction.hh"
#include "G4RotationMatrix.hh"
#include "MirdMapConstants.hh"
//#include "XRayBuilder.hh"
#include"XRAYBuilder.hh"
HumanPhantomRunAction::HumanPhantomRunAction()
{

// add new units for dose 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);  
  
 phantom = (HumanPhantomConstruction*)
             G4RunManager::GetRunManager()->GetUserDetectorConstruction(); 
 gun = (HumanPhantomPrimaryGeneratorAction*)
                G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
}

HumanPhantomRunAction::~HumanPhantomRunAction()
{

}

void HumanPhantomRunAction::fillPerEvent(G4double EDet)
{
  //accumulate statistic
  //
  sumEDet += EDet;  sum2EDet += EDet*EDet;  
}

void HumanPhantomRunAction::BeginOfRunAction(const G4Run* aRun)
{
	//to be removed..put this in JAVA...why? it is good here
	G4double scaleXY=phantom->GetScaleXY();
    G4double scaleZ=phantom->GetScaleZ();
	HumanPhantomPrimaryGeneratorAction::fieldCenterX=scaleXY*HumanPhantomPrimaryGeneratorAction::fieldCenterX;
	HumanPhantomPrimaryGeneratorAction::fieldCenterY=scaleZ*HumanPhantomPrimaryGeneratorAction::fieldCenterY;
	HumanPhantomPrimaryGeneratorAction::fieldWidth=scaleXY*HumanPhantomPrimaryGeneratorAction::fieldWidth;
	HumanPhantomPrimaryGeneratorAction::fieldHeight=scaleZ*HumanPhantomPrimaryGeneratorAction::fieldHeight;

	if(!HumanPhantomConstruction::isMammo)
		G4cout << " Phantom weight/height scale factors: scaleXY= "<<scaleXY<<"; scaleZ= "<<scaleZ<< G4endl;
	//WOOOOOOOOOOORKS!!!!!!!!!
	//====================================
	bool dummy=false;//true;//will be false when XRayClass is ready

	//SPECTRUM
	if(HumanPhantomPrimaryGeneratorAction::IsSpectrum()){
      if(dummy){
		//TESTED SPECTRUM...READ FILE!

	    int number_of_lines = 0;
        std::string line;
        std::ifstream myfile("Data/givenSpectra/xray80kv17ua25al.spectrum");

	    if (!myfile.is_open()){
		  //G4cout<<"Error opening file! \n";
		  return;
	    }

		ndata=0;
		int dataCounter=0;//in line data counter
	
		double temp;
		std::vector<double> enV;
		std::vector<double> pV;
		while (myfile.good()){
			number_of_lines++;
			std::getline(myfile, line);//read line and put it in line variable.
        
			char* str=new char[line.length()+1];
			strcpy(str,line.c_str());//copy line into str.
			//---------------------------
			if (number_of_lines==1){
				G4cout<<"\n DUMMY: "<<line<<"\n";//print some info to know what file we just opened!
			}
			//-----------------------------
			if (number_of_lines==2){

				char* pch;
				pch=strtok(str," ,");//space and comma
				while (pch != NULL)
				{
					dataCounter++;
					if(dataCounter==1){
						std::istringstream buffer(pch);
						buffer>>ndata;	
					}

					if(dataCounter==2){
						std::istringstream buffer(pch);
						buffer>>temp;
						enV.push_back(temp);
					}

					if(dataCounter==3){
						std::istringstream buffer(pch);
						buffer>>temp;
						//pV.push_back(temp);//no probab data here!!

						dataCounter=0;//reset data counter
					}

					//G4cout<<"\n"<<pch;
					pch = strtok (NULL, " ,");//next data and also break the loop if next data = null
				}			
				//G4cout <<"\n"<<ndata<<"\n";//<<"  "<<enV.back()<<"  "<<pV.back()<<"\n";//

			} else if (number_of_lines>2){
				char* pch;
				pch=strtok(str," ,");//space and comma

				while (pch != NULL)
				{
					dataCounter++;
					if(dataCounter==1){
						std::istringstream buffer(pch);
						buffer>>temp;
						enV.push_back(temp);
					}
					if(dataCounter==2){
						std::istringstream buffer(pch);
						buffer>>temp;
						pV.push_back(temp);

						dataCounter=0;//reset data counter
					}

					pch = strtok (NULL, " ,");//next data and also break the loop if next data = null
				}
				//G4cout <<enV.back()<<"  "<<pV.back()<<"\n";
			}
		
		}//while (myfile.good())
		//G4cout << "\nNdata= " << ndata<<"\n";
		//teste
		//double a= 5.2;
		//int ia=a;
		//double b=5.8;
		//int ib=b;
		//G4cout << "\a: "<<a <<" ia: "<<ia<< " b: "<<b<<" ib: "<<ib <<"\n";
		//G4cout << "\nNumber of lines in text file: " << number_of_lines<<"\n";

		enV_array=new double[ndata];
		for (int i=0;i<ndata;i++)
			enV_array[i]=enV[i];
		prepareAliasSampling(ndata,pV);

		myfile.close();
	  } else {
		  xraybuilder = new XRayBuilder(
			  HumanPhantomPrimaryGeneratorAction::kv,//80.0
			  HumanPhantomPrimaryGeneratorAction::filtration,//2.5
			  HumanPhantomPrimaryGeneratorAction::anodAngle,//17,"W"
			  HumanPhantomPrimaryGeneratorAction::ripple,//0
			  HumanPhantomPrimaryGeneratorAction::anodMaterial//"W"
			  );
		  enV_array=xraybuilder->GetEnergyArray();
		  ndata = xraybuilder->GetBins();
		  double* prob=xraybuilder->GetProbArray();
		  photonsPerDAP=xraybuilder->Get_photons_per_mm2_per_uGy();
		  kermaPer_mas_at75cm=xraybuilder->GetAirKerma_per_mAs_at75cm();
		  prepareAliasSampling(ndata,prob);
	  }
	}//END SPECTRUM

	if(!HumanPhantomConstruction::isMammo){
	//==============
	G4double xcenter=0.0*cm;
	G4double ycenter=50*cm;
	G4double xbox=20*cm;
	G4double ybox=30*cm;
		
	//rotation--------=>AP,PA,LLAT,RLAT.-------------show FIELD setup
	G4double arot = HumanPhantomPrimaryGeneratorAction::GetAzimuthalRotationAngle();
	G4double prot = HumanPhantomPrimaryGeneratorAction::GetPolarRotationAngle();
	//---------------
	xcenter=HumanPhantomPrimaryGeneratorAction::GetFieldCenterX();//xcenter=xcenter*(phantom->GetScaleXY());//YESS
	ycenter=HumanPhantomPrimaryGeneratorAction::GetFieldCenterY();//ycenter=ycenter*(phantom->GetScaleZ());
	xbox=HumanPhantomPrimaryGeneratorAction::GetFieldWidth()*0.5;//xbox=xbox*(phantom->GetScaleXY());
	ybox=HumanPhantomPrimaryGeneratorAction::GetFieldHeight()*0.5;//ybox=ybox*(phantom->GetScaleZ());//half size!
	G4double zbox=std::min(xbox,ybox);
	zbox=zbox/100.0;//flat<=>make an rectangle
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if (pVVisManager) { 
      G4Box box("xrayRectangle", xbox, ybox, zbox);
	  G4VisAttributes va1(G4Colour::Red());
      va1.SetForceSolid(false);
	  pVVisManager->Draw(box,va1,G4Translate3D(xcenter,ycenter,0.0*cm)*G4RotateY3D(arot)*G4RotateX3D(prot));
  }

  G4int nbEventInRun = aRun->GetNumberOfEventToBeProcessed();
  HumanPhantomPrimaryGeneratorAction::SetNumberOfEvents(nbEventInRun);
  //======================================================================
  if(!HumanPhantomPrimaryGeneratorAction::isNumberOfEventsValid()){
	  G4cout << "Incomplete simulation: Number of events must be increased! "<< G4endl
	  <<"Number of events per scan: "<<HumanPhantomPrimaryGeneratorAction::GetNumberOfEventsPerScan()<< G4endl
	  <<"Minimum number of events required (for 1 number of event/scan): "<<HumanPhantomPrimaryGeneratorAction::minimumNbOfEvents<< G4endl;
  } else{
	  if (HumanPhantomPrimaryGeneratorAction::isCTScan){
	    G4cout << G4endl<<"Number of events per scan: "<<HumanPhantomPrimaryGeneratorAction::GetNumberOfEventsPerScan()<< G4endl;
		G4cout <<"Slice thickness: "<<G4BestUnit(HumanPhantomPrimaryGeneratorAction::GetSliceThickness(),"Length")<< G4endl;
		G4cout <<"Angle increment: "<<G4BestUnit(HumanPhantomPrimaryGeneratorAction::GetAngleIncrement(),"Angle")<< G4endl;
		G4cout <<"Is fan beam? : "<<HumanPhantomPrimaryGeneratorAction::GetFanBeam()<< G4endl;
		G4cout <<"Is helical scan? : "<<HumanPhantomPrimaryGeneratorAction::GetHelicalScan()<< G4endl;
		G4cout <<"Pitch factor : "<<HumanPhantomPrimaryGeneratorAction::pitch<< G4endl;
		G4cout <<"Is assymmetric half-field scan? : "<<HumanPhantomPrimaryGeneratorAction::GetHalfField()<< G4endl;
		G4cout <<"Is dental panoramic half-rotation scan? : "<<HumanPhantomPrimaryGeneratorAction::GetDentalPanoramic()<< G4endl;
	  }
	  else
		G4cout << G4endl<<"Number of events: "<<nbEventInRun<< G4endl;
  }
  double fca=HumanPhantomPrimaryGeneratorAction::GetFCA();
  double w=HumanPhantomPrimaryGeneratorAction::GetFieldWidth();
  double h=HumanPhantomPrimaryGeneratorAction::GetFieldHeight();
  G4cout << G4endl<<"Focus to phantom symmetry axis distance: "
	  <<G4BestUnit(HumanPhantomPrimaryGeneratorAction::GetFCA(),"Length") <<" [mm]: "<<fca<<G4endl;
  G4cout<<"Field width [mm]: "<<w<<" ;Field Height (or StartScan-EndScan Distance for CT) [mm]: "<<h<<G4endl;
  if (HumanPhantomPrimaryGeneratorAction::isCTScan){
	  G4cout <<"CTDI [uGy]: "
	  <<HumanPhantomPrimaryGeneratorAction::GetCTDI()<< G4endl;
  //else//nothing here: DAP will be RE-EVALUATED
	  double DLP=h*HumanPhantomPrimaryGeneratorAction::GetCTDI()/10.0/1000.0;
	  G4cout<<"DLP, Dose Length Product index (CTDI x length of scan) [mGy x cm]: "<<DLP<<G4endl;
  }
  ////////////////////////
  if(!HumanPhantomPrimaryGeneratorAction::IsSpectrum()){
	G4int mode = HumanPhantomPrimaryGeneratorAction::GetIQUANTA();

	if (mode==HumanPhantomPrimaryGeneratorAction::IACTIVITY){
		//activity based calculations
		double activity = HumanPhantomPrimaryGeneratorAction::GetActivity_Bq();//-1 default!!
		G4cout<<"Radiation source activity [Bq or quantas/second]: "<<activity<<G4endl;
		double time = HumanPhantomPrimaryGeneratorAction::GetExposureTime_s();
		G4cout<<"Exposure time [s]: "<<time<<G4endl;
		double yield = HumanPhantomPrimaryGeneratorAction::GetRadiationYield();
		G4cout<<"Radiation emmission probability (yield): "<<yield<<G4endl;
		G4bool isIsotropicSource =  HumanPhantomPrimaryGeneratorAction::GetIsotropicSource();
		G4cout<<"Is isotropic source?: "<<isIsotropicSource<<G4endl;
		if (!isIsotropicSource){
			double beamArea = HumanPhantomPrimaryGeneratorAction::GetEmmissionArea_mm2();
			G4cout<<"Beam area for non-isotropic source [mm2]: "<<beamArea<<G4endl;
		}
		
	} else if (mode==HumanPhantomPrimaryGeneratorAction::IKERMA){
		double kerma = HumanPhantomPrimaryGeneratorAction::GetKerma_Gy();//joule/kg
		G4cout<<"Kerma measured/computed at entry surface [Gy]: "<<kerma<<G4endl;
		//G4cout 
		// <<"Incident energy: "<<G4BestUnit(HumanPhantomPrimaryGeneratorAction::incidentEnergy,"Energy") << G4endl;
		//here incident energy is not updated!!!!!!
		
		double mentrCoeff = HumanPhantomPrimaryGeneratorAction::GetMassEnergyTransferCoefficient_cm2Perg();
		G4cout<<"Mass-energy transfer coefficient for photons at incident energy in surrounding environment material [cm2/g]: "<<mentrCoeff<<G4endl;
		
	} else if (mode==HumanPhantomPrimaryGeneratorAction::INUMBER){
		double particleCounter = HumanPhantomPrimaryGeneratorAction::GetParticleCounter();
		G4cout<<"Number of quantas hitting the target: "<<particleCounter<<G4endl;
	} else{
		G4cout<<"Number of quantas hitting the target equals the simulation number of runs!"<<G4endl;
	}
  }

  ///////////////////
  }else{
	  //is mamo
	  G4double FSD = phantom->GetPointSourceOrFrontalBeamToDetectorDistance();
	  G4double height = phantom->GetWorldSizeZ();  //total
      G4double diam = phantom->GetWorldSizeRadius();
	  diam=2.0*diam;
	  height=height-FSD;
	   G4cout << G4endl<<"Focus to breast distance: "
	  <<G4BestUnit(FSD,"Length")<< G4endl;
	   G4cout <<"Breast diameter: "
	  <<G4BestUnit(diam,"Length")<< G4endl;
	   G4cout <<"Breast height: "
	  <<G4BestUnit(diam,"Length")<< G4endl;

	   if(!HumanPhantomPrimaryGeneratorAction::IsSpectrum()){
			G4int mode = HumanPhantomPrimaryGeneratorAction::GetIQUANTA();

			if (mode==HumanPhantomPrimaryGeneratorAction::IACTIVITY){
				//activity based calculations
				double activity = HumanPhantomPrimaryGeneratorAction::GetActivity_Bq();//-1 default!!
				G4cout<<"Radiation source activity [Bq or quantas/second]: "<<activity<<G4endl;
				double time = HumanPhantomPrimaryGeneratorAction::GetExposureTime_s();
				G4cout<<"Exposure time [s]: "<<time<<G4endl;
				double yield = HumanPhantomPrimaryGeneratorAction::GetRadiationYield();
				G4cout<<"Radiation emmission probability (yield): "<<yield<<G4endl;
				G4bool isIsotropicSource =  HumanPhantomPrimaryGeneratorAction::GetIsotropicSource();
				G4cout<<"Is isotropic source?: "<<isIsotropicSource<<G4endl;
				if (!isIsotropicSource){
					double beamArea = HumanPhantomPrimaryGeneratorAction::GetEmmissionArea_mm2();
					G4cout<<"Beam area for non-isotropic source [mm2]: "<<beamArea<<G4endl;
				}
		
			} else if (mode==HumanPhantomPrimaryGeneratorAction::IKERMA){
				double kerma = HumanPhantomPrimaryGeneratorAction::GetKerma_Gy();//joule/kg
				G4cout<<"Kerma measured/computed at entry surface [Gy]: "<<kerma<<G4endl;
				//G4cout 
				// <<"Incident energy: "<<G4BestUnit(HumanPhantomPrimaryGeneratorAction::incidentEnergy,"Energy") << G4endl;
				double mentrCoeff = HumanPhantomPrimaryGeneratorAction::GetMassEnergyTransferCoefficient_cm2Perg();
				G4cout<<"Mass-energy transfer coefficient for photons at incident energy in surrounding environment material [cm2/g]: "<<mentrCoeff<<G4endl;
		
			} else if (mode==HumanPhantomPrimaryGeneratorAction::INUMBER){
				double particleCounter = HumanPhantomPrimaryGeneratorAction::GetParticleCounter();
				G4cout<<"Number of quantas hitting the target: "<<particleCounter<<G4endl;
			} else{
				G4cout<<"Number of quantas hitting the target equals the simulation number of runs!"<<G4endl;
			}
  }

  ///////////////////
	}
  G4int run_number = aRun->GetRunID();
  G4cout << "### Run " << run_number << " start." << G4endl;

 energyTotal["logicalHead"]=0.;
 energyTotal["logicalSkull"]=0.;
 energyTotal["logicalBrain"]=0.;

 energyTotal["logicalTrunk"]=0.;
 energyTotal["logicalLeftLeg"]=0.;
 energyTotal["logicalRightLeg"]=0.; 
 energyTotal["logicalLeftArmBone"]=0.;
 energyTotal["logicalRightArmBone"]=0.;
 energyTotal["logicalLeftLegBone"]=0.;
 energyTotal["logicalRightLegBone"]=0.;
 energyTotal["logicalUpperSpine"]=0.;
 energyTotal["logicalLeftScapula"]=0.; 
 energyTotal["logicalRightScapula"]=0.; 
 energyTotal["logicalLeftAdrenal"]=0.; 
 energyTotal["logicalRightAdrenal"]=0.; 
 energyTotal["logicalMiddleLowerSpine"]=0.;
 energyTotal["logicalPelvis"]=0.;
 energyTotal["logicalStomach"]=0.;
 energyTotal["logicalUpperLargeIntestine"]=0.;
 energyTotal["logicalLowerLargeIntestine"]=0.;
 energyTotal["logicalRibCage"]=0.;
 energyTotal["logicalSpleen"]=0.;
 energyTotal["logicalPancreas"]=0.;
 
 energyTotal["logicalLeftKidney"]=0.;
 energyTotal["logicalRightKidney"]=0.;
 energyTotal["logicalUrinaryBladder"]=0.;
 
 energyTotal["logicalHeart"]=0.0; //!!
 energyTotal["logicalLeftTesticle"]=0.0;//!!!!!!!!!!!!!
 energyTotal["logicalRightTesticle"]=0.0;//!!!!!!!!!!!!!
 energyTotal["logicalFacialBones"]=0.;//////////////
 energyTotal["logicalLeftClavicle"]=0.0;
 energyTotal["logicalRightClavicle"]=0.0;
 energyTotal["logicalGallBladder"]=0.0;
 energyTotal["logicalEsophagus"]=0.0;//NA
 energyTotal["logicalSmallIntestine"]=0.0;
 energyTotal["logicalThymus"]=0.0;

 energyTotal["logicalLeftLung"]=0.;
 energyTotal["logicalRightLung"]=0.;
 energyTotal["logicalThyroid"]=0.0;//!!!
 energyTotal["logicalLiver"]=0.0;//!!

 energyTotal["logicalLeftOvary"]=0.;
 energyTotal["logicalRightOvary"]=0.;
 energyTotal["logicalUterus"]=0.;
 energyTotal["logicalLeftBreast"]=0.;
 energyTotal["logicalRightBreast"]=0.; 
 
 energyTotal["logicalmaleGenitalia"]=0.0;//

 energy = energyTotal;//========
 energy2 = energyTotal;//=========
 dose = energyTotal;//=========
 bodypartMassMap=energyTotal;//=========
 totalPhantomEnergy=0.0;
 totalPhantomEnergy2=0.0;
 bodypartWeightMap=energyTotal;//=========0.0 Initialization to match the names!!!!

 sumEDet = sum2EDet = 0.;
}

//this will return the actual number of quantas of radiation field hitting the target as in real case scenario
//for example: number of photons hitting the target
double HumanPhantomRunAction::getScaleFactor(int numberOfEvents){
	double pi = 3.14159265359;
	G4double incidentEnergy = HumanPhantomPrimaryGeneratorAction::incidentEnergy;//not a spectrum

	//G4cout <<"Incident energy from scaleFactor: "<<G4BestUnit(HumanPhantomPrimaryGeneratorAction::incidentEnergy,"Energy") << G4endl;
	//yes, here it is the corect value!

	G4double FCA=1.0;
	G4double exposedArea=1.0;	

	if(HumanPhantomConstruction::isMammo){
		FCA=phantom->GetPointSourceOrFrontalBeamToDetectorDistance();
		double radius = phantom->GetWorldSizeRadius();//mm
		exposedArea = pi*radius*radius;//in mm2
	} else {
		FCA=HumanPhantomPrimaryGeneratorAction::GetFCA();//mm
		double w=HumanPhantomPrimaryGeneratorAction::GetFieldWidth();
		double h=HumanPhantomPrimaryGeneratorAction::GetFieldHeight();
		//no rotation-like CT here. It works with rectangle field!!!!
		exposedArea=w*h;//in mm2
	}
	
	G4int mode = HumanPhantomPrimaryGeneratorAction::GetIQUANTA();

	if (mode==HumanPhantomPrimaryGeneratorAction::IACTIVITY){
		//activity based calculations
		double activity = HumanPhantomPrimaryGeneratorAction::GetActivity_Bq();//-1 default!!
		double time = HumanPhantomPrimaryGeneratorAction::GetExposureTime_s();
		double yield = HumanPhantomPrimaryGeneratorAction::GetRadiationYield();
		G4bool isIsotropicSource =  HumanPhantomPrimaryGeneratorAction::GetIsotropicSource();
		double beamArea = 4.0*pi* FCA*FCA; //isotropic in mm2
		if (!isIsotropicSource){
			beamArea = HumanPhantomPrimaryGeneratorAction::GetEmmissionArea_mm2();
			//it is supposed to cover the targetArea (exposedArea)
		}
		//fluence = Axyxt/(4pir^2); Number of quantas = fluence x Area
		return activity*time*yield*exposedArea/beamArea;
	} else if (mode==HumanPhantomPrimaryGeneratorAction::IKERMA){
		double kerma = HumanPhantomPrimaryGeneratorAction::GetKerma_Gy();//joule/kg
		incidentEnergy = incidentEnergy/joule;//in Joule
		double mentrCoeff = HumanPhantomPrimaryGeneratorAction::GetMassEnergyTransferCoefficient_cm2Perg();
		mentrCoeff = mentrCoeff*10*10/0.001;//in mm2/kg
		//K = fluence (N/Area) x E x (mutr/rho)
		return kerma*exposedArea/(incidentEnergy*mentrCoeff);
	} else if (mode==HumanPhantomPrimaryGeneratorAction::INUMBER){
		return HumanPhantomPrimaryGeneratorAction::GetParticleCounter();
	}
	//if here just return the number of events as simulation equals real-case scenario
	return numberOfEvents;
}

void HumanPhantomRunAction::EndOfRunAction(const G4Run* aRun)
{
  G4cout << "Number of events = " << aRun->GetNumberOfEvent() << G4endl;

  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;

  if (HumanPhantomPrimaryGeneratorAction::isCTScan){
	G4cout << " Rotations: " << HumanPhantomPrimaryGeneratorAction::GetRotations() << G4endl;
  }

  if(!HumanPhantomConstruction::isMammo){
  totalRunEnergyDeposit(NbOfEvents);
  }
  else{
  G4double err = basicAnalyze(sumEDet, sum2EDet, NbOfEvents);
  sumEDet /= NbOfEvents; sum2EDet /= NbOfEvents;
  G4double rmsEAbs =err;

  G4double mass = phantom->GetDetector_logical()->GetMass();
  G4double dose = sumEDet/mass;
  G4double rmsDose = rmsEAbs;//%///mass;

  double scaleFactor=1.0;
  double pi = 3.14159265359;
  double DAP=HumanPhantomPrimaryGeneratorAction::GetDAP_uGymm2();
	    
	  //photonsPerDAP means photons_per_mm2_per_uGy; DAP is uGy*mm2=>we get PHOTONS
	  if(HumanPhantomPrimaryGeneratorAction::IsSpectrum())
	  {		
		if(HumanPhantomPrimaryGeneratorAction::useMASforDAP){
            double mAs=HumanPhantomPrimaryGeneratorAction::GetMAS();
			double fca=phantom->GetPointSourceOrFrontalBeamToDetectorDistance();
			double zk=750.0;//mm
			DAP=kermaPer_mas_at75cm*mAs*zk*zk/(fca*fca);//kerma at FCA
			
			double radius = phantom->GetWorldSizeRadius();//mm
			DAP=DAP*pi*radius*radius;
		}
		scaleFactor=photonsPerDAP*DAP;//energyperevent means energy/photon
	  }
	  else
		scaleFactor=getScaleFactor(NbOfEvents);//NbOfEvents;//total photons

	  //============================
	  //=======WR	  
	  G4ParticleGun* particleGun =gun->GetParticleGun();
	  G4String particleName = particleGun->GetParticleDefinition()->GetParticleName();
	  double wr = 10.0;//for generic neutrons and the rest!!
	  if (particleName=="gamma" || particleName=="e-" || particleName=="e+" || particleName=="mu-"|| particleName=="mu+") wr=1.0;
	  else if (particleName=="alpha"|| particleName=="GenericIon" || particleName=="He3"|| particleName=="deuteron"|| particleName=="triton") wr=20.0;//or heavy ions
	  else if (particleName=="proton"||particleName=="pi+"||particleName=="pi-") wr=2.0;//protons, charged pions
	  else if (particleName=="neutron") wr=10.0;//generic neutron
	  else wr=1.0;
	  //CANCER RISK MORTALITY ALL CANCERS=BEIR7 p311 MIDDLE 40 age ADULT!
	  double age[7];
	  age[0]=20.;
	  age[1]=30.;
      age[2]=40.;
      age[3]=50.;
      age[4]=60.;
      age[5]=70.;
      age[6]=80.;
      double cancerMortality_breast[7];
      cancerMortality_breast[0]=101.;
      cancerMortality_breast[1]=61.;
      cancerMortality_breast[2]=35.;
      cancerMortality_breast[3]=19.;
      cancerMortality_breast[4]=9.;
      cancerMortality_breast[5]=5.;
      cancerMortality_breast[6]=2.;
	  double ageValue = phantom->GetPhantomAge();
	  bool lowerThanValue = true;
	  int ageSize=7;
	  int lowerIndex = findNearestDataIndexFromArray(age, ageSize, ageValue, lowerThanValue);
	  if (lowerIndex==ageSize-1){
        //scaned array does not encompass value
        if (age[0]>ageValue){
            lowerIndex = 0;
        } else if (age[ageSize-1]<ageValue){
            lowerIndex = ageSize-2;
        }
     }     
     int higherIndex = lowerIndex + 1;
	 double breastCancer = linInt(age[lowerIndex],cancerMortality_breast[lowerIndex],
                               age[higherIndex],cancerMortality_breast[higherIndex],ageValue);

	  double mortalityRiskAge40=breastCancer;//35.0;//for FEMALE BREAST, per 100 mGy per 100000 cases//WARNING LOW LET ONLY!!!NOT GOOD FOR ALPHAs except we take into account wr!!
	  G4double mGy = 1.e-3*gray;
	  G4double mgy=dose*scaleFactor/mGy;
	  int cases=10.0*mgy*wr*mortalityRiskAge40/100.0;//*10 for million and /100 to 1 mGy
	  //=============================
     
	  string resultS="";	

	  if(HumanPhantomPrimaryGeneratorAction::IsSpectrum())
		G4cout << G4endl<< "DAP[uGy*mm2]: " <<DAP<< G4endl<<G4endl; 
     
	  ostringstream s;
	  s<<DAP;
	  G4String temp=s.str();
	  if(HumanPhantomPrimaryGeneratorAction::IsSpectrum())
		resultS=resultS+"DAP[uGy*mm2]" + " = " + temp   + "\n"+ "\n";
	  s.str("");
	  s.clear();

	  G4cout
     << G4endl<< "--------------------End of Run------------------------------"<< G4endl
     << G4endl<< " Energy in Breast : " << G4BestUnit(sumEDet*scaleFactor,"Energy")
     << " +- "                          <<rmsEAbs<<" %"
     << G4endl;

	  s<<rmsEAbs;
	  temp=G4BestUnit(sumEDet*scaleFactor,"Energy");
	  G4String temp2=s.str();
	  resultS=resultS+"Energy in Breast"+ " = " + temp  + " +- " +temp2+" %"  + "\n";
	  s.str("");
	  s.clear();

     G4cout
     << G4endl<< " Dose in Breast : " << G4BestUnit(dose*scaleFactor,"Dose")
     << " +- "                          <<rmsDose<<" %" 
	 << G4endl
	 
	 << G4endl
		  <<"----- Lifetime fatal cancer risk [cases/1 million population]: "<<cases<< " ;Age: "<<ageValue<<" -----"
	 <<G4endl

	 << G4endl<< "------------------------------------------------------------"<< G4endl
     << G4endl;
	 
	 temp=G4BestUnit(dose*scaleFactor,"Dose");
	 resultS=resultS+"Dose in Breast"+ " = " + temp  + " +- " +temp2+" %"  + "\n";

	 s<<cases;
	 temp=s.str();
	
	 ostringstream s111;s111<<ageValue;
	 G4String temp22=s111.str();
	 
	 //s<<ageValue;
	 //G4String temp22=s.str();
	 
	 resultS=resultS+"----- Lifetime fatal cancer risk [cases/1 million population]: "+ temp +" ;Age: "+ temp22 + " -----" + "\n";
	 s.str("");
	 s.clear();

	 //============print into file========================
	 std::ofstream outfile ("result.txt");
	 outfile <<resultS<<std::endl;
	 outfile.close();

	 G4cout << "Result file created in program folder! "<< G4endl<< G4endl;
  }//if is Mammo!!!
  
}


void HumanPhantomRunAction::Fill(G4String bodypartName, 
				      G4double energyDeposit)

{
 energyTotal[bodypartName] += energyDeposit;
 energy[bodypartName] += energyDeposit;
 energy2[bodypartName] += energyDeposit*energyDeposit;

 if (bodypartName=="logicalmaleGenitalia"){ 
	 energyTotal[bodypartName]=0.0;//
	 energy[bodypartName]=0.0;//
	 energy2[bodypartName]=0.0;//
	 bodypartMassMap[bodypartName] =1.0*kg;//virtual
 }else
 bodypartMassMap[bodypartName]=phantom->GetBodypartMass(bodypartName); 
 
 dose[bodypartName] += energyDeposit/bodypartMassMap[bodypartName];
 if (bodypartName=="logicalmaleGenitalia"){ 
	 dose[bodypartName]=0.0;//
 }
 if (bodypartName=="logicalmaleGenitalia"){ 
	 bodypartWeightMap[bodypartName] =0.0;//virtual
 }else
 bodypartWeightMap[bodypartName]=phantom->GetBodypartWeight(bodypartName);

 if (bodypartName!="logicalmaleGenitalia"){
  totalPhantomEnergy += energyDeposit;
  totalPhantomEnergy2 += energyDeposit*energyDeposit;
 }
}

//==============
bool HumanPhantomRunAction::isGoodEnoughToPrint(G4String bodyname){
	G4String sex = phantom->GetPhantomSex();
	if(sex=="Female"){
		if (bodyname=="LeftTesticle" || bodyname=="RightTesticle")
			return false;
	}
	if(sex=="Male"){
		if (bodyname=="LeftOvary" || bodyname=="RightOvary"
			||bodyname=="LeftBreast" || bodyname=="RightBreast"
			||bodyname=="Uterus"
			)
			return false;
	}
	return true;
}
//================

void HumanPhantomRunAction::totalRunEnergyDeposit(G4int n) 
{
  string resultS="";	
  std::map<std::string,G4double>::iterator i = energyTotal.begin();
  std::map<std::string,G4double>::iterator end = energyTotal.end();

  std::map<std::string,G4double>::iterator ie = energy.begin();
  std::map<std::string,G4double>::iterator ie2 = energy2.begin();
  std::map<std::string,G4double>::iterator id = dose.begin();
  std::map<std::string,G4double>::iterator im = bodypartMassMap.begin();
  std::map<std::string,double>::iterator iw = bodypartWeightMap.begin();

  G4double effectiveDose =0.;
  //G4double totalOrganDose =0.;//according to BEIR VII, the cancer coefficient is for low-LET radiation (X,gamma, electrons etc) and must be multiplied 
  //with dose in Gy (or equivalent dose to take into account high-LET radiation). Still, effecive dose IS A USEFUL indicator since:
  //lets say we have cases 2 and dose 4; case 3 and dose 5. Naive (2+3)(4+5) is meaningles since is not equal with correct approach:
  //(2+3)(4+5) != 2 x 4 + 3 x 5. So total dose as sum over all organ doses is meaningless!!!! Total dose as total energy over phantom mass is the
  //correct total dose. Its equivalent dose (weighted by radiation type) or effective dose (weighted further by organ sensitivities) are correct
  //quantities total cancer cases should be multipy!
  //G4double toerror = 0.;

  G4double totalEnergyDepositInPhantom =0.;
  G4int k=0;
  double scaleFactor=1.0;
  while(i!=end)
    {
		
      G4String bodypart = i->first;
      G4double energyDep = i->second;
      
	  G4double energyPerEvent = ie->second;
	  G4double energyPerEvent2 = ie2->second;
	  G4double dosePerEvent = id->second;
	  G4double massPerEvent = im->second;//weird calling...HA HA HA
	  double weight = iw->second;
//G4cout<<"WEIGHT for eff dose= "<<weight<<G4endl;
	  //=============
	  G4double err = basicAnalyze(energyPerEvent, energyPerEvent2, n);
      energyPerEvent /= n; energyPerEvent2 /= n;
	  G4double rmsEnergyPerEvent =err;
	  //===============
	  double DAP=HumanPhantomPrimaryGeneratorAction::GetDAP_uGymm2();
	  if (HumanPhantomPrimaryGeneratorAction::isCTScan){
		  double CTDI=HumanPhantomPrimaryGeneratorAction::GetCTDI();//uGy
		  double w=HumanPhantomPrimaryGeneratorAction::GetSliceThickness();//mm default units...for pencil beam
		  double h=w;//mm. always h=sliceThickness
		  if (HumanPhantomPrimaryGeneratorAction::isFanBeam){
			  w=HumanPhantomPrimaryGeneratorAction::GetFieldWidth();//mm default units!
		  }
		  DAP=h*w*CTDI;//per rotation, because CTDI means per rotation
		  DAP=DAP*HumanPhantomPrimaryGeneratorAction::GetRotations();
	  }
	  
	  //photonsPerDAP means photons_per_mm2_per_uGy; DAP is uGy*mm2=>we get PHOTONS
	  if(HumanPhantomPrimaryGeneratorAction::IsSpectrum())
	  {		
		if((HumanPhantomPrimaryGeneratorAction::useMASforDAP)&&
			(!HumanPhantomPrimaryGeneratorAction::isCTScan)){
            double mAs=HumanPhantomPrimaryGeneratorAction::GetMAS();
			double fca=HumanPhantomPrimaryGeneratorAction::GetFCA();//mm
			double zk=750.0;//mm
			DAP=kermaPer_mas_at75cm*mAs*zk*zk/(fca*fca);//kerma at FCA
			double w=HumanPhantomPrimaryGeneratorAction::GetFieldWidth();
			double h=HumanPhantomPrimaryGeneratorAction::GetFieldHeight();
			DAP=DAP*w*h;
		}
		scaleFactor=photonsPerDAP*DAP;//energyperevent means energy/photon
	  }
	  else
		scaleFactor=getScaleFactor(n);//n;//total photons
	  //=============================
	  
	  int length=bodypart.length();
	  G4String newBodypart=bodypart.substr(7,length);
	  
	  G4String temp="";
	  G4String temp2="";
	  if (isGoodEnoughToPrint(newBodypart) && newBodypart!="maleGenitalia")
	  {
	  G4cout << "Energy in " <<newBodypart << " = "  
	     << G4BestUnit(energyPerEvent*scaleFactor,"Energy") 
		  << " +- " <<rmsEnergyPerEvent<<" %"
	     << G4endl;

	  ostringstream s;
	  s<<rmsEnergyPerEvent;
	  temp=G4BestUnit(energyPerEvent*scaleFactor,"Energy");
	  temp2=s.str();
	  resultS=resultS+"Energy in "+newBodypart + " = " + temp  + " +- " +temp2+" %"  + "\n";

	  //=============
	  dosePerEvent=dosePerEvent/n;
	  //===============
  	  G4cout		  
		 // << "Dose in " <<bodypart << " = "  
		  << "Dose in " <<newBodypart << " = "  
	     << G4BestUnit(dosePerEvent*scaleFactor,"Dose") 
		<< " +- " <<rmsEnergyPerEvent<<" %"
		<<" Organ Mass= "<<G4BestUnit(massPerEvent,"Mass")
	     << G4endl;

	  temp=G4BestUnit(dosePerEvent*scaleFactor,"Dose");
	  temp2=s.str();
	  G4String temp3=G4BestUnit(massPerEvent,"Mass");
	  resultS=resultS+" Dose in "+newBodypart + " = " + temp  + " +- " +temp2+" %"
		  +" ;organ mass= "+temp3+ "\n";
	  //========================
	  //effectiveDose=effectiveDose+energyPerEvent*weight;//at the end, we will multiply with scale factor and divide with total mass!!!
	  effectiveDose=effectiveDose+dosePerEvent*weight;//at the end, we will multiply with scale factor ACCORDING TO DEFINITION OF EFFECTIVE DOSE
	  //G4cout << "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWW: " <<weight<< G4endl; //ok!!
	  //================
	  //================
	  //totalOrganDose = totalOrganDose+dosePerEvent;//at the end we will multiply with scaling factor
	  //toerror = toerror+dosePerEvent*dosePerEvent*rmsEnergyPerEvent*rmsEnergyPerEvent;//stdev2 x 100^2 = sum of d2 x p2; p = 100  x stdev/d, i.e
	  //toerror = sum of stdevs squared divided by 10000.
	  //============
	  }
      i++;
	  ie++;ie2++;id++;im++;iw++;//==========================
      k++;
      totalEnergyDepositInPhantom += energyDep;

	  if((HumanPhantomPrimaryGeneratorAction::IsSpectrum()) && (i==end))
	  {
		G4cout << G4endl<< "DAP[uGy*mm2]: " <<DAP<< G4endl<<G4endl; 
		
		ostringstream s;
		s<<DAP;
		temp=s.str();
	    
		resultS=resultS+"\n"+ "DAP[uGy*mm2]: " + temp   + "\n"+ "\n";

	  }
    }
	  //==================
	  G4String ageGroup = phantom->GetAgeGroup();
	  G4double scaleZ = phantom->GetScaleZ();
	  MirdMapConstants* mmc = new MirdMapConstants();
	  mmc->initTrunk();
	  std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();
	  G4double ct = scaleZ*trk[ageGroup].Ct;
	  mmc->initLegs();
      std::map<std::string,legs_struct> lgs = mmc->GetLegsMap(); 
      G4double cl=scaleZ*lgs[ageGroup].Cl;
	  mmc->initHead();
	  std::map<std::string,head_struct> hed = mmc->GetHeadMap();  
      G4double ch0 = scaleZ*hed[ageGroup].Ch0;
      G4double ch1 = scaleZ*hed[ageGroup].Ch1;
      G4double ch2 = scaleZ*hed[ageGroup].Ch2;

	  G4double geomHeight = ct+cl+ch0+ch1+ch2;
 	  //=============
	  G4double err = basicAnalyze(totalPhantomEnergy, totalPhantomEnergy2, n);
      totalPhantomEnergy /= n; totalPhantomEnergy2 /= n;
	  G4double rmsTotalEnergyPerEvent =err;
	  //===============
  	  G4cout << "Total Energy deposit in the body is: " 
	     << G4BestUnit(totalPhantomEnergy*scaleFactor,"Energy") 
		<< " +- " <<rmsTotalEnergyPerEvent<<" %"
	     << G4endl;
	  
	  ostringstream s;
	  s<<rmsTotalEnergyPerEvent;
	  G4String temp=G4BestUnit(totalPhantomEnergy*scaleFactor,"Energy");
	  G4String temp2=s.str();
	  resultS=resultS+"Total Energy deposit in the body is: " + temp  + " +- " +temp2+" %"  + "\n";

	  G4double tmass=phantom->GetPhantomTotalMass();
	  G4cout << "Total phantom mass (from geometry scaling routines): " 
	     << G4BestUnit(tmass,"Mass") 
	     << G4endl;

	  temp=G4BestUnit(tmass,"Mass");
	  resultS=resultS+"Total phantom mass (from geometry scaling routines): " + temp  + "\n";
	  //============	  
	  G4cout << "Phantom height (from geometry scaling routines) [cm]: " 
	     << geomHeight/cm 
	     << G4endl;
	  
	  ostringstream sS;
	  sS<<geomHeight/cm;
	  temp2=sS.str();
	  resultS=resultS+"Phantom height (from geometry scaling routines) [cm]: " + temp2  + "\n";
	  //=============
	  G4double tdose=totalPhantomEnergy/tmass;
	  G4cout << "Total dose deposited in the body (absorbed dose = total deposited energy/phantom mass): " 
	     << G4BestUnit(tdose*scaleFactor,"Dose") 
		<< " +- " <<rmsTotalEnergyPerEvent<<" %"
	     << G4endl;

	  //ostringstream s;
	  //s<<rmsTotalEnergyPerEvent;
	  temp=G4BestUnit(tdose*scaleFactor,"Dose");
	  temp2=s.str();
	  resultS=resultS+"Total dose deposited in the body (absorbed dose = total deposited energy/phantom mass): " + temp  + " +- " +temp2+" %"  + "\n";
	  //WR=========================
	  G4ParticleGun* particleGun =gun->GetParticleGun();
	  G4String particleName = particleGun->GetParticleDefinition()->GetParticleName();
	  double wr = 10.0;//for generic neutrons and the rest!!
	  if (particleName=="gamma" || particleName=="e-" || particleName=="e+" || particleName=="mu-"|| particleName=="mu+") wr=1.0;
	  else if (particleName=="alpha"|| particleName=="GenericIon" || particleName=="He3"|| particleName=="deuteron"|| particleName=="triton") wr=20.0;//or heavy ions
	  else if (particleName=="proton"||particleName=="pi+"||particleName=="pi-") wr=2.0;//protons, charged pions
	  else if (particleName=="neutron") wr=10.0;//generic neutron
	  else wr=1.0;
	  //CANCER RISK MORTALITY ALL CANCERS=BEIR7 p311 MIDDLE 40 age ADULT!
	  double age[11];//array declaration
    age[0]=0.;
	age[1]=5.;
	age[2]=10.;
	age[3]=15.;
	age[4]=20.;
    age[5]=30.;
    age[6]=40.;
    age[7]=50.;
    age[8]=60.;
    age[9]=70.;
    age[10]=80.;
    double cancerMortality_male[11];
    cancerMortality_male[0]=1099.;
	cancerMortality_male[1]=852.;
	cancerMortality_male[2]=712.;
	cancerMortality_male[3]=603.;
	cancerMortality_male[4]=511.;
    cancerMortality_male[5]=381.;
    cancerMortality_male[6]=377.;
    cancerMortality_male[7]=360.;
    cancerMortality_male[8]=319.;
    cancerMortality_male[9]=250.;
    cancerMortality_male[10]=153.;
    double cancerMortality_female[11];
    cancerMortality_female[0]=1770.;
	cancerMortality_female[1]=1347.;
	cancerMortality_female[2]=1104.;
	cancerMortality_female[3]=914.;
	cancerMortality_female[4]=762.;
    cancerMortality_female[5]=542.;
    cancerMortality_female[6]=507.;
    cancerMortality_female[7]=469.;
    cancerMortality_female[8]=409.;
    cancerMortality_female[9]=317.;
    cancerMortality_female[10]=190.;
	double ageValue = phantom->GetPhantomAge();
    bool lowerThanValue = true;
    int ageSize=7;
    int lowerIndex = findNearestDataIndexFromArray(age, ageSize, ageValue, lowerThanValue);
    if (lowerIndex==ageSize-1){
        //scaned array does not encompass value
        if (age[0]>ageValue){
            lowerIndex = 0;
        } else if (age[ageSize-1]<ageValue){
            lowerIndex = ageSize-2;
        }
    }
    //lowerThanValue = false;
    int higherIndex = lowerIndex + 1;//findNearestDataIndexFromArray(age, ageSize, ageValue, lowerThanValue);
    double maleCancer = linInt(age[lowerIndex],cancerMortality_male[lowerIndex],
                               age[higherIndex],cancerMortality_male[higherIndex],ageValue);
    double femaleCancer = linInt(age[lowerIndex],cancerMortality_female[lowerIndex],
                               age[higherIndex],cancerMortality_female[higherIndex],ageValue);
	  double mortalityRiskAge40=maleCancer;//377.0;//for MALE, per 100 mGy per 100000 cases
      G4String sex=phantom->GetPhantomSex();
	  if (sex =="Female") mortalityRiskAge40=femaleCancer;//507.;
	  G4double mGy = 1.e-3*gray;
	  G4double mgy=effectiveDose*scaleFactor/mGy;//totalOrganDose*scaleFactor/mGy;//effectiveDose*scaleFactor/mGy;//tdose*scaleFactor/mGy;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!whole body!!!
	  int cases=10.0*mgy*wr*mortalityRiskAge40/100.0;//*10 for million and /100 to 1 mGyWARNING: Low LET only, not good for alphas except multiplying with wr!!!!!
	  double pcases = 100*(10.0*mgy*wr*mortalityRiskAge40/100.0)/1000000;//%
	  G4cout //<< "mGy Dose= "<<mgy
		  << G4endl<<"----- Lifetime fatal cancer risk [cases/1 million population]: "<<cases<<" (or "<<pcases<<" %)"
		  <<" ;Age: "<<ageValue<<" ;Exposed to: "<<particleName<<" -----"<<G4endl;
	  //G4cout<<G4endl<<"LAR is computed from the sum of all equivalent organ doses (which for Low-LET radiation is simply the sum of all absorbed organ doses) using BEIR VII model."<<G4endl;BULSHIT!
	  ostringstream s1;
	  s1<<cases;
	  temp=s1.str();
	  ostringstream s111;s111<<ageValue;
	  G4String temp22=s111.str();
	  ostringstream ps1;
	  ps1<<pcases;
	  G4String ptemp = ps1.str();
	  resultS=resultS+"\n"+"----- Lifetime fatal cancer risk [cases/1 million population]: " + temp+" (or "+ptemp+" %)"
		  +" ;Age: "+temp22  + " -----"  + "\n";
	  //resultS=resultS+"LAR is computed from the sum of all equivalent organ doses (which for Low-LET radiation is simply the sum of all absorbed organ doses) using BEIR VII model."+"\n";
	  //-----------------
	 // G4double todose = totalOrganDose;
	 // if (toerror!=0.0)
	//	toerror = toerror/(todose*todose);
	//  if (toerror>=0.0)
	//	toerror = sqrt(toerror);
	//  ostringstream errors;
	//  errors<<toerror;
	//  G4cout << "Sum of all organ doses (absorbed dose) [Gy]: " 
	//     << G4BestUnit(todose*scaleFactor,"Dose") 
	//	<< " +- " <<toerror<<" %"//rmsTotalEnergyPerEvent<<" %"
	 //    << G4endl<< G4endl;
	//  temp=G4BestUnit(todose*scaleFactor,"Dose");
	//  temp2=errors.str();//s.str();
	 // resultS=resultS+"Sum of all organ doses (absorbed dose) [Gy]: " + temp  + " +- " +temp2+" %"  + "\n";
	  //------------------------------
	  G4double edose=effectiveDose;
	  /*G4cout << "Effective dose in body [Gy<->Sv]: " 
	     << G4BestUnit(edose*scaleFactor/tmass,"Dose") 
		<< " +- " <<rmsTotalEnergyPerEvent<<" %"
	     << G4endl<< G4endl;*/
	  G4cout << "Effective dose in body [Gy<->Sv]: " 
	     << G4BestUnit(edose*scaleFactor,"Dose") 
		<< " +- " <<rmsTotalEnergyPerEvent<<" %"
	     << G4endl<< G4endl;
	  temp=G4BestUnit(edose*scaleFactor,"Dose");
	  temp2=s.str();
	  resultS=resultS+"Effective dose in body [Gy<->Sv]: " + temp  + " +- " +temp2+" %"  + "\n";
	  //============print into file========================
	  std::ofstream outfile ("result.txt");
	  outfile <<resultS<<std::endl;
	  outfile.close();

	  //G4cout << "Result file created in program folder! "<< G4endl<< G4endl;
}

void HumanPhantomRunAction::prepareAliasSampling(int nsbin, std::vector<double> fs_array){
	// "================double[] ws_array, int[] ibin_array)=================
		// "
		// " inputs: nsbin: number of bins in the histogram
		// " fs_array: bin probabilities
		// "
		// " Note that we don't need the bin limits at this point, they
		// " are needed for the actual sampling (in alias_sample)
		// "
		// " outputs: ws_array, ibin_array: alias table ready for sampling
		// "
		// "====================================================================

		// ;Copyright NRC;

		// $INTEGER nsbin,ibin_array(nsbin);
		// $REAL fs_array(nsbin),ws_array(nsbin);
	    ws_array=new double[nsbin];
		ibin_array=new int[nsbin];

		int i = 0;
		int j_l = 0;
		int j_h = 0;
		double sum = 0.0;
		double aux = 0.0;

		bool exit1b = false;
		bool exit2b = false;

		sum = 0;
		for (i = 1; i <= nsbin; i++) {
			if (fs_array[i - 1] < 1.e-30) {
				fs_array[i - 1] = 1.e-30;
			}
			ws_array[i - 1] = -fs_array[i - 1];
			ibin_array[i - 1] = 1;
			sum = sum + fs_array[i - 1];
		}
		sum = sum / nsbin;

		// DO i=1,nsbin-1 [
		for (i = 1; i <= nsbin - 1; i++) {
			exit1b = false;
			exit2b = false;
			for (j_h = 1; j_h <= nsbin; j_h++) {
				if (ws_array[j_h - 1] < 0) {
					if (abs(ws_array[j_h - 1]) > sum) {
						// GOTO :AT_EXIT1:;
						exit1b = true;
						break;
					}
				}
			}
			if (!exit1b)
				j_h = nsbin;
			// :AT_EXIT1:

			for (j_l = 1; j_l <= nsbin; j_l++) {
				if (ws_array[j_l - 1] < 0) {
					if (abs(ws_array[j_l - 1]) < sum) {
						// GOTO :AT_EXIT2:;
						exit2b = true;
						break;
					}
				}
			}
			if (!exit2b)
				j_l = nsbin;
			// :AT_EXIT2:

			aux = sum - abs(ws_array[j_l - 1]);
			ws_array[j_h - 1] = ws_array[j_h - 1] + aux;
			ws_array[j_l - 1] = -ws_array[j_l - 1] / sum;
			ibin_array[j_l - 1] = j_h;

			if (i == nsbin - 1) {
				ws_array[j_h - 1] = 1;
			}

		}

		HumanPhantomPrimaryGeneratorAction::SetSpectrumAliasSamplingData(ws_array,ibin_array, enV_array, ndata);
}

void HumanPhantomRunAction::prepareAliasSampling(int nsbin, double* fs_array){
	// "================double[] ws_array, int[] ibin_array)=================
		// "
		// " inputs: nsbin: number of bins in the histogram
		// " fs_array: bin probabilities
		// "
		// " Note that we don't need the bin limits at this point, they
		// " are needed for the actual sampling (in alias_sample)
		// "
		// " outputs: ws_array, ibin_array: alias table ready for sampling
		// "
		// "====================================================================

		// ;Copyright NRC;

		// $INTEGER nsbin,ibin_array(nsbin);
		// $REAL fs_array(nsbin),ws_array(nsbin);
	    ws_array=new double[nsbin];
		ibin_array=new int[nsbin];

		int i = 0;
		int j_l = 0;
		int j_h = 0;
		double sum = 0.0;
		double aux = 0.0;

		bool exit1b = false;
		bool exit2b = false;

		sum = 0;
		for (i = 1; i <= nsbin; i++) {
			if (fs_array[i - 1] < 1.e-30) {
				fs_array[i - 1] = 1.e-30;
			}
			ws_array[i - 1] = -fs_array[i - 1];
			ibin_array[i - 1] = 1;
			sum = sum + fs_array[i - 1];
		}
		sum = sum / nsbin;

		// DO i=1,nsbin-1 [
		for (i = 1; i <= nsbin - 1; i++) {
			exit1b = false;
			exit2b = false;
			for (j_h = 1; j_h <= nsbin; j_h++) {
				if (ws_array[j_h - 1] < 0) {
					if (abs(ws_array[j_h - 1]) > sum) {
						// GOTO :AT_EXIT1:;
						exit1b = true;
						break;
					}
				}
			}
			if (!exit1b)
				j_h = nsbin;
			// :AT_EXIT1:

			for (j_l = 1; j_l <= nsbin; j_l++) {
				if (ws_array[j_l - 1] < 0) {
					if (abs(ws_array[j_l - 1]) < sum) {
						// GOTO :AT_EXIT2:;
						exit2b = true;
						break;
					}
				}
			}
			if (!exit2b)
				j_l = nsbin;
			// :AT_EXIT2:

			aux = sum - abs(ws_array[j_l - 1]);
			ws_array[j_h - 1] = ws_array[j_h - 1] + aux;
			ws_array[j_l - 1] = -ws_array[j_l - 1] / sum;
			ibin_array[j_l - 1] = j_h;

			if (i == nsbin - 1) {
				ws_array[j_h - 1] = 1;
			}

		}

		HumanPhantomPrimaryGeneratorAction::SetSpectrumAliasSamplingData(ws_array,ibin_array, enV_array, ndata);
}

double HumanPhantomRunAction::linInt(double x1, double y1, double x2, double y2,double x){
	double result = -1.0;
	double mn[2];
    // insucces
	mn[0] = -1.0;// m
	mn[1] = -1.0;// n
	double num = x1 - x2;
	if (num != 0.0) {
		mn[0] = (y1 - y2) / num;
		mn[1] = (x1 * y2 - y1 * x2) / num;
		result = mn[0] * x + mn[1];
	}
	return result;
}

int HumanPhantomRunAction::findNearestDataIndexFromArray(double* a, int arraySize, double value, bool lowerThanValue){
	bool b = true;
	int ip = 0;

	if (arraySize > 1) {
		while (b) {
			if (lowerThanValue) {
                if ((a[ip] <= value) && (a[ip + 1] > value)) {
					break;
				}
			} else {
				if (ip > 0)
					if ((a[ip] >= value) && (a[ip - 1] < value)) {
						break;
					}
			}

			ip++;
			if (ip == arraySize - 1) {
				b = false;
				break;
			}
        }
			//nearestposition = ip;// ----------------
        return ip;//a[ip];
    } else {
			//nearestposition = 0;// ---------------
        return 0;//a[0];
    }
}
