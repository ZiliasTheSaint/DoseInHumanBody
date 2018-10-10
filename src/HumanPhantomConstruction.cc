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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// D. Fulea, National Institute of Public Health Bucharest, Cluj-Napoca Regional Center, Romania
//
#include "G4UnitsTable.hh"
#include "globals.hh"
#include <map>

#include "HumanPhantomConstruction.hh"
#include "HumanPhantomRunAction.hh"
#include "HumanPhantomSD.hh"
#include "G4SDManager.hh"

#include "PhantomBuilder.hh"
#include "FemaleBuilder.hh"
#include "MaleBuilder.hh"
#include "PhantomHeadBuilder.hh"

#include "G4RunManager.hh"
#include "HumanPhantomMaterial.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PVPlacement.hh"
#include "G4Tubs.hh"

#include "MirdMapConstants.hh"

G4bool HumanPhantomConstruction::isMammo=false;

HumanPhantomConstruction::HumanPhantomConstruction()
{
	//static final kind variables:
	actualFemalePhantomMass = 71.2969*kg;//73.3201*kg;
	actualMalePhantomMass = 69.7241*kg;//71.7475*kg;
	actualPhantomHeight = 178.60 *cm;//174*cm;//cl+ct+ch0+ch1+ch2

	phantomHeight = 200*cm;//174*cm;//updated from messenger
	phantomMass = 100*kg;//73.3201*kg;//updated from messenger
	ageGroup="adult";//updated from messenger
 messenger = new HumanPhantomMessenger(this);
 material = new HumanPhantomMaterial();
 motherWorld=0;
 ftotal_diameter =14.0*cm;//cm
 ftotal_height=5.0*cm;//cm
 fUpperAir_height=40.0*cm;//cm
}

HumanPhantomConstruction::~HumanPhantomConstruction()
{
 delete material;
 delete messenger;
}

G4VPhysicalVolume* HumanPhantomConstruction::Construct()
{
  material -> DefineMaterials();
  G4VPhysicalVolume* result;

  if (model!="MAMMO"){
	  ///setting actual mass and height according to user choice!!!
	  if (ageGroup=="adult"){
		  actualFemalePhantomMass = 71.4614*kg;//73.3201*kg;
		  actualMalePhantomMass = 69.8885*kg;//71.7475*kg;
		  actualPhantomHeight = 178.60 *cm;
	  }else if (ageGroup=="age_15"){
		  actualFemalePhantomMass = 55.0503*kg;//73.3201*kg;
		  actualMalePhantomMass = 53.726*kg;//71.7475*kg;
		  actualPhantomHeight = 168.07 *cm;
	  }else if (ageGroup=="age_10"){
		  actualFemalePhantomMass = 31.9078*kg;//73.3201*kg;
		  actualMalePhantomMass = 31.3504*kg;//71.7475*kg;
		  actualPhantomHeight = 139.77 *cm;
	  }else if (ageGroup=="age_5"){
		  actualFemalePhantomMass = 19.1091*kg;//73.3201*kg;
		  actualMalePhantomMass = 18.7723*kg;//71.7475*kg;
		  actualPhantomHeight = 109.11 *cm;
	  }else if (ageGroup=="age_1"){
		  actualFemalePhantomMass = 9.28103*kg;//73.3201*kg;
		  actualMalePhantomMass = 9.11146*kg;//71.7475*kg;
		  actualPhantomHeight = 74.41 *cm;
	  }else if (ageGroup=="newborn"){
		  actualFemalePhantomMass = 3.50105*kg;//73.3201*kg;
		  actualMalePhantomMass = 3.43447*kg;//71.7475*kg;
		  actualPhantomHeight = 50.96 *cm;
	  }
	  ///////////SCALING all coordonates and dimension...PCXMC method
	  if (sex =="Female") 
		scaleXY=std::sqrt(actualPhantomHeight*phantomMass/(phantomHeight*actualFemalePhantomMass));
	  else
		scaleXY=std::sqrt(actualPhantomHeight*phantomMass/(phantomHeight*actualMalePhantomMass));
	  scaleZ=phantomHeight/actualPhantomHeight;
	//temp  for field assessement
	  //scaleZ=1.0;//////////////////////////////TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
	  //scaleXY=1.0;//TTTTTTTTTTTTTTTTTTTTTEEEEEEEEEEEEEEEEEEEEEEEEEEMMMMMMMMMMMMMMMMMMMMMMMMMPPPPPPPPPPPPPPPPP
	  //===========================
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String bodypartSD = "BodyPartSD";
  HumanPhantomSD* userPhantomSD = new HumanPhantomSD(bodypartSD);
  SDman->AddNewDetector(userPhantomSD);

  BasePhantomBuilder*  builder = 0;

  //////////////
  //SetPhantomModel("MIRDHead");//OUTSIDE CALL FROM run.mac
  //SetPhantomSex("Male");//OUTSIDE CALL FROM run.mac

  SetBodyPartSensitivity("Head",true);
  SetBodyPartSensitivity("Trunk",true);
  
  //Skeletal system
  SetBodyPartSensitivity("UpperSpine",true);
  SetBodyPartSensitivity("MiddleLowerSpine",true);
  SetBodyPartSensitivity("Skull",true);
  SetBodyPartSensitivity("FacialBones",true);//!!!!!!!!!!!
  SetBodyPartSensitivity("Pelvis",true);
  
  //Organs
  SetBodyPartSensitivity("Brain",true);
  SetBodyPartSensitivity("UrinaryBladder",true);
  SetBodyPartSensitivity("Spleen",true);
  SetBodyPartSensitivity("Pancreas",true);
  SetBodyPartSensitivity("Thyroid",true);//?????????????????????????????????????

  //Gastrointestinal System
  SetBodyPartSensitivity("Stomach",true);
  SetBodyPartSensitivity("UpperLargeIntestine",true);
  SetBodyPartSensitivity("LowerLargeIntestine",true);
    
  SetBodyPartSensitivity("LeftLeg",true);
  SetBodyPartSensitivity("RightLeg",true);
  SetBodyPartSensitivity("LeftLegBone",true);
  SetBodyPartSensitivity("RightLegBone",true);
  SetBodyPartSensitivity("LeftArmBone",true);
  SetBodyPartSensitivity("RightArmBone",true);
  SetBodyPartSensitivity("LeftLung",true);
  SetBodyPartSensitivity("RightLung",true);
  SetBodyPartSensitivity("LeftKidney",true);
  SetBodyPartSensitivity("RightKidney",true);
  SetBodyPartSensitivity("Heart",true);////??????????????NOT AVAILABLE YET!
  SetBodyPartSensitivity("Liver",true);//????????????????????????????????

  SetBodyPartSensitivity("RibCage",true);
  SetBodyPartSensitivity("LeftScapula",true);
  SetBodyPartSensitivity("RightScapula",true);
  SetBodyPartSensitivity("LeftAdrenal",true);
  SetBodyPartSensitivity("RightAdrenal",true);
  
  SetBodyPartSensitivity("LeftBreast",true);
  SetBodyPartSensitivity("RightBreast",true);
  //Genitalia
  SetBodyPartSensitivity("LeftOvary",true);
  SetBodyPartSensitivity("RightOvary",true);
  SetBodyPartSensitivity("Uterus",true);

  //36 organs and regions so far!!!

  SetBodyPartSensitivity("LeftTesticle",true);
  SetBodyPartSensitivity("RightTesticle",true);
  SetBodyPartSensitivity("LeftClavicle",true);
  SetBodyPartSensitivity("RightClavicle",true);
  SetBodyPartSensitivity("Esophagus",true);
  SetBodyPartSensitivity("GallBladder",true);
  SetBodyPartSensitivity("SmallIntestine",true);
  SetBodyPartSensitivity("Thymus",true);

  if (model == "MIRDHead")
    { 
      G4cout << "HeadBuilder instantiated" << G4endl;
      builder = new PhantomHeadBuilder;
      if (model ==  "MIRDHead") builder->SetModel("MIRD");
      //else if (model ==  "ORNLHead") builder->SetModel("ORNLMale");
    }
  else
  {  
      if (sex =="Female") 
	{ 
	  builder = new FemaleBuilder;
	  builder ->SetScaleXY(scaleXY);
	  builder ->SetScaleZ(scaleZ);
	  builder ->SetAgeGroup(ageGroup);
	  builder->SetModel(model);
	  G4cout <<model << " "<< sex << G4endl;
	}
      else if (sex == "Male") 
	{
	  builder = new MaleBuilder;
	  builder ->SetScaleXY(scaleXY);
	  builder ->SetScaleZ(scaleZ);
	  builder ->SetAgeGroup(ageGroup);
	  builder->SetModel(model);
	  G4cout <<model << " "<< sex << G4endl;
	}
  }
  //=====================================================
  motherWorld=ConstructWorld();
  builder->SetMotherVolume(motherWorld);//ConstructWorld());
  //======================================================
  builder->BuildHead("yellow", false, sensitivities["Head"]);  
  builder->BuildSkull("grey", false,sensitivities["Skull"]);
  builder->BuildBrain("red", true,sensitivities["Brain"]);
  builder->BuildFacial("red", true,sensitivities["FacialBones"]);

  if (model != "MIRDHead")
  { 
      builder->BuildTrunk("yellow", false, sensitivities["Trunk"]);  

      builder->BuildLeftLeg("yellow", false,sensitivities["LeftLeg"]);
      builder->BuildRightLeg("yellow", false,sensitivities["RightLeg"]);
      builder->BuildLeftArmBone("grey", true,sensitivities["LeftArmBone"]);
      builder->BuildRightArmBone("grey", true, sensitivities["RightArmBone"]);        
      builder->BuildLeftLegBone("grey", true,sensitivities["LeftLegBone"]);
      builder ->BuildRightLegBone("grey", true,sensitivities["RightLegBone"]);
      builder->BuildUpperSpine("grey", true,sensitivities["UpperSpine"]);   

      if (model == "MIRD") 
	  {
	   builder->BuildLeftScapula("grey", true, sensitivities["LeftScapula"]); 
	   builder->BuildRightScapula("grey", true, sensitivities["RightScapula"]);

	   builder->BuildLeftAdrenal("yellow", true, sensitivities["LeftAdrenal"]);
	   builder->BuildRightAdrenal("yellow", true, sensitivities["RightAdrenal"]);
	  }
  
      builder->BuildMiddleLowerSpine("grey", true,sensitivities["MiddleLowerSpine"]);  
     builder->BuildPelvis("grey", true,sensitivities["Pelvis"]); 	    
      builder->BuildStomach("orange", true,sensitivities["Stomach"]);
     builder->BuildUpperLargeIntestine("lightBlue", true,sensitivities["UpperLargeIntestine"]);
     builder->BuildLowerLargeIntestine("lightBlue", true,sensitivities["LowerLargeIntestine"]);
	builder->BuildRibCage("grey", true,sensitivities["RibCage"]);    
      builder->BuildSpleen("green", true,sensitivities["Spleen"]);
     builder->BuildPancreas("brown", true,sensitivities["Pancreas"]);       
	  
      builder->BuildLeftKidney("green", true,sensitivities["LeftKidney"]);
      builder->BuildRightKidney("green", true,sensitivities["RightKidney"]);

      builder->BuildUrinaryBladder("green", true,sensitivities["UrinaryBladder"]);       
      builder->BuildLeftLung("blue", true,sensitivities["LeftLung"]);
      builder->BuildRightLung("blue", true,sensitivities["RightLung"]);
	  
	  if(sex=="Female")
	  {
		builder->BuildLeftOvary("purple", true,sensitivities["LeftOvary"]);
		builder->BuildRightOvary("purple", true,sensitivities["RightOvary"]);
		builder->BuildUterus("purple", true,sensitivities["Uterus"]);
				
		builder->BuildLeftBreast("purple", true,sensitivities["LeftBreast"]); 
		builder->BuildRightBreast("purple", true,sensitivities["RightBreast"]);	    
      }
       
	  builder->BuildHeart("red", true,sensitivities["Heart"]);//	!NA yet! 
	  builder->BuildLeftClavicle("red", true,sensitivities["LeftClavicle"]);//	!NA yet!
	  builder->BuildRightClavicle("red", true,sensitivities["RightClavicle"]);//	!NA yet!
	  builder->BuildGallBladder("red", true,sensitivities["GallBladder"]);//	!NA yet!
	  builder->BuildSmallIntestine("red", true,sensitivities["SmallIntestine"]);//	!NA yet!
	  builder->BuildThymus("red", true,sensitivities["Thymus"]);//	!NA yet!
	  builder->BuildEsophagus("red", true,sensitivities["Esophagus"]);//	!NA yet!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
	  builder->BuildThyroid("orange", true,sensitivities["Thyroid"]); //	!NA yet!  
      builder->BuildLiver("orange", true,sensitivities["Liver"]); //	!NA yet! 
	  
	  if(sex=="Male"){
		  BuildMaleGenitalia();
	    builder->BuildLeftTesticle("red", true,sensitivities["LeftTesticle"]);//	!NA yet!
	    builder->BuildRightTesticle("red", true,sensitivities["RightTesticle"]);//	!NA yet! 
	  }
      
    }
  ////ICRP2007
  double wGonads=0.08;double wBones=0.13;//0.01 bone surface + 0.12 ActiveBoneMarrow
  double wIntestines=0.12;double wLungs=0.12;
  double wStomach=0.12;double wBreasts=0.12;
  double wBladder=0.04;double wBrain=0.01;
  double wLiver=0.04;double wThyroid=0.04;
  double wRemainder=0.18;//0.12+0.04 oesophagus+0.01skin+0.01salivaryGlands
  
  G4double mGonads=0.0;
  if(sex=="Female"){
	  mGonads = builder->getMassOfLeftOvary()+builder->getMassOfRightOvary();
  } else {
	  mGonads = builder->getMassOfLeftTesticle()+builder->getMassOfRightTesticle();
  }
  G4double mBones=builder->getMassOfSkull()+builder->getMassOfFacial()+builder->getMassOfLeftArmBone()+builder->getMassOfRightArmBone()+
	  builder->getMassOfLeftLegBone()+builder->getMassOfRightLegBone()+builder->getMassOfUpperSpine()+
	  builder->getMassOfLeftScapula()+builder->getMassOfRightScapula()+builder->getMassOfMiddleLowerSpine()+
	  builder->getMassOfPelvis()+builder->getMassOfRibCage()+builder->getMassOfLeftClavicle()+builder->getMassOfRightClavicle();
  G4double mIntestines=builder->getMassOfUpperLargeIntestine()+builder->getMassOfLowerLargeIntestine()+builder->getMassOfSmallIntestine();
  G4double mLungs=builder->getMassOfLeftLung()+builder->getMassOfRightLung();
  G4double mBreasts=builder->getMassOfLeftBreast()+builder->getMassOfRightBreast();
  G4double mRemainder=builder->getMassOfHead()+builder->getMassOfTrunk()+builder->getMassOfLeftLeg()+builder->getMassOfRightLeg()+
	  builder->getMassOfLeftAdrenal()+builder->getMassOfRightAdrenal()+builder->getMassOfSpleen()+builder->getMassOfPancreas()+
	  builder->getMassOfLeftKidney()+builder->getMassOfRightKidney()+builder->getMassOfHeart()+
	  builder->getMassOfEsophagus()+builder->getMassOfGallBladder()+builder->getMassOfThymus();
  if(sex=="Female"){
	  mRemainder=mRemainder+builder->getMassOfUterus();
  }
  //TOTAL=1!!
  //////////////AFTER BUILD ALL, we get mass in order to subtract!!!
  bodypartMass["logicalHead"]=builder->getMassOfHead();bodypartWeight["logicalHead"]=wRemainder*bodypartMass["logicalHead"]/mRemainder;//remainder
  bodypartMass["logicalSkull"]=builder->getMassOfSkull();bodypartWeight["logicalSkull"]=wBones*bodypartMass["logicalSkull"]/mBones;//skeletal;bone surface: ICRP2007:0.01; Bonemarrow =0.12=>BONE=0.13
  bodypartMass["logicalFacialBones"]=builder->getMassOfFacial();bodypartWeight["logicalFacialBones"]=wBones*bodypartMass["logicalFacialBones"]/mBones;
  bodypartMass["logicalBrain"]=builder->getMassOfBrain();bodypartWeight["logicalBrain"]=wBrain;//brain;ICRP2007:0.01

  bodypartMass["logicalEsophagus"]=builder->getMassOfEsophagus();bodypartWeight["logicalEsophagus"]=wRemainder*bodypartMass["logicalEsophagus"]/mRemainder;//remainder
  bodypartMass["logicalGallBladder"]=builder->getMassOfGallBladder();bodypartWeight["logicalGallBladder"]=wRemainder*bodypartMass["logicalGallBladder"]/mRemainder;//remainder
  bodypartMass["logicalThymus"]=builder->getMassOfThymus();bodypartWeight["logicalThymus"]=wRemainder*bodypartMass["logicalThymus"]/mRemainder;//remainder
  bodypartMass["logicalLeftClavicle"]=builder->getMassOfLeftClavicle();bodypartWeight["logicalLeftClavicle"]=wBones*bodypartMass["logicalLeftClavicle"]/mBones;//skeletal;ICRP2007:0.01+0.12
  bodypartMass["logicalRightClavicle"]=builder->getMassOfRightClavicle();bodypartWeight["logicalRightClavicle"]=wBones*bodypartMass["logicalRightClavicle"]/mBones;

  bodypartMass["logicalTrunk"]=builder->getMassOfTrunk();bodypartWeight["logicalTrunk"]=wRemainder*bodypartMass["logicalTrunk"]/mRemainder;//remainder;ICRP2007:0.12
  bodypartMass["logicalLeftLeg"]=builder->getMassOfLeftLeg();bodypartWeight["logicalLeftLeg"]=wRemainder*bodypartMass["logicalLeftLeg"]/mRemainder;//remainder;ICRP2007:0.12
  bodypartMass["logicalRightLeg"]=builder->getMassOfRightLeg();bodypartWeight["logicalRightLeg"]=wRemainder*bodypartMass["logicalRightLeg"]/mRemainder;;//remainder;ICRP2007:0.12
  bodypartMass["logicalLeftArmBone"]=builder->getMassOfLeftArmBone();bodypartWeight["logicalLeftArmBone"]=wBones*bodypartMass["logicalLeftArmBone"]/mBones;//skeletal;ICRP2007:0.01+0.12
  bodypartMass["logicalRightArmBone"]=builder->getMassOfRightArmBone();bodypartWeight["logicalRightArmBone"]=wBones*bodypartMass["logicalRightArmBone"]/mBones;//skeletal;ICRP2007:0.01+0.12
  bodypartMass["logicalLeftLegBone"]=builder->getMassOfLeftLegBone();bodypartWeight["logicalLeftLegBone"]=wBones*bodypartMass["logicalLeftLegBone"]/mBones;//skeletal;ICRP2007:0.01+0.12
  bodypartMass["logicalRightLegBone"]=builder->getMassOfRightLegBone();bodypartWeight["logicalRightLegBone"]=wBones*bodypartMass["logicalRightLegBone"]/mBones;//skeletal;ICRP2007:0.01+0.12
  bodypartMass["logicalUpperSpine"]=builder->getMassOfUpperSpine();bodypartWeight["logicalUpperSpine"]=wBones*bodypartMass["logicalUpperSpine"]/mBones;//skeletal;ICRP2007:0.01+0.12
  bodypartMass["logicalLeftScapula"]=builder->getMassOfLeftScapula();bodypartWeight["logicalLeftScapula"]=wBones*bodypartMass["logicalLeftScapula"]/mBones;//skeletal;ICRP2007:0.01+0.12
  bodypartMass["logicalRightScapula"]=builder->getMassOfRightScapula();bodypartWeight["logicalRightScapula"]=wBones*bodypartMass["logicalRightScapula"]/mBones;//skeletal;ICRP2007:0.01+0.12
  bodypartMass["logicalLeftAdrenal"]=builder->getMassOfLeftAdrenal();bodypartWeight["logicalLeftAdrenal"]=wRemainder*bodypartMass["logicalLeftAdrenal"]/mRemainder;//remainder:ICRP2007:0.12
  bodypartMass["logicalRightAdrenal"]=builder->getMassOfRightAdrenal();bodypartWeight["logicalRightAdrenal"]=wRemainder*bodypartMass["logicalRightAdrenal"]/mRemainder;//remainder:ICRP2007:0.12
  bodypartMass["logicalMiddleLowerSpine"]=builder->getMassOfMiddleLowerSpine();bodypartWeight["logicalMiddleLowerSpine"]=wBones*bodypartMass["logicalMiddleLowerSpine"]/mBones;//skeletal;ICRP2007:0.01+0.12
  bodypartMass["logicalPelvis"]=builder->getMassOfPelvis();bodypartWeight["logicalPelvis"]=wBones*bodypartMass["logicalPelvis"]/mBones;//skeletal;ICRP2007:0.01+0.12
  bodypartMass["logicalStomach"]=builder->getMassOfStomach();bodypartWeight["logicalStomach"]=wStomach;//stomach:ICRP2007:0.12
  
  bodypartMass["logicalUpperLargeIntestine"]=builder->getMassOfUpperLargeIntestine();bodypartWeight["logicalUpperLargeIntestine"]=wIntestines*bodypartMass["logicalUpperLargeIntestine"]/mIntestines;//remainder:ICRP2007:0.12
  bodypartMass["logicalLowerLargeIntestine"]=builder->getMassOfLowerLargeIntestine();bodypartWeight["logicalLowerLargeIntestine"]=wIntestines*bodypartMass["logicalLowerLargeIntestine"]/mIntestines;//remainder:ICRP2007:0.12
  bodypartMass["logicalSmallIntestine"]=builder->getMassOfSmallIntestine();bodypartWeight["logicalSmallIntestine"]=wIntestines*bodypartMass["logicalSmallIntestine"]/mIntestines;//remainder:ICRP2007:0.12
  
  bodypartMass["logicalRibCage"]=builder->getMassOfRibCage();bodypartWeight["logicalRibCage"]=wBones*bodypartMass["logicalRibCage"]/mBones;//skeletal;ICRP2007:0.01+0.12
  bodypartMass["logicalSpleen"]=builder->getMassOfSpleen();bodypartWeight["logicalSpleen"]=wRemainder*bodypartMass["logicalSpleen"]/mRemainder;//remainder:ICRP2007:0.12
  bodypartMass["logicalPancreas"]=builder->getMassOfPancreas();bodypartWeight["logicalPancreas"]=wRemainder*bodypartMass["logicalPancreas"]/mRemainder;//remainder:ICRP2007:0.12
  
  bodypartMass["logicalLiver"]=builder->getMassOfLiver();bodypartWeight["logicalLiver"]=wLiver;//liver:ICRP2007:0.04//TBD
  
  bodypartMass["logicalLeftKidney"]=builder->getMassOfLeftKidney();bodypartWeight["logicalLeftKidney"]=wRemainder*bodypartMass["logicalLeftKidney"]/mRemainder;//remainder:ICRP2007:0.12
  bodypartMass["logicalRightKidney"]=builder->getMassOfRightKidney();bodypartWeight["logicalRightKidney"]=wRemainder*bodypartMass["logicalRightKidney"]/mRemainder;//remainder:ICRP2007:0.12
  bodypartMass["logicalUrinaryBladder"]=builder->getMassOfUrinaryBladder();bodypartWeight["logicalUrinaryBladder"]=wBladder;//bladder:ICRP2007:0.04
  
  bodypartMass["logicalHeart"]=builder->getMassOfHeart();bodypartWeight["logicalHeart"]=wRemainder*bodypartMass["logicalHeart"]/mRemainder;//remainder:ICRP2007:0.12//TBD
  
  bodypartMass["logicalLeftLung"]=builder->getMassOfLeftLung();bodypartWeight["logicalLeftLung"]=wLungs*bodypartMass["logicalLeftLung"]/mLungs;//lung:ICRP2007:0.12
  bodypartMass["logicalRightLung"]=builder->getMassOfRightLung();bodypartWeight["logicalRightLung"]=wLungs*bodypartMass["logicalRightLung"]/mLungs;//lung:ICRP2007:0.12
  
  bodypartMass["logicalThyroid"]=builder->getMassOfThyroid();bodypartWeight["logicalThyroid"]=wThyroid;//thyroid:ICRP2007:0.04//TBD
  
  if(sex=="Female")
  {
	  bodypartMass["logicalLeftOvary"]=builder->getMassOfLeftOvary();bodypartWeight["logicalLeftOvary"]=wGonads*bodypartMass["logicalLeftOvary"]/mGonads;//gonads:ICRP2007:0.08
	  bodypartMass["logicalRightOvary"]=builder->getMassOfRightOvary();bodypartWeight["logicalRightOvary"]=wGonads*bodypartMass["logicalRightOvary"]/mGonads;//gonads:ICRP2007:0.08
	  bodypartMass["logicalUterus"]=builder->getMassOfUterus();bodypartWeight["logicalUterus"]=wRemainder*bodypartMass["logicalUterus"]/mRemainder;//remainder:ICRP2007:0.12
	  bodypartMass["logicalLeftBreast"]=builder->getMassOfLeftBreast();bodypartWeight["logicalLeftBreast"]=wBreasts*bodypartMass["logicalLeftBreast"]/mBreasts;//breasts:ICRP2007:0.12
	  bodypartMass["logicalRightBreast"]=builder->getMassOfRightBreast();bodypartWeight["logicalRightBreast"]=wBreasts*bodypartMass["logicalRightBreast"]/mBreasts;//breasts:ICRP2007:0.12
  } else{
	  //male
	  if(mGonads>0.0){
		  //testes are defined!!!
	   bodypartMass["logicalLeftTesticle"]=builder->getMassOfLeftTesticle();bodypartWeight["logicalLeftTesticle"]=wGonads*bodypartMass["logicalLeftTesticle"]/mGonads;
	   bodypartMass["logicalRightTesticle"]=builder->getMassOfRightTesticle();bodypartWeight["logicalRightTesticle"]=wGonads*bodypartMass["logicalRightTesticle"]/mGonads;;
	  }
  }

  result=builder->GetPhantom(); 
  delete builder;

  //================getting mass==========
  G4double mass=0.0;  
  //=========================
  mass=0.0;
  std::map<std::string,G4double>::iterator i = bodypartMass.begin();
  std::map<std::string,G4double>::iterator end = bodypartMass.end();
  while(i!=end)
  {		
      G4String bodypart = i->first;
      G4double masspart = i->second;

	  i++;
	  mass=mass+masspart;
  }
  //====================
  SetPhantomTotalMass(mass);
  //================  
  } else {//MAMMOGRAPHY CASE
	  result=ConstructMAMMO();
  }

  return result;
}

void  HumanPhantomConstruction::SetBodyPartSensitivity(G4String bodyPartName, G4bool bodyPartSensitivity)
{
  sensitivities[bodyPartName] = bodyPartSensitivity;
  if(bodyPartSensitivity==true) 
    G4cout << " >>> " << bodyPartName << " added as sensitive volume." << G4endl;
}

G4VPhysicalVolume* HumanPhantomConstruction::ConstructWorld()
{
  G4Material* air = material -> GetMaterial("Air");

  // World Volume
  G4double worldSize = 2. *m ;//increase to 4m total size
  //save its size
  fWorldSize=2.0*worldSize;//Total size!!

  G4Box* world = new G4Box("world", worldSize, worldSize, worldSize);

  G4LogicalVolume* logicWorld = new G4LogicalVolume(world, 
						    air, 
						    "logicalWorld", 0, 0,0);

  G4VPhysicalVolume* motherVolume = new G4PVPlacement(0,G4ThreeVector(),
						      "physicalWorld",
						      logicWorld,
						      0,
						      false,
						      0);

  // Visualization Attributes
  G4VisAttributes* WorldVisAtt = new G4VisAttributes(G4Colour(0.94,0.5,0.5));
    
  WorldVisAtt->SetForceSolid(false);
  //logicWorld->SetVisAttributes(G4VisAttributes::Invisible);
 
  return motherVolume;
}

G4VPhysicalVolume* HumanPhantomConstruction::BuildMaleGenitalia()
{
  G4Material* soft = material -> GetMaterial("soft_tissue");

  MirdMapConstants* mmc = new MirdMapConstants();  
  mmc->initTrunk();
  std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  G4double at = scaleXY*trk[ageGroup].At;
  mmc->initLegs();
  std::map<std::string,legs_struct> lgs = mmc->GetLegsMap();  
  G4double cpl = scaleZ*lgs[ageGroup].Cpl;
  mmc->initTestes();
  std::map<std::string,testes_struct> tst = mmc->GetTestesMap();    
  G4double c = scaleZ*tst[ageGroup].c;
  mmc->initSkin();
  std::map<std::string,skin_struct> skn = mmc->GetSkinMap();    
  G4double S = scaleZ*skn[ageGroup].S;//z here=>goes in height
  delete mmc;
  G4double z1 = 2.0*c+S;//abs value = 4.80 
  G4double rmax = 0.5*at;//0.5*at*(1.0+z1/cpl);//10.48 but z1 is <0=>max is 0.5*at
  //here Y=Z!!
  G4double maleXSize = 2.0*rmax;//scaleXY*2.*10.48 *cm ;
  G4double maleYSize = z1;//scaleZ*2.*2.40 *cm ;
  G4double maleZSize = rmax;//scaleXY*2.*5.24 *cm ;
  

  G4Box* maleGenitalia = new G4Box("MaleGenitalia", maleXSize/2., maleYSize/2., maleZSize/2.);

  G4LogicalVolume* logicmaleGenitalia = new G4LogicalVolume(maleGenitalia, 
						    soft, 
						    "logicalmaleGenitalia");

  G4ThreeVector pos = //G4ThreeVector(0*cm, -scaleZ*2.4*cm, scaleXY*5.24*cm);
	  G4ThreeVector(0*cm, -0.5*z1, 0.5*rmax);
  G4VPhysicalVolume* maleGenitaliaVolume = new G4PVPlacement(0,
						  pos,
						  "maleGenitaliaVolume",
						  logicmaleGenitalia,
						  motherWorld,
						  false,
						  0, false);//true);virtual region just to define softtissue there!
  
  //====================MUST ADD IT TO SD FOR A NICE SCORING ALORITHM!!!!!!!!!!
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  logicmaleGenitalia->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  //====================
  
  // Visualization Attributes
  G4VisAttributes* maleGenitaliaVisAtt = new G4VisAttributes(G4Colour(0.8,0.2,0.3));//red    
  maleGenitaliaVisAtt->SetForceSolid(true);
  //logicmaleGenitalia->SetVisAttributes(maleGenitaliaVisAtt);//G4VisAttributes::Invisible);
  logicmaleGenitalia->SetVisAttributes(G4VisAttributes::Invisible);

  //CONCLUSION: xMIRD=x;yMIRD=-z;zMIRD=y;
 
  return maleGenitaliaVolume;
}

void HumanPhantomConstruction::SetPhantomSex(G4String newSex)
{
  sex=newSex;

  if (sex == "Male")
    {
      G4cout << ">> Male Phantom will be built." << G4endl;
    }
  if (sex == "Female")
    {
      G4cout << ">> Female Phantom will be built." << G4endl;
    }
  if ((sex != "Female") && (sex != "Male"))
    G4cout << sex << " can not be defined!" << G4endl;
}

void HumanPhantomConstruction::SetPhantomModel(G4String newModel)
{
  model = newModel;
  HumanPhantomConstruction::isMammo=false;
  if (model == "MIRD")
    {
      G4cout<<" >> Phantom " << model << " will be built."<<G4endl;
    }
  
  if (model == "MIRDHead")
    {
      G4cout<<" >> Phantom " << model << " will be built."<<G4endl;
    }
   
  if (model == "MAMMO")
    {
	  HumanPhantomConstruction::isMammo=true;
      G4cout<<" >> Phantom " << model << " will be built."<<G4endl;
    }
}

G4VPhysicalVolume* HumanPhantomConstruction::ConstructMAMMO(){
    //The world is a surrounding cylinder.
	bool fCheckOverlaps=true;	

	G4double factive_diameter=ftotal_diameter;
	G4double factive_height=ftotal_height;

  G4double world_diameter = 1.0*(ftotal_diameter);
  G4double world_Z  = 1.0*(ftotal_height+fUpperAir_height);
  
  G4Tubs* solidWorldD=0;
  G4LogicalVolume* logicWorldD=0;
  G4VPhysicalVolume* physicalWorldD=0;
  //world solid  
  solidWorldD=new G4Tubs("World_solid", 0.0, 0.5*world_diameter, 0.5*world_Z, 0.0, 2.0*pi);
  
  //world logic
  G4Material* air = material -> GetMaterial("Air");
  G4Material* breast = material -> GetMaterial("adipose_glandular");

  logicWorldD = new G4LogicalVolume(
						solidWorldD,          //its solid
                        air,           //its material
                       "World_logic");            //its name
  //world physical
  physicalWorldD = new G4PVPlacement(0,                     //no rotation
                  G4ThreeVector(),       //at (0,0,0) CENTER PLACED AT THESE COORDINATES!
                  logicWorldD,            //its logical volume
                      "World_phys",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      fCheckOverlaps);        //overlaps checking

  
  
  //---some initializations:
  G4double deltah=world_Z;
  G4double zoffset=0.5*deltah;//redundant
  G4ThreeVector pos = G4ThreeVector(0, 0, -zoffset);//redundant
  G4double deltah_save=0.0;
  //---end initialisation

  //------------------Zero level. Source to detector layer
  G4Tubs* solidUpperAir=0;G4LogicalVolume* logicUpperAir=0;G4VPhysicalVolume* physicalUpperAir=0;
  deltah=deltah-fUpperAir_height;
  zoffset=0.5*deltah;
  pos = G4ThreeVector(0, 0, -zoffset);//its position ...at top, toward negative z-axis.
  solidUpperAir=new G4Tubs("SourceTop_solid", 0, 0.5*ftotal_diameter, 0.5*fUpperAir_height, 0.0, 2.0*pi);
  logicUpperAir =new G4LogicalVolume(solidUpperAir, air, "SourceTop_logic");                       
  physicalUpperAir=new G4PVPlacement(0, pos, logicUpperAir, "SourceTop_phys", logicWorldD, false, 0, fCheckOverlaps);          

  deltah=deltah-fUpperAir_height;//next zofsset take into account the last half layer
  //////////////////////
  G4Tubs* solidDetector=0;
  logicDetector=0;
  physicalDetector=0;

  deltah=deltah-factive_height;
  zoffset=0.5*deltah;
  pos = G4ThreeVector(0, 0, -zoffset);//its position ...at top, toward negative z-axis.  
  solidDetector=new G4Tubs("Detector_solid", 0.0, 0.5*factive_diameter, 0.5*factive_height, 0.0, 2.0*pi);                
  logicDetector = new G4LogicalVolume(solidDetector, breast, "Detector_logic");  
  physicalDetector = new G4PVPlacement(0, pos, logicDetector, "Detector_phys", logicWorldD, false, 0, fCheckOverlaps); 
  
  deltah=deltah-factive_height;//next zofsset take into account the last half layer

  G4VisAttributes* VisAtt3= new G4VisAttributes(G4Colour(0.8,0.2,0.3));//red!
  VisAtt3->SetVisibility(true);
  VisAtt3->SetForceSolid(true);
  logicDetector->SetVisAttributes(VisAtt3);
  
  //always return the physical World  
  return physicalWorldD;
}