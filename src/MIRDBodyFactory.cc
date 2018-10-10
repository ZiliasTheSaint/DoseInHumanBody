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
//D. Fulea, National Institute of Public Health Bucharest, Cluj-Napoca Regional Center, Romania
//
#include "MIRDBodyFactory.hh"

#include "MIRDHead.hh"
#include "MIRDSkull.hh"
#include "MIRDBrain.hh"

#include "MIRDTrunk.hh"
#include "MIRDLeftLeg.hh"
#include "MIRDRightLeg.hh"
#include "MIRDLeftArmBone.hh"
#include "MIRDRightArmBone.hh"
#include "MIRDLeftLegBone.hh"
#include "MIRDRightLegBone.hh"
#include "MIRDUpperSpine.hh"
#include "MIRDLeftScapula.hh"
#include "MIRDRightScapula.hh"
#include "MIRDLeftAdrenal.hh"
#include "MIRDRightAdrenal.hh"
#include "MIRDMiddleLowerSpine.hh"
#include "MIRDPelvis.hh"
#include "MIRDStomach.hh"
#include "MIRDUpperLargeIntestine.hh"
#include "MIRDLowerLargeIntestine.hh"
#include "MIRDRibCage.hh"
#include "MIRDPancreas.hh"
#include "MIRDSpleen.hh"

#include "MIRDLiver.hh"//!

#include "MIRDLeftKidney.hh"
#include "MIRDRightKidney.hh"
#include "MIRDUrinaryBladder.hh"

#include "MIRDLeftLung.hh"
#include "MIRDRightLung.hh"
#include "MIRDHeart.hh"//!
#include "MIRDThyroid.hh"//!

#include "MIRDUterus.hh"
#include "MIRDLeftBreast.hh"
#include "MIRDRightBreast.hh"
#include "MIRDRightOvary.hh"
#include "MIRDLeftOvary.hh"

#include "MIRDLeftTesticle.hh"//NA
#include "MIRDRightTesticle.hh"
#include "MIRDLeftClavicle.hh"
#include "MIRDRightClavicle.hh"
#include "MIRDGallBladder.hh"
#include "MIRDEsophagus.hh"
#include "MIRDSmallIntestine.hh"
#include "MIRDThymus.hh"
#include "MIRDFacial.hh"

MIRDBodyFactory::MIRDBodyFactory()
{
  // Map with name of the organ and pointer to the MIRDOrgan class
  organ["Head"] = new MIRDHead();
  organ["Skull"] = new MIRDSkull();
  organ["Brain"] = new MIRDBrain(); 

  organ["Trunk"] = new MIRDTrunk(); 
  organ["LeftLeg"] = new MIRDLeftLeg();
  organ["RightLeg"] = new MIRDRightLeg();  
  organ["LeftArmBone"] = new MIRDLeftArmBone();
  organ["RightArmBone"] = new MIRDRightArmBone();  
  organ["LeftLegBone"] = new MIRDLeftLegBone();
  organ["RightLegBone"] = new MIRDRightLegBone();

  organ["UpperSpine"] = new MIRDUpperSpine();  
  organ["LeftScapula"]= new MIRDLeftScapula(); 
  organ["RightScapula"]= new MIRDRightScapula(); 
  organ["LeftAdrenal"]= new MIRDLeftAdrenal();  
  organ["RightAdrenal"]= new MIRDRightAdrenal(); 
  organ["MiddleLowerSpine"] = new MIRDMiddleLowerSpine();
  organ["Pelvis"]= new MIRDPelvis();
  organ["Stomach"] = new MIRDStomach();
  organ["UpperLargeIntestine"] = new MIRDUpperLargeIntestine();
  organ["LowerLargeIntestine"] = new MIRDLowerLargeIntestine();
  organ["RibCage"] = new MIRDRibCage(); 
  organ["Spleen"] = new MIRDSpleen(); 
  organ["Pancreas"] = new MIRDPancreas();
  
  organ["Liver"] = new MIRDLiver();

  organ["LeftKidney"] = new MIRDLeftKidney();
  organ["RightKidney"] = new MIRDRightKidney();
  organ["UrinaryBladder"] = new MIRDUrinaryBladder();

  organ["LeftLung"]= new MIRDLeftLung();
  organ["RightLung"] = new MIRDRightLung();
  
  organ["Heart"] = new MIRDHeart();
  organ["Thyroid"] = new MIRDThyroid();

  organ["Uterus"] = new MIRDUterus(); 
  organ["LeftOvary"] = new MIRDLeftOvary();
  organ["RightOvary"] = new MIRDRightOvary();  
  organ["LeftBreast"] = new MIRDLeftBreast();
  organ["RightBreast"] = new MIRDRightBreast();

  organ["LeftTesticle"] = new MIRDLeftTesticle();
  organ["RightTesticle"] = new MIRDRightTesticle();
  organ["LeftClavicle"] = new MIRDLeftClavicle();
  organ["RightClavicle"] = new MIRDRightClavicle();
  organ["GallBladder"] = new MIRDGallBladder();
  organ["Esophagus"] = new MIRDEsophagus();
  organ["SmallIntestine"] = new MIRDSmallIntestine();
  organ["Thymus"] = new MIRDThymus();
  organ["FacialBones"] = new MIRDFacial();
}

MIRDBodyFactory::~MIRDBodyFactory()
{
	delete organ["Head"]; organ["Head"]=0;
	delete organ["Brain"]; organ["Brain"]=0;
	delete organ["Skull"]; organ["Skull"] =0;
  delete organ["RightAdrenal"]; organ["RightAdrenal"]=0;
  delete organ["LeftAdrenal"]; organ["LeftAdrenal"]=0;
  delete organ["RightScapula"];organ["RightScapula"] =0;
  delete organ["LeftScapula"];organ["LeftScapula"] =0;
  delete organ["LeftBreast"]; organ["LeftBreast"]=0;
  delete organ["RightBreast"]; organ["RightBreast"]=0;
  delete organ["RightLegBone"]; organ["RightLegBone"]=0;
  delete organ["LeftLegBone"]; organ["LeftLegBone"]=0;
  delete organ["RightOvary"]; organ["RightOvary"]=0;
  delete organ["LeftOvary"]; organ["LeftOvary"]=0;
  delete organ["RightLung"]; organ["RightLung"] =0;
  delete organ["LeftLung"]; organ["LeftLung"]=0;
  delete organ["Uterus"]; organ["Uterus"]=0;
  delete organ["UrinaryBladder"]; organ["UrinaryBladder"]=0;
  delete organ["RightKidney"]; organ["RightKidney"] =0;  
  delete organ["LeftKidney"]; organ["LeftKidney"] =0; 
  delete organ["Pancreas"]; organ["Pancreas"] =0;  
  delete organ["Spleen"]; organ["Spleen"] =0; 
  delete organ["RibCage"]; organ["RibCage"] =0;  
  delete organ["LowerLargeIntestine"]; organ["LowerLargeIntestine"] =0; 
  delete organ["UpperLargeIntestine"]; organ["UpperLargeIntestine"] =0; 
  delete organ["Stomach"]; organ["Stomach"] =0;  
  delete organ["Pelvis"]; organ["Pelvis"] =0;  
  delete organ["MiddleLowerSpine"]; organ["MidlleLowerSpine"]=0;
  delete organ["UpperSpine"]; organ["UpperSpine"]=0;  
  delete organ["RightArmBone"]; organ["RightArmBone"] =0;
  delete organ["LeftArmBone"]; organ["LeftArmBone"] =0;  
  delete organ["Trunk"]; organ["Trunk"]=0;
  delete organ["RightLeg"]; organ["RightLeg"]=0;
  delete organ["LeftLeg"]; organ["LeftLeg"]=0;

  delete organ["Liver"];organ["Liver"]=0;
  delete organ["Heart"];organ["Heart"]=0;
  delete organ["Thyroid"];organ["Thyroid"]=0;

  delete organ["LeftTesticle"];organ["LeftTesticle"]=0;
  delete organ["RightTesticle"];organ["RightTesticle"]=0;
  delete organ["LeftClavicle"];organ["LeftClavicle"]=0;
  delete organ["RightClavicle"];organ["RightClavicle"]=0;
  delete organ["GallBladder"];organ["GallBladder"]=0;
  delete organ["Esophagus"];organ["Esophagus"]=0;
  delete organ["SmallIntestine"];organ["SmallIntestine"]=0;
  delete organ["Thymus"];organ["Thymus"]=0;
  delete organ["FacialBones"]; organ["FacialBones"]=0;
}



G4VPhysicalVolume* MIRDBodyFactory::CreateOrgan(const G4String& organ_name,G4VPhysicalVolume* motherVolume,
						  const G4String& colourName, G4bool visAttribute,
						  G4bool sensitivity)
{
	organ[organ_name]->SetScaleXY(scaleXY);
	organ[organ_name]->SetScaleZ(scaleZ);
	organ[organ_name]->SetAgeGroup(ageGroup);
 G4VPhysicalVolume* pV=organ[organ_name] -> Construct(organ_name,motherVolume,colourName, visAttribute, sensitivity);
 bodyVOL[organ_name]=organ[organ_name]->getCubicVolume();
 bodyRHO[organ_name]=organ[organ_name]->getDensity();
 
 return pV;
}

void MIRDBodyFactory::SetScaleXY(G4double val){
	scaleXY=val;
}

void MIRDBodyFactory::SetScaleZ(G4double val){
	scaleZ=val;
}

void MIRDBodyFactory::SetAgeGroup(G4String val){
	ageGroup=val;
}

