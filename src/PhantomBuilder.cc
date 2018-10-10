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
#include "PhantomBuilder.hh"
#include "VBodyFactory.hh"
#include "MIRDBodyFactory.hh"

PhantomBuilder::PhantomBuilder(): model("MIRD")
{  
  // sex can be "female" or "male"
  body = 0;
  motherVolume = 0;
  headVolume = 0;
  trunkVolume = 0;
  leftLegVolume =0;
  rightLegVolume =0;
  maleGenitaliaVolume = 0;  

  scaleXY=1.0;
  scaleZ=1.0;
  ageGroup="newborn";//init
}

PhantomBuilder::~PhantomBuilder()
{
} 
void PhantomBuilder::BuildTrunk(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (motherVolume == 0)
    G4Exception("PhantomBuilder::BuildTrunk()", "human_phantom0014", FatalException, "The world volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  motherVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  trunkVolume = body -> CreateOrgan("Trunk", motherVolume, colourName, solidVis, sensitivity);
  volOfTrunk=body->GetBodyVOL("Trunk");
  rhoOfTrunk=body->GetBodyRHO("Trunk");

}

void PhantomBuilder::BuildLeftLeg(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (motherVolume == 0)
    G4Exception("PhantomBuilder::BuildLeftLeg()", "human_phantom0015", FatalException, "The world volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  motherVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  leftLegVolume = body -> CreateOrgan("LeftLeg", motherVolume, colourName, solidVis, sensitivity);
  volOfLeftLeg=body->GetBodyVOL("LeftLeg");
  rhoOfLeftLeg=body->GetBodyRHO("LeftLeg");

}

void PhantomBuilder::BuildRightLeg(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (motherVolume == 0)
    G4Exception("PhantomBuilder::BuildRightLeg()", "human_phantom0016", FatalException, "The world volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  motherVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  rightLegVolume = body -> CreateOrgan("RightLeg", motherVolume, colourName, solidVis, sensitivity);
  volOfRightLeg=body->GetBodyVOL("RightLeg");
  rhoOfRightLeg=body->GetBodyRHO("RightLeg");

}

void PhantomBuilder::BuildLeftLegBone(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (leftLegVolume == 0)
    G4Exception("PhantomBuilder::BuildLeftLegBone()", "human_phantom0017", FatalException, "The left leg volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  leftLegVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan("LeftLegBone", leftLegVolume,colourName, solidVis, sensitivity);
  volOfLeftLegBone=body->GetBodyVOL("LeftLegBone");
  rhoOfLeftLegBone=body->GetBodyRHO("LeftLegBone");

}

void PhantomBuilder::BuildRightLegBone(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("PhantomBuilder::BuildRightLegBone()", "human_phantom0018", FatalException, "The right leg volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " << rightLegVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan("RightLegBone", rightLegVolume, colourName, solidVis, sensitivity);
  volOfRightLegBone=body->GetBodyVOL("RightLegBone");
  rhoOfRightLegBone=body->GetBodyRHO("RightLegBone");

}

void PhantomBuilder::BuildLeftArmBone(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("PhantomBuilder::BuildLeftArmBone()", "human_phantom0019", FatalException, "The world volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  trunkVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan("LeftArmBone" ,trunkVolume,colourName,solidVis, sensitivity);
  volOfLeftArmBone=body->GetBodyVOL("LeftArmBone");
  rhoOfLeftArmBone=body->GetBodyRHO("LeftArmBone");

}

void PhantomBuilder::BuildRightArmBone(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("PhantomBuilder::BuildRightArmBone()", "human_phantom0020", FatalException, "The trunk volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  trunkVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan("RightArmBone",trunkVolume,colourName,solidVis, sensitivity);
  volOfRightArmBone=body->GetBodyVOL("RightArmBone");
  rhoOfRightArmBone=body->GetBodyRHO("RightArmBone");

}

void PhantomBuilder::BuildLeftScapula(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("G4PhantomBuilder::BuildLeftScapula()", "human_phantom0021", FatalException, "The trunk volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  trunkVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan("LeftScapula",trunkVolume,colourName,solidVis, sensitivity);
  volOfLeftScapula=body->GetBodyVOL("LeftScapula");
  rhoOfLeftScapula=body->GetBodyRHO("LeftScapula");

}

void PhantomBuilder::BuildRightScapula(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("PhantomBuilder::BuildRightScapula()", "human_phantom0022", FatalException, "The trunk volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  trunkVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan("RightScapula",trunkVolume,colourName,solidVis, sensitivity);
  volOfRightScapula=body->GetBodyVOL("RightScapula");
  rhoOfRightScapula=body->GetBodyRHO("RightScapula");

}

void PhantomBuilder::BuildHead(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (motherVolume == 0)
    G4Exception("PhantomBuilder::BuildHead()", "human_phantom0023", FatalException, "The mother volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  motherVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  headVolume = body -> CreateOrgan("Head",motherVolume, colourName, solidVis, sensitivity);
  volOfHead=body->GetBodyVOL("Head");
  rhoOfHead=body->GetBodyRHO("Head");
}

void PhantomBuilder::BuildFacial(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (headVolume == 0)
    G4Exception("PhantomBuilder::BuildSkull()", "human_phantom0024", FatalException, "The head volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  headVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan( "FacialBones",headVolume, colourName, solidVis, sensitivity);
  volOfFacial=body->GetBodyVOL("FacialBones");
  rhoOfFacial=body->GetBodyRHO("FacialBones");
}

void PhantomBuilder::BuildSkull(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (headVolume == 0)
    G4Exception("PhantomBuilder::BuildSkull()", "human_phantom0024", FatalException, "The head volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  headVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan( "Skull",headVolume, colourName, solidVis, sensitivity);
  volOfSkull=body->GetBodyVOL("Skull");
  rhoOfSkull=body->GetBodyRHO("Skull");
}

void PhantomBuilder::BuildUpperSpine(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (headVolume == 0)
    G4Exception("PhantomBuilder::BuildUpperSpine()", "human_phantom0025", FatalException, "The head volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  headVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan("UpperSpine",headVolume,colourName, solidVis, sensitivity);
  volOfUpperSpine=body->GetBodyVOL("UpperSpine");
  rhoOfUpperSpine=body->GetBodyRHO("UpperSpine");

}

void PhantomBuilder::BuildMiddleLowerSpine(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("PhantomBuilder::BuildMiddleLowerSpine()", "human_phantom0026", FatalException, "The trunk volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  trunkVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan("MiddleLowerSpine",trunkVolume, colourName, solidVis, sensitivity);
  volOfMiddleLowerSpine=body->GetBodyVOL("MiddleLowerSpine");
  rhoOfMiddleLowerSpine=body->GetBodyRHO("MiddleLowerSpine");

}

void PhantomBuilder::BuildPelvis(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("PhantomBuilder::BuildPelvis()", "human_phantom0027", FatalException, "The trunk volume is missing !!!!!");

    body -> CreateOrgan( "Pelvis",trunkVolume, colourName, solidVis, sensitivity);
  volOfPelvis=body->GetBodyVOL("Pelvis");
  rhoOfPelvis=body->GetBodyRHO("Pelvis");

}



void PhantomBuilder::BuildLeftClavicle(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("PhantomBuilder", "human_phantom0028", FatalException, "The trunk volume is missing !!!!!");

   body -> CreateOrgan( "LeftClavicle",trunkVolume, colourName, solidVis, sensitivity);
   volOfLeftClavicle=body->GetBodyVOL("LeftClavicle");
   rhoOfLeftClavicle=body->GetBodyRHO("LeftClavicle");
}
void PhantomBuilder::BuildRightClavicle(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("PhantomBuilder", "human_phantom0028", FatalException, "The trunk volume is missing !!!!!");

   body -> CreateOrgan( "RightClavicle",trunkVolume, colourName, solidVis, sensitivity);
   volOfRightClavicle=body->GetBodyVOL("RightClavicle");
   rhoOfRightClavicle=body->GetBodyRHO("RightClavicle");
}


void PhantomBuilder::BuildBrain(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
 if (headVolume == 0)
   G4Exception("PhantomBuilder::BuildBrain()", "human_phantom0029", FatalException, "The head volume is missing !!!!!");

    body -> CreateOrgan("Brain",headVolume, colourName, solidVis, sensitivity);
	volOfBrain=body->GetBodyVOL("Brain");
    rhoOfBrain=body->GetBodyRHO("Brain");
}

void PhantomBuilder::BuildHeart(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
    if (trunkVolume == 0)
      G4Exception("PhantomBuilder::BuildHeart()", "human_phantom0030", FatalException, "The trunk volume is missing !!!!!");
    body -> CreateOrgan("Heart", trunkVolume,colourName, solidVis, sensitivity);
  volOfHeart=body->GetBodyVOL("Heart");
  rhoOfHeart=body->GetBodyRHO("Heart");

}

void PhantomBuilder::BuildLeftLung(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
   if (trunkVolume == 0)
     G4Exception("PhantomBuilder::BuildLeftLung()", "human_phantom0031", FatalException, "The trunk volume is missing !!!!!");

    body -> CreateOrgan("LeftLung",trunkVolume,colourName,solidVis, sensitivity);
  volOfLeftLung=body->GetBodyVOL("LeftLung");
  rhoOfLeftLung=body->GetBodyRHO("LeftLung");

}

void PhantomBuilder::BuildRightLung(const G4String& colourName, G4bool solidVis, G4bool sensitivity )
{ 
   if (trunkVolume == 0)
     G4Exception("PhantomBuilder::BuildRightLung()", "human_phantom0032", FatalException, "The trunk volume is missing !!!!!");

    body -> CreateOrgan("RightLung",trunkVolume,colourName, solidVis, sensitivity);
  volOfRightLung=body->GetBodyVOL("RightLung");
  rhoOfRightLung=body->GetBodyRHO("RightLung");
  
}

void PhantomBuilder::BuildStomach(const G4String& colourName, G4bool solidVis, G4bool sensitivity )
{ 
  if (trunkVolume == 0)
    G4Exception("PhantomBuilder::BuildStomach()", "human_phantom0033", FatalException, "The trunk volume is missing !!!!!");

    body -> CreateOrgan("Stomach",trunkVolume,colourName, solidVis, sensitivity);
  volOfStomach=body->GetBodyVOL("Stomach");
  rhoOfStomach=body->GetBodyRHO("Stomach");

}

void PhantomBuilder::BuildRibCage(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
   if (trunkVolume == 0)
     G4Exception("PhantomBuilder::BuildRibCage()", "human_phantom0034", FatalException, "The trunk volume is missing !!!!!");

    body -> CreateOrgan("RibCage",trunkVolume,colourName, solidVis, sensitivity);
  volOfRibCage=body->GetBodyVOL("RibCage");
  rhoOfRibCage=body->GetBodyRHO("RibCage");

}

void PhantomBuilder::BuildSpleen(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
   if (trunkVolume == 0)
     G4Exception("PhantomBuilder::BuildSpleen()", "human_phantom0035", FatalException, "The trunk volume is missing !!!!!");

    body -> CreateOrgan("Spleen", trunkVolume,colourName, solidVis, sensitivity);
  volOfSpleen=body->GetBodyVOL("Spleen");
  rhoOfSpleen=body->GetBodyRHO("Spleen");

}

void PhantomBuilder::BuildUpperLargeIntestine(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
   if (trunkVolume == 0)
     G4Exception("PhantomBuilder::BuildUpperLargeIntestine()", "human_phantom0036", FatalException, "The trunk volume is missing !!!!!");

    body -> CreateOrgan("UpperLargeIntestine",trunkVolume, colourName, solidVis, sensitivity);
  volOfUpperLargeIntestine=body->GetBodyVOL("UpperLargeIntestine");
  rhoOfUpperLargeIntestine=body->GetBodyRHO("UpperLargeIntestine");

}

void PhantomBuilder::BuildLowerLargeIntestine(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("PhantomBuilder::BuildLowerLargeIntestine()", "human_phantom0037", FatalException, "The trunk volume is missing !!!!!");

   body -> CreateOrgan("LowerLargeIntestine", trunkVolume, colourName,solidVis, sensitivity);
  volOfLowerLargeIntestine=body->GetBodyVOL("LowerLargeIntestine");
  rhoOfLowerLargeIntestine=body->GetBodyRHO("LowerLargeIntestine");

}

void PhantomBuilder::BuildLeftKidney(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
 if (trunkVolume == 0)
   G4Exception("PhantomBuilder::BuildLeftKidney()", "human_phantom0039", FatalException, "The trunk volume is missing !!!!!");

    body -> CreateOrgan("LeftKidney", trunkVolume,colourName, solidVis, sensitivity);
  volOfLeftKidney=body->GetBodyVOL("LeftKidney");
  rhoOfLeftKidney=body->GetBodyRHO("LeftKidney");

}

void PhantomBuilder::BuildRightKidney(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
   if (trunkVolume == 0)
     G4Exception("PhantomBuilder::BuildRightKidney()", "human_phantom0040", FatalException, "The trunk volume is missing !!!!!");

    body -> CreateOrgan("RightKidney",trunkVolume,colourName, solidVis, sensitivity);
  volOfRightKidney=body->GetBodyVOL("RightKidney");
  rhoOfRightKidney=body->GetBodyRHO("RightKidney");

}

void PhantomBuilder::BuildLeftAdrenal(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("PhantomBuilder::BuildLeftAdrenal()", "human_phantom0041", FatalException, "The trunk volume is missing !!!!!");

    body -> CreateOrgan("LeftAdrenal", trunkVolume,colourName, solidVis, sensitivity);
  volOfLeftAdrenal=body->GetBodyVOL("LeftAdrenal");
  rhoOfLeftAdrenal=body->GetBodyRHO("LeftAdrenal");

}

void PhantomBuilder::BuildRightAdrenal(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("PhantomBuilder::BuildRightAdrenal()", "human_phantom0042", FatalException, "The trunk volume is missing !!!!!");

    body -> CreateOrgan("RightAdrenal", trunkVolume,colourName, solidVis, sensitivity);
  volOfRightAdrenal=body->GetBodyVOL("RightAdrenal");
  rhoOfRightAdrenal=body->GetBodyRHO("RightAdrenal");

}

void PhantomBuilder::BuildLiver(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("PhantomBuilder::BuildLiver()", "human_phantom0043", FatalException, "The trunk volume is missing !!!!!");

    body -> CreateOrgan("Liver", trunkVolume,colourName, solidVis, sensitivity);
  volOfLiver=body->GetBodyVOL("Liver");
  rhoOfLiver=body->GetBodyRHO("Liver");

}

void PhantomBuilder::BuildPancreas(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
   if (trunkVolume == 0)
     G4Exception("PhantomBuilder::BuildPancreas()", "human_phantom0044", FatalException, "The trunk volume is missing !!!!!");

    body -> CreateOrgan("Pancreas",trunkVolume,colourName, solidVis, sensitivity);
  volOfPancreas=body->GetBodyVOL("Pancreas");
  rhoOfPancreas=body->GetBodyRHO("Pancreas");

}

void PhantomBuilder::BuildUrinaryBladder(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (trunkVolume == 0)
    G4Exception("PhantomBuilder::BuildUrinaryBladder()", "human_phantom0045", FatalException, "The trunk volume is missing !!!!!");

    body -> CreateOrgan("UrinaryBladder",trunkVolume, colourName, solidVis, sensitivity);
  volOfUrinaryBladder=body->GetBodyVOL("UrinaryBladder");
  rhoOfUrinaryBladder=body->GetBodyRHO("UrinaryBladder");

}

void PhantomBuilder::BuildThyroid(const G4String& colourName, G4bool solidVis, G4bool sensitivity )
{ 
   if (headVolume == 0)
     G4Exception("PhantomBuilder::BuildThyroid()", "human_phantom0046", FatalException, "The trunk volume is missing !!!!!");

   body -> CreateOrgan("Thyroid",headVolume, colourName,solidVis, sensitivity);
  volOfThyroid=body->GetBodyVOL("Thyroid");
  rhoOfThyroid=body->GetBodyRHO("Thyroid");

}

void PhantomBuilder::BuildEsophagus(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("PhantomBuilder", "human_phantom0028", FatalException, "The trunk volume is missing !!!!!");

   body -> CreateOrgan( "Esophagus",trunkVolume, colourName, solidVis, sensitivity);
   volOfEsophagus=body->GetBodyVOL("Esophagus");
   rhoOfEsophagus=body->GetBodyRHO("Esophagus");
}

void PhantomBuilder::BuildGallBladder(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("PhantomBuilder", "human_phantom0028", FatalException, "The trunk volume is missing !!!!!");

   body -> CreateOrgan( "GallBladder",trunkVolume, colourName, solidVis, sensitivity);
   volOfGallBladder=body->GetBodyVOL("GallBladder");
   rhoOfGallBladder=body->GetBodyRHO("GallBladder");
}

void PhantomBuilder::BuildSmallIntestine(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("PhantomBuilder", "human_phantom0028", FatalException, "The trunk volume is missing !!!!!");

   body -> CreateOrgan( "SmallIntestine",trunkVolume, colourName, solidVis, sensitivity);
   volOfSmallIntestine=body->GetBodyVOL("SmallIntestine");
   rhoOfSmallIntestine=body->GetBodyRHO("SmallIntestine");
}

void PhantomBuilder::BuildThymus(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
   if (trunkVolume == 0)
   G4Exception("PhantomBuilder", "human_phantom0028", FatalException, "The trunk volume is missing !!!!!");

   body -> CreateOrgan( "Thymus",trunkVolume, colourName, solidVis, sensitivity);
   volOfThymus=body->GetBodyVOL("Thymus");
   rhoOfThymus=body->GetBodyRHO("Thymus");
}

G4VPhysicalVolume* PhantomBuilder::GetPhantom()
{
  return motherVolume;
}

void PhantomBuilder::SetMotherVolume(G4VPhysicalVolume* mother)
{
  motherVolume = mother;
}


void PhantomBuilder::SetModel(G4String modelFlag)
{
  model = modelFlag;

  if(model=="MIRD"){
	  body = new MIRDBodyFactory();
	  body->SetScaleXY(scaleXY);
	  body->SetScaleZ(scaleZ);
	  body->SetAgeGroup(ageGroup);
  }
}

void PhantomBuilder::SetScaleXY(G4double val){
	scaleXY=val;
}

void PhantomBuilder::SetScaleZ(G4double val){
	scaleZ=val;
}

void PhantomBuilder::SetAgeGroup(G4String val){
	ageGroup=val;
}