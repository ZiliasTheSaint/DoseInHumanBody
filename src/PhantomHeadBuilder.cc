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
#include "PhantomHeadBuilder.hh"
#include "VBodyFactory.hh"
#include "MIRDBodyFactory.hh"
#include "G4PVPlacement.hh"

PhantomHeadBuilder::PhantomHeadBuilder(): model("MIRD")
{  
  // sex can be "female" or "male"
  body = 0;
  motherVolume = 0;
  headVolume = 0;
}

PhantomHeadBuilder::~PhantomHeadBuilder()
{
  delete body;
} 

void PhantomHeadBuilder::BuildHead(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (motherVolume == 0)
    G4Exception("G4PhantomHeadBuilder::BuildHead()", "human_phantom0011", FatalException, "The moder Volume volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  motherVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  headVolume = body -> CreateOrgan("Head",motherVolume, colourName, solidVis, sensitivity);
  volOfHead=body->GetBodyVOL("Head");
  rhoOfHead=body->GetBodyRHO("Head");
}

G4double PhantomHeadBuilder::getMassOfHead(){
	G4double d = (volOfHead-volOfSkull-volOfBrain-volOfFacial)*rhoOfHead;
	return d;
}

void PhantomHeadBuilder::BuildFacial(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (headVolume == 0)
    G4Exception("PhantomBuilder::BuildSkull()", "human_phantom0024", FatalException, "The head volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  headVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan( "FacialBones",headVolume, colourName, solidVis, sensitivity);
  volOfFacial=body->GetBodyVOL("FacialBones");
  rhoOfFacial=body->GetBodyRHO("FacialBones");
}

void PhantomHeadBuilder::BuildSkull(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (headVolume == 0)
    G4Exception("G4PhantomHeadBuilder::BuildSkull()", "human_phantom0012", FatalException, "The head volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " <<  headVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  body -> CreateOrgan( "Skull",headVolume, colourName, solidVis, sensitivity);
  volOfSkull=body->GetBodyVOL("Skull");
  rhoOfSkull=body->GetBodyRHO("Skull");
}

G4double PhantomHeadBuilder::getMassOfSkull(){
	G4double d = volOfSkull*rhoOfSkull;
	return d;
}

void PhantomHeadBuilder::BuildBrain(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
 if (headVolume == 0)
   G4Exception("G4PhantomHeadBuilder::BuildBrain()", "human_phantom0013", FatalException, "The head volume is missing !!!!!");

    body -> CreateOrgan("Brain",headVolume, colourName, solidVis, sensitivity);
	volOfBrain=body->GetBodyVOL("Brain");
    rhoOfBrain=body->GetBodyRHO("Brain");
}

G4double PhantomHeadBuilder::getMassOfBrain(){
	G4double d = volOfBrain*rhoOfBrain;
	return d;
}

G4VPhysicalVolume* PhantomHeadBuilder::GetPhantom()
{
  return motherVolume;
}

void PhantomHeadBuilder::SetMotherVolume(G4VPhysicalVolume* mother)
{
  motherVolume = mother;
}


void PhantomHeadBuilder::SetModel(G4String modelFlag)
{
  model = modelFlag;

  if(model=="MIRD") body = new MIRDBodyFactory(); 

}

