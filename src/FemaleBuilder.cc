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

#include"FemaleBuilder.hh"
#include "VBodyFactory.hh"

FemaleBuilder::FemaleBuilder()
{  
}

FemaleBuilder::~FemaleBuilder()

{ 
  delete body;
} 

void FemaleBuilder::BuildUterus(const G4String& colourName, G4bool solidVis, G4bool sensitivity )

{
if (trunkVolume == 0)
  G4Exception("FemaleBuilder::BuildUterus()", "human_phantom0006", FatalException, "The trunk volume is missing !!!!!");

 body -> CreateOrgan("Uterus",trunkVolume, colourName, solidVis, sensitivity);
 volOfUterus=body->GetBodyVOL("Uterus");
 rhoOfUterus=body->GetBodyRHO("Uterus");
}

void FemaleBuilder::BuildLeftOvary(const G4String& colourName, G4bool solidVis, G4bool sensitivity )

{
if (trunkVolume == 0)
  G4Exception("FemaleBuilder::BuildLeftOvary()", "human_phantom0007", FatalException, "The trunk volume is missing !!!!!");

  body -> CreateOrgan("LeftOvary",trunkVolume, colourName, solidVis, sensitivity);  
  volOfLeftOvary=body->GetBodyVOL("LeftOvary");
  rhoOfLeftOvary=body->GetBodyRHO("LeftOvary");
}

void FemaleBuilder::BuildRightOvary(const G4String& colourName, G4bool solidVis, G4bool sensitivity )

{
 if (trunkVolume == 0)
   G4Exception("FemaleBuilder::BuildRightOvary()", "human_phantom0008", FatalException, "The trunk volume is missing !!!!!");

  body -> CreateOrgan("RightOvary",trunkVolume, colourName, solidVis, sensitivity);  
  volOfRightOvary=body->GetBodyVOL("RightOvary");
  rhoOfRightOvary=body->GetBodyRHO("RightOvary");
}

void FemaleBuilder::BuildLeftBreast(const G4String& colourName,
				      G4bool solidVis, G4bool sensitivity)
{ 
  if (motherVolume == 0)
    G4Exception("FemaleBuilder::BuildLeftBreast()", "human_phantom0009", FatalException, "The world volume is missing !!!!!");

  body -> CreateOrgan("LeftBreast",motherVolume, colourName, solidVis, sensitivity);  
  volOfLeftBreast=body->GetBodyVOL("LeftBreast");
  rhoOfLeftBreast=body->GetBodyRHO("LeftBreast");

}

void FemaleBuilder::BuildRightBreast(const G4String& colourName, 
				       G4bool solidVis, 
				       G4bool sensitivity)
{
 if (motherVolume == 0)
   G4Exception("FemaleBuilder::BuildRightBreast()", "human_phantom0010", FatalException, "The world volume is missing !!!!!");

  body -> CreateOrgan("RightBreast",motherVolume, colourName,  solidVis, sensitivity);  
  volOfRightBreast=body->GetBodyVOL("RightBreast");
  rhoOfRightBreast=body->GetBodyRHO("RightBreast");

}
