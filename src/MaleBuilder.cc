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
#include"MaleBuilder.hh"
#include "VBodyFactory.hh"

MaleBuilder::MaleBuilder()
{  
}

MaleBuilder::~MaleBuilder()
{
  delete body;
} 

void MaleBuilder::BuildLeftTesticle(const G4String& colourName, G4bool solidVis, G4bool sensitivity )

{
 body -> CreateOrgan("LeftTesticle",trunkVolume, colourName, solidVis, sensitivity);
 volOfLeftTesticle=body->GetBodyVOL("LeftTesticle");
 rhoOfLeftTesticle=body->GetBodyRHO("LeftTesticle");
}

void MaleBuilder::BuildRightTesticle(const G4String& colourName, G4bool solidVis, G4bool sensitivity )

{
//if (trunkVolume == 0)
  //G4Exception("FemaleBuilder::BuildUterus()", "human_phantom0006", FatalException, "The trunk volume is missing !!!!!");

 body -> CreateOrgan("RightTesticle",trunkVolume, colourName, solidVis, sensitivity);
 volOfRightTesticle=body->GetBodyVOL("RightTesticle");
 rhoOfRightTesticle=body->GetBodyRHO("RightTesticle");
}
