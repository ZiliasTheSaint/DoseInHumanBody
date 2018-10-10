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
#include "MIRDFacial.hh"

#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "HumanPhantomMaterial.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4VSolid.hh"
#include "G4SubtractionSolid.hh"
#include "HumanPhantomColour.hh"
#include "MirdMapConstants.hh"

MIRDFacial::MIRDFacial()
{
	VOL=54.7*cm3;//init
	
	G4double d = 1.4862*g/cm3;//skeleton
	RHO=d;
}

MIRDFacial::~MIRDFacial()
{
  ;
}

G4VPhysicalVolume* MIRDFacial::Construct(const G4String& volumeName, G4VPhysicalVolume* mother, 
					      const G4String& colourName, G4bool wireFrame,G4bool sensitivity)
{

  MirdMapConstants* mmc = new MirdMapConstants();
  mmc->initFacial();
  std::map<std::string,facial_struct> esp = mmc->GetFacialMap();    
  VOL = scaleXY*scaleXY*scaleZ*esp[ageGroup].volume; 
  return 0;

}
