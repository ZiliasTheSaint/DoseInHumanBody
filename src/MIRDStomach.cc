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
#include "MIRDStomach.hh"

#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Ellipsoid.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "HumanPhantomMaterial.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "HumanPhantomColour.hh"

#include "MirdMapConstants.hh"

MIRDStomach::MIRDStomach()
{
}

MIRDStomach::~MIRDStomach()
{
}

G4VPhysicalVolume* MIRDStomach::Construct(const G4String& volumeName,G4VPhysicalVolume* mother,
					         const G4String& colourName, G4bool wireFrame, G4bool sensitivity)
{

  G4cout << "Construct "<< volumeName << G4endl;
 
 HumanPhantomMaterial* material = new HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;

 MirdMapConstants* mmc = new MirdMapConstants();  
  mmc->initTrunk();
  std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  G4double ct = scaleZ*trk[ageGroup].Ct;//our zcenter is MIRD 0.5*ct z coordinate
  mmc->initStomach();
  std::map<std::string,stomach_struct> stm = mmc->GetStomachMap();  
  G4double a = scaleXY*stm[ageGroup].a;
  G4double b = scaleXY*stm[ageGroup].b;
  G4double c = scaleZ*stm[ageGroup].c;
  G4double x0 = scaleXY*stm[ageGroup].x0;  
  G4double y0 = scaleXY*stm[ageGroup].y0;  
  G4double z0 = scaleZ*stm[ageGroup].z0;    
  delete mmc;

 G4double ax = a;//scaleXY*4. * cm;
 G4double by= b;//scaleXY*3. * cm;
 G4double cz = c;//scaleZ*8. * cm;

  G4Ellipsoid* stomach_out = new G4Ellipsoid("stomach_out", 
					 ax, by, cz);

  G4LogicalVolume* logicStomach = new G4LogicalVolume(stomach_out, soft,
						      "logical" + volumeName, 0, 0, 0);
  
  // Define rotation and position here!
  G4VPhysicalVolume* physStomach = new G4PVPlacement(0,
	  //G4ThreeVector(scaleXY*8. *cm,-scaleXY*4. * cm, 0),
	  G4ThreeVector(x0,y0, 0.5*(0.0-ct)+z0),//0),
      			       "physicalStomach",
  			       logicStomach,
			       mother,
			       false,
			       0, checkOverlaps);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicStomach->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);

   G4VisAttributes* StomachVisAtt = new G4VisAttributes(colour);
  StomachVisAtt->SetForceSolid(wireFrame);
  logicStomach->SetVisAttributes(StomachVisAtt);

  G4cout << "Stomach created !!!!!!" << G4endl;

  // Testing Stomach Volume
  G4double StomachVol = logicStomach->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Stomach = " << StomachVol/cm3 << " cm^3" << G4endl;
  
  // Testing Stomach Material
  G4String StomachMat = logicStomach->GetMaterial()->GetName();
  G4cout << "Material of Stomach = " << StomachMat << G4endl;
  
  // Testing Density
  G4double StomachDensity = logicStomach->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << StomachDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double StomachMass = (StomachVol)*StomachDensity;
  G4cout << "Mass of Stomach = " << StomachMass/gram << " g" << G4endl;
  
  VOL=StomachVol;RHO=StomachDensity;
  return physStomach;
}
