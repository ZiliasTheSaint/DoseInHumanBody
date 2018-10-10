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
#include "MIRDUrinaryBladder.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Ellipsoid.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "HumanPhantomMaterial.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "HumanPhantomColour.hh"

#include "MirdMapConstants.hh"

MIRDUrinaryBladder::MIRDUrinaryBladder()
{
}

MIRDUrinaryBladder::~MIRDUrinaryBladder()
{
}

G4VPhysicalVolume* MIRDUrinaryBladder::Construct(const G4String& volumeName, G4VPhysicalVolume* mother,
							const G4String& colourName, G4bool wireFrame, G4bool sensitivity)
{
 
  G4cout << "Construct " << volumeName << G4endl;
  
 HumanPhantomMaterial* material = new HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;

 MirdMapConstants* mmc = new MirdMapConstants();  
  mmc->initTrunk();
  std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  G4double ct = scaleZ*trk[ageGroup].Ct;//our zcenter is MIRD 0.5*ct z coordinate
  mmc->initUrinaryBladder();
  std::map<std::string,urinaryBladder_struct> urn = mmc->GetUrinaryBladderMap();  
  G4double a = scaleXY*urn[ageGroup].a;
  G4double b = scaleXY*urn[ageGroup].b;
  G4double c = scaleZ*urn[ageGroup].c;
  G4double y0 = scaleXY*urn[ageGroup].y0;  
  G4double z0 = scaleZ*urn[ageGroup].z0;    
  delete mmc;

 G4double ax = a;//scaleXY*4.958*cm; 
 G4double by= b;//scaleXY*3.458 *cm;
 G4double cz= c;//scaleZ*3.458 *cm;
 
 G4Ellipsoid* bladder = new G4Ellipsoid("bladder_out",ax, by, cz);
 /*
 ax = scaleXY*4.706 * cm;
 by = scaleXY*3.206 * cm;
 cz = scaleZ*3.206 * cm;
 G4Ellipsoid* inner = new G4Ellipsoid("innerBladder", ax, by, cz);
 
 G4SubtractionSolid* totalBladder = new G4SubtractionSolid("bladder", bladder, inner);
 *///want total +content not just wall=>error by INFN

 G4LogicalVolume* logicUrinaryBladder = new G4LogicalVolume(bladder//totalBladder
							    , soft,"logical" + volumeName,
							    0, 0, 0);
  
  // Define rotation and position here!
  G4VPhysicalVolume* physUrinaryBladder = new G4PVPlacement(0,
	  //G4ThreeVector(0 *cm, -scaleXY*4.5 *cm,-scaleZ*27. *cm),
	  G4ThreeVector(0 *cm, y0,0.5*(0.0-ct)+z0),//-scaleZ*27. *cm),//y0<0
      			       "physicalUrinaryBladder",
  			       logicUrinaryBladder,
			       mother,
			       false,
			       0, checkOverlaps);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicUrinaryBladder->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* UrinaryBladderVisAtt = new G4VisAttributes(colour);
  //G4VisAttributes* UrinaryBladderVisAtt = new G4VisAttributes(G4Colour(0.85,0.65,0.125));

  UrinaryBladderVisAtt->SetForceSolid(wireFrame);
  logicUrinaryBladder->SetVisAttributes(UrinaryBladderVisAtt);

  G4cout << "UrinaryBladder created !!!!!!" << G4endl;

  // Testing UrinaryBladder Volume
  G4double UrinaryBladderVol = logicUrinaryBladder->GetSolid()->GetCubicVolume();
  G4cout << "Volume of UrinaryBladder = " << UrinaryBladderVol/cm3 << " cm^3" << G4endl;
  
  // Testing UrinaryBladder Material
  G4String UrinaryBladderMat = logicUrinaryBladder->GetMaterial()->GetName();
  G4cout << "Material of UrinaryBladder = " << UrinaryBladderMat << G4endl;
  
  // Testing Density
  G4double UrinaryBladderDensity = logicUrinaryBladder->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << UrinaryBladderDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double UrinaryBladderMass = (UrinaryBladderVol)*UrinaryBladderDensity;
  G4cout << "Mass of UrinaryBladder = " << UrinaryBladderMass/gram << " g" << G4endl;

  VOL=UrinaryBladderVol;RHO=UrinaryBladderDensity;
  return physUrinaryBladder;
}
