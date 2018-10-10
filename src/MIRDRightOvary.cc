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
#include "MIRDRightOvary.hh"
#include "G4SDManager.hh"

#include "globals.hh"

#include "G4VisAttributes.hh"
#include "G4Ellipsoid.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "HumanPhantomMaterial.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UnionSolid.hh"
#include "HumanPhantomColour.hh"

#include "MirdMapConstants.hh"

MIRDRightOvary::MIRDRightOvary()
{
}

MIRDRightOvary::~MIRDRightOvary()
{

}

G4VPhysicalVolume* MIRDRightOvary::Construct(const G4String& volumeName, G4VPhysicalVolume* mother, 
						    const G4String& colourName, G4bool wireFrame,G4bool sensitivity)
{ 
  G4cout << "Construct "<< volumeName << G4endl;
 
 HumanPhantomMaterial* material = new HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;
 
  MirdMapConstants* mmc = new MirdMapConstants();  
  mmc->initTrunk();
  std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  G4double ct = scaleZ*trk[ageGroup].Ct;//our zcenter is MIRD 0.5*ct z coordinate
  mmc->initOvaries();
  std::map<std::string,ovaries_struct> ovr = mmc->GetOvariesMap();  
  G4double a = scaleXY*ovr[ageGroup].a;
  G4double b = scaleXY*ovr[ageGroup].b;
  G4double c = scaleZ*ovr[ageGroup].c;
  G4double x0 = scaleXY*ovr[ageGroup].x0;  
  G4double z0 = scaleZ*ovr[ageGroup].z0;  
  delete mmc;

 G4double ax= a;//scaleXY*1. *cm;
 G4double by= b;//scaleXY*0.5*cm;
 G4double cz= c;//scaleZ*2.*cm;

 G4Ellipsoid* OneOvary = new G4Ellipsoid("OneOvary",
					 ax, by, cz);

  G4LogicalVolume* logicRightOvary = new G4LogicalVolume(OneOvary,
						    soft,
						    "logical" + volumeName,
						    0, 0, 0);
  
  // Define rotation and position here!
  G4VPhysicalVolume* physRightOvary = new G4PVPlacement(0,
			       //G4ThreeVector(scaleXY*6. *cm,0.0*cm, -scaleZ*20*cm),
				   G4ThreeVector(-x0,0.0*cm, 0.5*(0.0-ct)+z0),
      			       "physicalRightOvary",
  			       logicRightOvary,
			       mother,
			       false,
			       0, checkOverlaps);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicRightOvary->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  // G4VisAttributes* RightOvaryVisAtt = new G4VisAttributes(G4Colour(0.85,0.44,0.84));
  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* RightOvaryVisAtt = new G4VisAttributes(colour);
  RightOvaryVisAtt->SetForceSolid(wireFrame);
  logicRightOvary->SetVisAttributes(RightOvaryVisAtt);

  G4cout << "RightOvary created !!!!!!" << G4endl;

  // Testing RightOvary Volume
  G4double RightOvaryVol = logicRightOvary->GetSolid()->GetCubicVolume();
  G4cout << "Volume of RightOvary = " << RightOvaryVol/cm3 << " cm^3" << G4endl;
  
  // Testing RightOvary Material
  G4String RightOvaryMat = logicRightOvary->GetMaterial()->GetName();
  G4cout << "Material of RightOvary = " << RightOvaryMat << G4endl;
  
  // Testing Density
  G4double RightOvaryDensity = logicRightOvary->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << RightOvaryDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double RightOvaryMass = (RightOvaryVol)*RightOvaryDensity;
  G4cout << "Mass of RightOvary = " << RightOvaryMass/gram << " g" << G4endl;
  
  VOL=RightOvaryVol;RHO=RightOvaryDensity;
  return physRightOvary;
}
