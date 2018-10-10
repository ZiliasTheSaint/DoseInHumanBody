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

#include "MIRDRightArmBone.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "HumanPhantomMaterial.hh"
#include "G4EllipticalTube.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UnionSolid.hh"
#include "G4EllipticalCone.hh"
#include "HumanPhantomColour.hh"

#include "MirdMapConstants.hh"

MIRDRightArmBone::MIRDRightArmBone()
{
}

MIRDRightArmBone::~MIRDRightArmBone()
{
}

G4VPhysicalVolume* MIRDRightArmBone::Construct(const G4String& volumeName,G4VPhysicalVolume* mother,
						     const G4String& colourName,G4bool wireFrame, G4bool sensitivity)
{
  // Remind! the elliptical cone gives problems! Intersections of volumes, 
  // wrong calculation of the volume!
   
  HumanPhantomMaterial* material = new HumanPhantomMaterial();   
  G4cout << "Construct " << volumeName << G4endl;   
  G4Material* skeleton = material -> GetMaterial("skeleton");  
  delete material;

  MirdMapConstants* mmc = new MirdMapConstants();
  mmc->initTrunk();
  std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  G4double ct = scaleZ*trk[ageGroup].Ct;
  mmc->initArmBone();
  std::map<std::string,armbone_struct> arm = mmc->GetArmBoneMap();  
  G4double a = scaleXY*arm[ageGroup].a;
  G4double b = scaleXY*arm[ageGroup].b;
  G4double z2 =scaleZ*arm[ageGroup].z2;
  G4double x0 =scaleXY*arm[ageGroup].x0;

  G4double dx = a;//scaleXY*1.4 * cm;//a
  G4double dy = b;//scaleXY*2.7 * cm;//b
  
  //G4EllipticalTube* rightArm = new G4EllipticalTube("OneRightArmBone",dx,dy,scaleZ*34.5 *cm);
  G4EllipticalTube* rightArm = new G4EllipticalTube("OneRightArmBone",dx,dy,0.5*z2);
 
  G4LogicalVolume* logicRightArmBone = new G4LogicalVolume(rightArm,
						      skeleton,
						      "logical" + volumeName,
						      0, 0,0);

  G4RotationMatrix* matrix = new G4RotationMatrix();
  matrix -> rotateX(180. * degree);

  G4VPhysicalVolume* physRightArmBone = new G4PVPlacement(matrix,
			       //G4ThreeVector(-scaleXY*18.4 * cm, 0.0, -scaleZ*0.5*cm), //-x0
				   //G4ThreeVector(-x0, 0.0, 0.0),
				   G4ThreeVector(-x0, 0.0, 0.5*(z2-ct)),
      			       "physicalRightArmBone",
  			       logicRightArmBone,
			       mother,
			       false,0,checkOverlaps);
			      

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicRightArmBone->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }


  // Visualization Attributes
  //G4VisAttributes* RightArmBoneVisAtt = new G4VisAttributes(G4Colour(0.46,0.53,0.6));
 HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* RightArmBoneVisAtt = new G4VisAttributes(colour);

  RightArmBoneVisAtt->SetForceSolid(wireFrame);
  logicRightArmBone->SetVisAttributes(RightArmBoneVisAtt);

  G4cout << "RightArmBone created !!!!!!" << G4endl;
 
  // Testing RightArmBone Volume
  G4double RightArmBoneVol = logicRightArmBone->GetSolid()->GetCubicVolume();
  G4cout << "Volume of RightArmBone = " << RightArmBoneVol/cm3 << " cm^3" << G4endl;
  
  // Testing RightArmBone Material
  G4String RightArmBoneMat = logicRightArmBone->GetMaterial()->GetName();
  G4cout << "Material of RightArmBone = " << RightArmBoneMat << G4endl;
  
  // Testing Density
  G4double RightArmBoneDensity = logicRightArmBone->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << RightArmBoneDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double RightArmBoneMass = (RightArmBoneVol)*RightArmBoneDensity;
  G4cout << "Mass of RightArmBone = " << RightArmBoneMass/gram << " g" << G4endl;
  
  VOL=RightArmBoneVol;RHO=RightArmBoneDensity;
  return physRightArmBone;
}
