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

#include "MIRDRightLegBone.hh"

#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "HumanPhantomMaterial.hh"
#include "G4EllipticalTube.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Cons.hh"
#include "G4UnionSolid.hh"
#include "HumanPhantomColour.hh"

#include "MirdMapConstants.hh"

MIRDRightLegBone::MIRDRightLegBone()
{
}

MIRDRightLegBone::~MIRDRightLegBone()
{
}

G4VPhysicalVolume* MIRDRightLegBone::Construct(const G4String& volumeName,G4VPhysicalVolume* mother,
                                                      const G4String& colourName, G4bool wireFrame,G4bool sensitivity)
{
 
  HumanPhantomMaterial* material = new HumanPhantomMaterial();
  
  G4cout << "Construct " << volumeName << G4endl;
   
  G4Material* skeleton = material -> GetMaterial("skeleton");
  
  delete material;
 
  MirdMapConstants* mmc = new MirdMapConstants();
  mmc->initTrunk();
  std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  G4double at = scaleXY*trk[ageGroup].At;  
  mmc->initLegs();
  std::map<std::string,legs_struct> lgs = mmc->GetLegsMap(); 
  G4double cl=scaleZ*lgs[ageGroup].Cl;
  G4double cpl=scaleZ*lgs[ageGroup].Cpl;
  mmc->initSkin();
  std::map<std::string,skin_struct> skn = mmc->GetSkinMap();
  G4double s = scaleXY*skn[ageGroup].S;//scale!!??!!fat=>more skin!
  delete mmc;
  G4double R1=0.175*at;//adult=3.5
  G4double R2=(at/4.0)*(cpl-cl)/cpl;//adult=1

  G4double dz = cl-s;//scaleZ*79.8 * cm;
  G4double rmin1 = 0.0 * cm;
  G4double rmin2 = 0.0 * cm;
  G4double rmax1 = R2;//scaleXY*1. * cm;//mird z=-80
  G4double rmax2 = R1;//scaleXY*3.5 * cm;//mird z=0
  G4double startphi = 0. * degree;
  G4double deltaphi = 360. * degree;

  G4Cons* leg_bone = new G4Cons("OneRightLegBone",  
			   rmin1, rmax1, 
			   rmin2, rmax2, dz/2., 
			   startphi, deltaphi);

  
  G4LogicalVolume* logicRightLegBone = new G4LogicalVolume(leg_bone, skeleton,"logical" + volumeName,
						      0, 0, 0);


  // Define rotation and position here!
  G4VPhysicalVolume* physRightLegBone = new G4PVPlacement(0,
				//G4ThreeVector(0.0 * cm, 0.0, scaleZ*0.1*cm),
				G4ThreeVector(0.0, 0.0, 0.5*s),
      			       "physicalRightLegBone",
  			       logicRightLegBone,
			       mother,
			       false,
			       0, checkOverlaps);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicRightLegBone->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* RightLegBoneVisAtt = new G4VisAttributes(colour);

  RightLegBoneVisAtt->SetForceSolid(wireFrame);
  logicRightLegBone->SetVisAttributes(RightLegBoneVisAtt);

  G4cout << "RightLegBone created !!!!!!" << G4endl;

  // Testing RightLegBone Volume
  G4double RightLegBoneVol = logicRightLegBone->GetSolid()->GetCubicVolume();
  G4cout << "Volume of RightLegBone = " << RightLegBoneVol/cm3 << " cm^3" << G4endl;
  
  // Testing RightLegBone Material
  G4String RightLegBoneMat = logicRightLegBone->GetMaterial()->GetName();
  G4cout << "Material of RightLegBone = " << RightLegBoneMat << G4endl;
  
  // Testing Density
  G4double RightLegBoneDensity = logicRightLegBone->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << RightLegBoneDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double RightLegBoneMass = (RightLegBoneVol)*RightLegBoneDensity;
  G4cout << "Mass of RightLegBone = " << RightLegBoneMass/gram << " g" << G4endl;

  VOL=RightLegBoneVol;RHO=RightLegBoneDensity;
  return physRightLegBone;
}
