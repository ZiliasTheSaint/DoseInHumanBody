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
//D. Fulea, National Institute of Public Health Bucharest, Cluj-Napoca Regional Center, Romania
//

#include "MIRDRightAdrenal.hh"

#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "HumanPhantomMaterial.hh"
#include "G4SDManager.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4Ellipsoid.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include "G4EllipticalTube.hh"
#include "G4Box.hh"
#include "G4UnionSolid.hh"
#include "HumanPhantomColour.hh"

#include "MirdMapConstants.hh"

MIRDRightAdrenal::MIRDRightAdrenal()
{
}

MIRDRightAdrenal::~MIRDRightAdrenal()
{
}

G4VPhysicalVolume* MIRDRightAdrenal::Construct(const G4String& volumeName,G4VPhysicalVolume* mother,
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
  mmc->initAdrenals();
  std::map<std::string,adrenals_struct> adr = mmc->GetAdrenalsMap();  
  G4double a = scaleXY*adr[ageGroup].a;
  G4double b = scaleXY*adr[ageGroup].b;
  G4double c = scaleZ*adr[ageGroup].c;
  G4double x0 = scaleXY*adr[ageGroup].x0;
  G4double y0 = scaleXY*adr[ageGroup].y0;
  G4double z0 = scaleZ*adr[ageGroup].z0;
  G4double theta = adr[ageGroup].theta;
  delete mmc;

 G4double ax= a;//scaleXY*1.5 *cm; //a
 G4double by= b;//scaleXY*0.5 *cm; //b
 G4double cz= c;//scaleZ*5.0 *cm; //c
 
 G4VSolid* rightAdrenal = new G4Ellipsoid("OneRightAdrenal",ax, by, cz, 0. *cm, cz); 
 
 
 G4LogicalVolume* logicRightAdrenal = new G4LogicalVolume(rightAdrenal,
						     soft,
						     "logical" + volumeName,
						     0, 0, 0);

 //============================
  G4RotationMatrix* rm = new G4RotationMatrix();
  rm -> rotateZ(theta);
  //===================
  G4VPhysicalVolume* physRightAdrenal = new G4PVPlacement(rm//0 
	  ,G4ThreeVector(-x0,//-scaleXY*4.5*cm,  // xo
							y0,//	     scaleXY*6.5 *cm, //yo
								0.5*(0.0-ct)+z0),//    scaleZ*3. *cm),//z>z0!=>centerOfHalfEllips go down to MIRD 0 coord and then up by z0!
  			       "physicalRightAdrenal", logicRightAdrenal,
			       mother,
			       false,
			       0, checkOverlaps);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicRightAdrenal->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  //  G4VisAttributes* RightAdrenalVisAtt = new G4VisAttributes(G4Colour(0.72,0.52,0.04));
  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* RightAdrenalVisAtt = new G4VisAttributes(colour);
  RightAdrenalVisAtt->SetForceSolid(wireFrame);
  logicRightAdrenal->SetVisAttributes(RightAdrenalVisAtt);

  G4cout << "Right RightAdrenal created !!!!!!" << G4endl;

  // Testing RightAdrenal Volume
  G4double RightAdrenalVol = logicRightAdrenal->GetSolid()->GetCubicVolume();
  G4cout << "Volume of RightAdrenal = " << RightAdrenalVol/cm3 << " cm^3" << G4endl;
  
  // Testing RightAdrenal Material
  G4String RightAdrenalMat = logicRightAdrenal->GetMaterial()->GetName();
  G4cout << "Material of RightAdrenal = " << RightAdrenalMat << G4endl;
  
  // Testing Density
  G4double RightAdrenalDensity = logicRightAdrenal->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << RightAdrenalDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double RightAdrenalMass = (RightAdrenalVol)*RightAdrenalDensity;
  G4cout << "Mass of RightAdrenal = " << RightAdrenalMass/gram << " g" << G4endl;

  VOL=RightAdrenalVol;RHO=RightAdrenalDensity;
  return physRightAdrenal;
}
