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

#include "MIRDRightKidney.hh"

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

MIRDRightKidney::MIRDRightKidney()
{
}

MIRDRightKidney::~MIRDRightKidney()
{
}

G4VPhysicalVolume* MIRDRightKidney::Construct(const G4String& volumeName,
						     G4VPhysicalVolume* mother, 
						     const G4String& colourName,
						     G4bool wireFrame,G4bool sensitivity)
{
  G4cout << "Construct " << volumeName << G4endl;
 
 HumanPhantomMaterial* material = new HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;
 
 MirdMapConstants* mmc = new MirdMapConstants();  
  mmc->initTrunk();
  std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  G4double ct = scaleZ*trk[ageGroup].Ct;//our zcenter is MIRD 0.5*ct z coordinate
  mmc->initKidneys();
  std::map<std::string,kidneys_struct> kdn = mmc->GetKidneysMap();  
  G4double a = scaleXY*kdn[ageGroup].a;
  G4double b = scaleXY*kdn[ageGroup].b;
  G4double c = scaleZ*kdn[ageGroup].c;
  G4double x0 = scaleXY*kdn[ageGroup].x0;
  G4double y0 = scaleXY*kdn[ageGroup].y0;
  G4double z0 = scaleZ*kdn[ageGroup].z0;
  G4double x1 = scaleXY*kdn[ageGroup].x1;
  delete mmc;

 G4double ax= a;//scaleXY*4.5 *cm; //a
 G4double by= b;//scaleXY*1.5 *cm; //b
 G4double cz= c;//scaleZ*5.5 *cm; //c
 
 G4VSolid* oneRightKidney = new G4Ellipsoid("OneRightKidney",ax, by, cz); 
 
 G4double xx = 2.0*x1;//scaleXY*6. * cm; 
 G4double yy = ct;//scaleXY*12.00*cm; ////infinite plane=>max coord
 G4double zz = ct;//scaleZ*12.00*cm;////infinite plane=>max coord
 G4VSolid* subtrRightKidney = new G4Box("SubtrRightKidney",xx/2., yy/2., zz/2.);
 
 G4SubtractionSolid* kidney = new G4SubtractionSolid("RightKidney",
						     oneRightKidney,
						     subtrRightKidney,
						     0, 
						     G4ThreeVector(//scaleXY*6. *cm, //// x1 =>shift with 2x1 but x1 semiaxes=>cut at -2x1+x1=-x1
							 2.0*x1,
								   0.0 *cm,
								   0.0 * cm));

  G4LogicalVolume* logicRightKidney = new G4LogicalVolume(kidney,
						     soft,
						     "logical" + volumeName,
						     0, 0, 0);

  G4VPhysicalVolume* physRightKidney = new G4PVPlacement(0 ,G4ThreeVector(-x0,//-scaleXY*6.*cm,  // xo
								     y0,//scaleXY*6. *cm, //yo
								     0.5*(0.0-ct)+z0),//-scaleZ*2.50 *cm),//zo
  			       "physicalRightKidney", logicRightKidney,
			       mother,
			       false,
			       0, checkOverlaps);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicRightKidney->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  //G4VisAttributes* RightKidneyVisAtt = new G4VisAttributes(G4Colour(0.72,0.52,0.04));
  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* RightKidneyVisAtt = new G4VisAttributes(colour);
  RightKidneyVisAtt->SetForceSolid(wireFrame);
  logicRightKidney->SetVisAttributes(RightKidneyVisAtt);

  G4cout << "RightKidney created !!!!!!" << G4endl;

  // Testing RightKidney Volume
  G4double RightKidneyVol = logicRightKidney->GetSolid()->GetCubicVolume();
  G4cout << "Volume of RightKidney = " << RightKidneyVol/cm3 << " cm^3" << G4endl;
  
  // Testing RightKidney Material
  G4String RightKidneyMat = logicRightKidney->GetMaterial()->GetName();
  G4cout << "Material of RightKidney = " << RightKidneyMat << G4endl;
  
  // Testing Density
  G4double RightKidneyDensity = logicRightKidney->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << RightKidneyDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double RightKidneyMass = (RightKidneyVol)*RightKidneyDensity;
  G4cout << "Mass of RightKidney = " << RightKidneyMass/gram << " g" << G4endl;

  VOL=RightKidneyVol;RHO=RightKidneyDensity;
  return physRightKidney;
}
