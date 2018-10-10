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
#include "MIRDUpperLargeIntestine.hh"

#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4EllipticalTube.hh"
#include "G4UnionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "HumanPhantomMaterial.hh"
#include "HumanPhantomColour.hh"

#include "MirdMapConstants.hh"

MIRDUpperLargeIntestine::MIRDUpperLargeIntestine()
{
}

MIRDUpperLargeIntestine::~MIRDUpperLargeIntestine()
{
}

G4VPhysicalVolume* MIRDUpperLargeIntestine::Construct(const G4String& volumeName,
							     G4VPhysicalVolume* mother,
							     const G4String& colourName
							     , G4bool wireFrame,G4bool sensitivity)
{
  G4cout << "Construct " << volumeName<< G4endl;
  
 HumanPhantomMaterial* material = new HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;

 MirdMapConstants* mmc = new MirdMapConstants();  
  mmc->initTrunk();
  std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  G4double ct = scaleZ*trk[ageGroup].Ct;//our zcenter is MIRD 0.5*ct z coordinate
  mmc->initUpperLargeIntestine();
  std::map<std::string,upperLargeIntestine_struct> uli = mmc->GetUpperLargeIntestineMap();  
  G4double a = scaleXY*uli[ageGroup].a;
  G4double b = scaleXY*uli[ageGroup].b;
  G4double d = scaleXY*uli[ageGroup].d;
  G4double x0 = scaleXY*uli[ageGroup].x0;  
  G4double y0 = scaleXY*uli[ageGroup].y0;  
  G4double z1 = scaleZ*uli[ageGroup].z1;    
  G4double z2 = scaleZ*uli[ageGroup].z2; 
  mmc->initTransverseColon();
  std::map<std::string,transverseColon_struct> trc = mmc->GetTransverseColonMap();  
  G4double bC = scaleXY*trc[ageGroup].b;
  G4double cC = scaleZ*trc[ageGroup].c;
  G4double dC = scaleXY*trc[ageGroup].d;
  G4double y0C = scaleXY*trc[ageGroup].y0;  
  G4double z0C = scaleZ*trc[ageGroup].z0;  
  G4double x1C = scaleXY*trc[ageGroup].x1;      
  delete mmc;
  G4double zcenter = 0.5*(z1+z2);
 G4double dx = a;//scaleXY*2.5 * cm; // aU
 G4double dy = b;//scaleXY*2.5* cm; //bU
 G4double dz = 0.5*(z2-z1);//scaleZ*4.775 * cm; //dzU

  G4VSolid* AscendingColonUpperLargeIntestine = new G4EllipticalTube("AscendingColon",dx, dy, dz);
 
  dx = bC;//scaleXY*2.5 * cm;//bt
  dy = cC;//scaleXY*1.5 *cm;//ct
  dz = x1C;//scaleZ*10.5* cm;//x1t

  G4VSolid* TraverseColonUpperLargeIntestine = new G4EllipticalTube("TraverseColon",dx, dy, dz);

  G4RotationMatrix* relative_rm =  new G4RotationMatrix();
  relative_rm -> rotateX(90. * degree);
  relative_rm -> rotateZ(0. * degree);
  relative_rm -> rotateY(90. * degree);//THAT's IT!!
  G4UnionSolid* upperLargeIntestine = new G4UnionSolid("UpperLargeIntestine",
						      AscendingColonUpperLargeIntestine,
						      TraverseColonUpperLargeIntestine,
						      relative_rm, 
						      //G4ThreeVector(scaleXY*8.50 *cm, 0.0,scaleZ*6.275 * cm)); //,0,dzU + ct transverse
							  G4ThreeVector(-x0, 0.0,0.5*(z2-z1)+cC)); //,0,dzU + ct transverse
  
  //The translation is applied first to move the origin of coordinates.
  //Then the rotation is used to rotate the coordinate system of the second solid to the coordinate system of the first. 
  G4LogicalVolume* logicUpperLargeIntestine = new G4LogicalVolume(upperLargeIntestine, soft,
								  "logical" + volumeName, 
								  0, 0, 0);
 
  G4VPhysicalVolume* physUpperLargeIntestine = new G4PVPlacement(0,
								 //G4ThreeVector(-scaleXY*8.50 * cm, -scaleXY*2.36 *cm,-scaleZ*15.775 *cm),
								 G4ThreeVector(x0, y0,0.5*(0.0-ct)+zcenter),//-scaleZ*15.775 *cm),
								 "physicalUpperLargeIntestine",                 //xo, yo, zo ascending colon
  			       logicUpperLargeIntestine,
			       mother,
			       false,
			       0, checkOverlaps);

  // Initialise a single volume, positioned in a frame which is rotated by
  // *pRot and traslated by tlate, relative to the coordinate system of the
  // mother volume pMotherLogical.

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicUpperLargeIntestine->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  //  G4VisAttributes* UpperLargeIntestineVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* UpperLargeIntestineVisAtt = new G4VisAttributes(colour);
  UpperLargeIntestineVisAtt->SetForceSolid(wireFrame);
  logicUpperLargeIntestine->SetVisAttributes(UpperLargeIntestineVisAtt);

  G4cout << "UpperLargeIntestine created !!!!!!" << G4endl;

  // Testing UpperLargeIntestine Volume
  G4double UpperLargeIntestineVol = logicUpperLargeIntestine->GetSolid()->GetCubicVolume();
  G4cout << "Volume of UpperLargeIntestine = " << UpperLargeIntestineVol/cm3 << " cm^3" << G4endl;
  
  // Testing UpperLargeIntestine Material
  G4String UpperLargeIntestineMat = logicUpperLargeIntestine->GetMaterial()->GetName();
  G4cout << "Material of UpperLargeIntestine = " << UpperLargeIntestineMat << G4endl;
  
  // Testing Density
  G4double UpperLargeIntestineDensity = logicUpperLargeIntestine->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << UpperLargeIntestineDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double UpperLargeIntestineMass = (UpperLargeIntestineVol)*UpperLargeIntestineDensity;
  G4cout << "Mass of UpperLargeIntestine = " << UpperLargeIntestineMass/gram << " g" << G4endl;

  VOL=UpperLargeIntestineVol;RHO=UpperLargeIntestineDensity;
  return physUpperLargeIntestine;
}
