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
#include "MIRDSkull.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "HumanPhantomMaterial.hh"
#include "G4EllipticalTube.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Ellipsoid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Box.hh"
#include "G4UnionSolid.hh"
#include "G4VSolid.hh"
#include "HumanPhantomColour.hh"
#include "MirdMapConstants.hh"

MIRDSkull::MIRDSkull()
{
}

MIRDSkull::~MIRDSkull()
{

}

G4VPhysicalVolume* MIRDSkull::Construct(const G4String& volumeName,G4VPhysicalVolume* mother,
					       const G4String& colourName,
					       G4bool wireFrame,G4bool sensitivity)
{
  
  HumanPhantomMaterial* material = new HumanPhantomMaterial();
   
  G4cout << "Construct "<<volumeName <<G4endl;
   
  G4Material* skeleton = material -> GetMaterial("skeleton");
 
  delete material;

  MirdMapConstants* mmc = new MirdMapConstants();
  mmc->initTrunk();
  std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  G4double bt = scaleXY*trk[ageGroup].Bt;
  mmc->initHead();
  std::map<std::string,head_struct> hed = mmc->GetHeadMap();    
  G4double rh = scaleXY*hed[ageGroup].Rh;
  G4double ch0 = scaleZ*hed[ageGroup].Ch0;
  G4double ch1 = scaleZ*hed[ageGroup].Ch1;
  G4double ch2 = scaleZ*hed[ageGroup].Ch2;
  mmc->initBrain();
  std::map<std::string,brain_struct> brn = mmc->GetBrainMap();    
  G4double a = scaleXY*brn[ageGroup].a;
  G4double b = scaleXY*brn[ageGroup].b;
  G4double c = scaleZ*brn[ageGroup].c;
  //G4double delta = scaleZ*brn[ageGroup].delta;
  mmc->initSkull();
  std::map<std::string,skull_struct> skl = mmc->GetSkullMap();    
  G4double d = scaleXY*skl[ageGroup].d;//fat=>thicker
  delete mmc;

  G4double deltaRh2=bt-rh;//move neck in the back instead....TO WORK!!!
  // Outer cranium
  G4double ax =a+d;//scaleXY* 6.8 * cm;//a out skull
  G4double by =b+d;//scaleXY* 9.8 * cm; // bout
  G4double cz =c+d;//scaleZ* 8.3 * cm; //cout
 
  G4Ellipsoid* craniumOut =  new G4Ellipsoid("CraniumOut", ax, by, cz);

  ax = a;//+0.1*d;//scaleXY*6. * cm; //a in
  by = b;//+0.1*d;//scaleXY*9. * cm; //b in 
  cz= c;//+0.1*d;//scaleZ*6.5 * cm; // cin
 
  G4Ellipsoid* craniumIn =  new G4Ellipsoid("CraniumIn", ax, by, cz);
 

  G4SubtractionSolid* cranium =  new G4SubtractionSolid("Cranium",
						      craniumOut,
							craniumIn,0,
						      //G4ThreeVector(0.0, 0.0,scaleZ*1. * cm));
							  G4ThreeVector(0.0, 0.0, 0.0));

  G4LogicalVolume* logicSkull = new G4LogicalVolume(cranium, skeleton, 
						    "logical" + volumeName,
						    0, 0, 0);
  
  // Define rotation and position here!
  G4VPhysicalVolume* physSkull = new G4PVPlacement(0,
						   //G4ThreeVector(0., 0.,scaleZ*7.75 * cm),
						  
						   //G4ThreeVector(0.0, 0.0,0.5*(2.0*c-ch0)+ch1),
						   G4ThreeVector(0.0, -deltaRh2,0.5*(0.0-ch0)+ch1+ch0),//They forget the neck!!!!!!!
						   "physicalSkull",
						   logicSkull,
						   mother,
						   false,
						   0, checkOverlaps);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicSkull->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  //G4VisAttributes* SkullVisAtt = new G4VisAttributes(G4Colour(0.46,0.53,0.6));
  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* SkullVisAtt = new G4VisAttributes(colour);
  SkullVisAtt->SetForceSolid(wireFrame); 
  SkullVisAtt->SetLineWidth(4.* mm);
  logicSkull->SetVisAttributes(SkullVisAtt);

  G4cout << "Skull created !!!!!!" << G4endl;


  // Testing Skull Volume
  G4double SkullVol = logicSkull->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Skull = " << SkullVol/cm3 << " cm^3" << G4endl;
  
  // Testing Skull Material
  G4String SkullMat = logicSkull->GetMaterial()->GetName();
  G4cout << "Material of Skull = " << SkullMat << G4endl;
  
  // Testing Density
  G4double SkullDensity = logicSkull->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << SkullDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double SkullMass = (SkullVol)*SkullDensity;
  G4cout << "Mass of Skull = " << SkullMass/gram << " g" << G4endl;

  VOL=SkullVol;RHO=SkullDensity;
  return physSkull;
}
