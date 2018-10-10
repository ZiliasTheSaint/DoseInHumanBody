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
#include "MIRDMiddleLowerSpine.hh"
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
#include "HumanPhantomColour.hh"

#include "MirdMapConstants.hh"

MIRDMiddleLowerSpine::MIRDMiddleLowerSpine()
{
}

MIRDMiddleLowerSpine::~MIRDMiddleLowerSpine()
{
}

G4VPhysicalVolume* MIRDMiddleLowerSpine::Construct(const G4String& volumeName,
							  G4VPhysicalVolume* mother,
							  const G4String& colourName
							  , G4bool wireFrame, G4bool sensitivity)
{
 HumanPhantomMaterial* material = new HumanPhantomMaterial();
   
 G4cout << "Construct "<< volumeName  <<G4endl;   
  G4Material* skeleton = material -> GetMaterial("skeleton"); 
  delete material;

  MirdMapConstants* mmc = new MirdMapConstants();
  mmc->initTrunk();
  std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  G4double bt = scaleXY*trk[ageGroup].Bt;
  G4double ct = scaleZ*trk[ageGroup].Ct;//
  mmc->initHead();
  std::map<std::string,head_struct> hed = mmc->GetHeadMap();  
  G4double rh = scaleXY*hed[ageGroup].Rh;
  G4double ch0 = scaleZ*hed[ageGroup].Ch0;
  mmc->initSpine();
  std::map<std::string,spine_struct> spn = mmc->GetSpineMap();  
  G4double a = scaleXY*spn[ageGroup].a;
  G4double b = scaleXY*spn[ageGroup].b;
  G4double y0 = scaleXY*spn[ageGroup].y0;
  G4double z3 = scaleZ*spn[ageGroup].z3;
  G4double z1 = scaleZ*spn[ageGroup].z1;
  G4double spineHeight=z3-z1;//70-22 = 48 adult
  G4double spineCenter=0.5*(z3+z1);
  delete mmc;
 
  //G4double deltaRh2=bt-rh;

  G4double dx = a;//scaleXY*2. *cm;
  G4double dy = b;//scaleXY*2.5 *cm;
  G4double dz = 0.5*spineHeight;//scaleZ*24. *cm;

  G4VSolid* middleLowerSpine = new G4EllipticalTube("MiddleLowerSpine",dx, dy, dz);

  G4LogicalVolume* logicMiddleLowerSpine = new G4LogicalVolume( middleLowerSpine, skeleton,
								"logical" + volumeName,
								0, 0, 0);   
  // Define rotation and position here!
  G4VPhysicalVolume* physMiddleLowerSpine = new G4PVPlacement(0,
	  //G4ThreeVector(0.0 *cm, scaleXY*5.5 * cm,scaleZ*11. * cm),
	  G4ThreeVector(0.0 *cm, y0,0.5*(0.0-ct)+spineCenter),//scaleZ*11. * cm),
							      "physicalMiddleLowerSpine",
							      logicMiddleLowerSpine,
							      mother,
							      false,
							      0, checkOverlaps);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicMiddleLowerSpine->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  // G4VisAttributes* MiddleLowerSpineVisAtt = new G4VisAttributes(G4Colour(0.46,0.53,0.6));
 
  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* MiddleLowerSpineVisAtt = new G4VisAttributes(colour);
  MiddleLowerSpineVisAtt->SetForceSolid(wireFrame);
  logicMiddleLowerSpine->SetVisAttributes(MiddleLowerSpineVisAtt);

  G4cout << "MiddleLowerSpine created !!!!!!" << G4endl;

  // Testing MiddleLowerSpine Volume
  G4double MiddleLowerSpineVol = logicMiddleLowerSpine->GetSolid()->GetCubicVolume();
  G4cout << "Volume of MiddleLowerSpine = " << MiddleLowerSpineVol/cm3 << " cm^3" << G4endl;
  
  // Testing MiddleLowerSpine Material
  G4String MiddleLowerSpineMat = logicMiddleLowerSpine->GetMaterial()->GetName();
  G4cout << "Material of MiddleLowerSpine = " << MiddleLowerSpineMat << G4endl;
  
  // Testing Density
  G4double MiddleLowerSpineDensity = logicMiddleLowerSpine->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << MiddleLowerSpineDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double MiddleLowerSpineMass = (MiddleLowerSpineVol)*MiddleLowerSpineDensity;
  G4cout << "Mass of MiddleLowerSpine = " << MiddleLowerSpineMass/gram << " g" << G4endl;

  VOL=MiddleLowerSpineVol;RHO=MiddleLowerSpineDensity;
  return physMiddleLowerSpine;
}
