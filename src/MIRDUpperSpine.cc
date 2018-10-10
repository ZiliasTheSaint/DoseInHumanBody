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
#include "MIRDUpperSpine.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4VisAttributes.hh"
#include "HumanPhantomMaterial.hh"
#include "G4EllipticalTube.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "HumanPhantomColour.hh"

#include "MirdMapConstants.hh"

MIRDUpperSpine::MIRDUpperSpine()
{
}

MIRDUpperSpine::~MIRDUpperSpine()
{
 
}

G4VPhysicalVolume* MIRDUpperSpine::Construct(const G4String& volumeName,
						    G4VPhysicalVolume* mother, 
						    const G4String& colourName
						    , G4bool wireFrame,G4bool sensitivity)
{
  HumanPhantomMaterial* material = new HumanPhantomMaterial();
   
  G4cout << "Construct " <<volumeName <<G4endl;
   
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
  G4double z4 = scaleZ*spn[ageGroup].z4;
  G4double spineHeight=z4-z3;//10.54 adult
  delete mmc;

  G4double deltaRh2=bt-rh;//move neck in the back instead....TO WORK!!!

  G4EllipticalTube* upperSpine = new G4EllipticalTube("UpperSpine",a, b, 0.5*spineHeight);
  /*
  G4double dx = scaleXY*2. *cm;
  G4double dy = scaleXY*2.5 *cm;
  G4double dz = scaleZ*4.25*cm;

 G4EllipticalTube* upperSpine = new G4EllipticalTube("UpperSpine",dx, dy, dz);

 G4double xx = scaleXY*20. * cm;
 G4double yy = scaleXY*10. * cm;
 G4double zz = scaleZ*5. * cm;

 G4Box* subtraction = new G4Box("box", xx/2., yy/2., zz/2.);

 G4RotationMatrix* matrix = new G4RotationMatrix();
 matrix -> rotateX(-25.* deg); 

 G4SubtractionSolid* upper_spine = new G4SubtractionSolid("upperspine",upperSpine, subtraction,
							  matrix, G4ThreeVector(0., -scaleXY*2.5 * cm, scaleZ*5.5* cm));
*/
 G4LogicalVolume* logicUpperSpine = //new G4LogicalVolume(upper_spine, skeleton, 
	 new G4LogicalVolume(upperSpine, skeleton, 
							"logical" + volumeName,
							0, 0, 0);  
  // Define rotation and position here!
  G4VPhysicalVolume* physUpperSpine = new G4PVPlacement(0,
			       // G4ThreeVector(0.0, scaleXY*5.5 *cm, -scaleZ*3.5 *cm),
				    //G4ThreeVector(0.0, y0, 0.0),//-0.5*headHeight),
					G4ThreeVector(0.0, y0-deltaRh2,0.5*(spineHeight-ch0)+0.0),
      			       "physicalUpperSpine",
  			       logicUpperSpine,
			       mother,
			       false,
			       0, checkOverlaps);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicUpperSpine->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  //G4VisAttributes* UpperSpineVisAtt = new G4VisAttributes(G4Colour(0.46,0.53,0.6));
  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* UpperSpineVisAtt = new G4VisAttributes(colour);
 
  UpperSpineVisAtt->SetForceSolid(wireFrame);
  logicUpperSpine->SetVisAttributes(UpperSpineVisAtt);

  G4cout << "UpperSpine created !!!!!!" << G4endl;
 
  // Testing UpperSpine Volume
  G4double UpperSpineVol = logicUpperSpine->GetSolid()->GetCubicVolume();
  G4cout << "Volume of UpperSpine = " << UpperSpineVol/cm3 << " cm^3" << G4endl;
  
  // Testing UpperSpine Material
  G4String UpperSpineMat = logicUpperSpine->GetMaterial()->GetName();
  G4cout << "Material of UpperSpine = " << UpperSpineMat << G4endl;
  
  // Testing Density
  G4double UpperSpineDensity = logicUpperSpine->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << UpperSpineDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double UpperSpineMass = (UpperSpineVol)*UpperSpineDensity;
  G4cout << "Mass of UpperSpine = " << UpperSpineMass/gram << " g" << G4endl;

  VOL=UpperSpineVol;RHO=UpperSpineDensity; 
  return physUpperSpine;
}
