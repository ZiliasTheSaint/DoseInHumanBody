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
#include "MIRDTrunk.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "HumanPhantomMaterial.hh"
#include "G4EllipticalTube.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "HumanPhantomColour.hh"

#include "MirdMapConstants.hh"

MIRDTrunk::MIRDTrunk()
{
}

MIRDTrunk::~MIRDTrunk()
{

}

G4VPhysicalVolume* MIRDTrunk::Construct(const G4String& volumeName, G4VPhysicalVolume* mother, 
					      const G4String& colourName, G4bool wireFrame,G4bool sensitivity)
{

  HumanPhantomMaterial* material = new HumanPhantomMaterial();
   
  G4cout << "Construct " << volumeName << G4endl;
   
  G4Material* soft = material -> GetMaterial("soft_tissue");
 
  delete material;
  //=================================================
  MirdMapConstants* mmc = new MirdMapConstants();
  mmc->initTrunk();
  std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();
  G4double at = trk[ageGroup].At;
  G4double bt = trk[ageGroup].Bt;
  G4double ct = trk[ageGroup].Ct;
  //G4cout << "AAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTT = " << at  <<"; age group= "<<ageGroup<<G4endl;//WORKS!!
  

  // MIRD Male trunk
  ct=0.5*ct;//half
  G4double dx = scaleXY*at;//20. * cm;
  G4double dy = scaleXY*bt;//10. * cm;
  G4double dz = scaleZ*ct;//35. * cm;

  G4EllipticalTube* trunk = new G4EllipticalTube("Trunk",dx, dy, dz);

  G4LogicalVolume* logicTrunk = new G4LogicalVolume(trunk, soft, 
	 					    "logical" + volumeName,
						    0, 0, 0);
  G4RotationMatrix* rm = new G4RotationMatrix();
  rm -> rotateX(90.* degree);

  // Define rotation and position here!
  G4VPhysicalVolume* physTrunk = new G4PVPlacement(rm,
				 //G4ThreeVector(0.* cm, scaleZ*35.0 *cm, 0.*cm),
				 G4ThreeVector(0.* cm, scaleZ*ct, 0.*cm),
      			       "physicalTrunk",
  			       logicTrunk,
			       mother,
			       false,
			       0, checkOverlaps);

  // Sensitive Body Part
  if (sensitivity == true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicTrunk->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  delete mmc;

  // Visualization Attributes
  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* TrunkVisAtt = new G4VisAttributes(colour);
  TrunkVisAtt->SetForceSolid(wireFrame);
  TrunkVisAtt->SetVisibility(true);//false);//true);
  logicTrunk->SetVisAttributes(TrunkVisAtt);

  G4cout << "Trunk created !!!!!!" << G4endl;

  // Testing Trunk Volume
  G4double TrunkVol = logicTrunk->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Trunk = " << TrunkVol/cm3 << " cm^3" << G4endl;
  
  // Testing Trunk Material
  G4String TrunkMat = logicTrunk->GetMaterial()->GetName();
  G4cout << "Material of Trunk = " << TrunkMat << G4endl;
  
  // Testing Density
  G4double TrunkDensity = logicTrunk->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << TrunkDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double TrunkMass = (TrunkVol)*TrunkDensity;
  G4cout << "Mass of Trunk = " << TrunkMass/gram << " g" << G4endl;

  VOL=TrunkVol;RHO=TrunkDensity;
  return physTrunk;
}
