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
#include "MIRDSpleen.hh"

#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "HumanPhantomMaterial.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Ellipsoid.hh"
#include "HumanPhantomColour.hh"

#include "MirdMapConstants.hh"

MIRDSpleen::MIRDSpleen()
{
}

MIRDSpleen::~MIRDSpleen()
{

}

G4VPhysicalVolume* MIRDSpleen::Construct(const G4String& volumeName,G4VPhysicalVolume* mother,
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
  mmc->initSpleen();
  std::map<std::string,spleen_struct> spl = mmc->GetSpleenMap();  
  G4double a = scaleXY*spl[ageGroup].a;
  G4double b = scaleXY*spl[ageGroup].b;
  G4double c = scaleZ*spl[ageGroup].c;
  G4double x0 = scaleXY*spl[ageGroup].x0;  
  G4double y0 = scaleXY*spl[ageGroup].y0;  
  G4double z0 = scaleZ*spl[ageGroup].z0;    
  delete mmc;

 G4double ax= a;//scaleXY*3.5 *cm;
 G4double by= b;//scaleXY*2. *cm;
 G4double cz= c;//scaleZ*6. * cm; 

 G4Ellipsoid* spleen = new G4Ellipsoid("spleen", ax, by, cz);


  G4LogicalVolume* logicSpleen = new G4LogicalVolume(spleen, soft,
						     "logical" + volumeName,
						      0, 0, 0);
  
  // Define rotation and position here!
  G4VPhysicalVolume* physSpleen = new G4PVPlacement(0,
						    //G4ThreeVector(scaleXY*11. *cm, scaleXY*3. *cm, scaleZ*2.*cm), // ztrans = half trunk lenght - z0
				G4ThreeVector(x0, y0, 0.5*(0.0-ct)+z0),//scaleZ*2.*cm), // ztrans = half trunk lenght - z0
      			       "physicalSpleen",
  			       logicSpleen,
			       mother,
			       false,
			       0, checkOverlaps);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicSpleen->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  // G4VisAttributes* SpleenVisAtt = new G4VisAttributes(G4Colour(0.41,0.41,0.41));
  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* SpleenVisAtt = new G4VisAttributes(colour);
  SpleenVisAtt->SetForceSolid(wireFrame);
  logicSpleen->SetVisAttributes(SpleenVisAtt);

  G4cout << "Spleen created !!!!!!" << G4endl;

  // Testing Spleen Volume
  G4double SpleenVol = logicSpleen->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Spleen = " << SpleenVol/cm3 << " cm^3" << G4endl;
  
  // Testing Spleen Material
  G4String SpleenMat = logicSpleen->GetMaterial()->GetName();
  G4cout << "Material of Spleen = " << SpleenMat << G4endl;
  
  // Testing Density
  G4double SpleenDensity = logicSpleen->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << SpleenDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double SpleenMass = (SpleenVol)*SpleenDensity;
  G4cout << "Mass of Spleen = " << SpleenMass/gram << " g" << G4endl;


  VOL=SpleenVol;RHO=SpleenDensity;
  return physSpleen;
}
