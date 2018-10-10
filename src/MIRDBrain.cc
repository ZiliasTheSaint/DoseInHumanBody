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
#include "MIRDBrain.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Ellipsoid.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "HumanPhantomMaterial.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "HumanPhantomColour.hh"
#include "MirdMapConstants.hh"

MIRDBrain::MIRDBrain()
{
	VOL=0.;RHO=0.;//initialization
}

MIRDBrain::~MIRDBrain()
{

}

G4VPhysicalVolume* MIRDBrain::Construct(const G4String& volumeName,G4VPhysicalVolume* mother,  
					       const G4String& colourName, G4bool wireFrame,G4bool sensitivity)
{
  G4cout << "Construct " << volumeName << G4endl;
 HumanPhantomMaterial* material = new HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
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
  delete mmc;

  G4double deltaRh2=bt-rh;//move neck in the back instead....TO WORK!!!

 G4double ax =a;//scaleXY* 6. * cm;
 G4double by= b;//scaleXY*9. * cm;
 G4double cz = c;//scaleZ*6.5 * cm;

 G4Ellipsoid* brain = new G4Ellipsoid("Brain", ax, by, cz);
 

  G4LogicalVolume* logicBrain =  new G4LogicalVolume(brain, soft, 
						     "logical" + volumeName,
						     0, 0, 0);
  
  // Define rotation and position here!
  G4VPhysicalVolume* physBrain = new G4PVPlacement(0,
						  // G4ThreeVector(0.*cm, 0.*cm,scaleZ* 8.75 * cm),

						   //G4ThreeVector(0.0, 0.0,0.5*(2.0*c-ch0)+ch1),
						   G4ThreeVector(0.0, -deltaRh2,0.5*(0.0-ch0)+ch1+ch0),////They forget the neck!!!!!!!
						   
						   "physicalBrain",
						   logicBrain,
						   mother,
						   false,
						   0, checkOverlaps);
  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicBrain->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }


  // Visualization Attributes
  // G4VisAttributes* BrainVisAtt = new G4VisAttributes(G4Colour(0.41,0.41,0.41));
  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  
  G4VisAttributes* BrainVisAtt = new G4VisAttributes(colour);
  BrainVisAtt->SetForceSolid(wireFrame);
  BrainVisAtt->SetLineWidth(0.7* mm);
  logicBrain->SetVisAttributes(BrainVisAtt);

  // Testing Brain Volume
  G4double BrainVol = logicBrain->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Brain = " << BrainVol/cm3 << " cm^3" << G4endl;
  
  // Testing Brain Material
  G4String BrainMat = logicBrain->GetMaterial()->GetName();
  G4cout << "Material of Brain = " << BrainMat << G4endl;
  
  // Testing Density
  G4double BrainDensity = logicBrain->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << BrainDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double BrainMass = (BrainVol)*BrainDensity;
  G4cout << "Mass of Brain = " << BrainMass/gram << " g" << G4endl;

  VOL=BrainVol;RHO=BrainDensity;
  return physBrain;
}
