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
#include "MIRDUterus.hh"
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

MIRDUterus::MIRDUterus()
{
}

MIRDUterus::~MIRDUterus()
{
}

G4VPhysicalVolume* MIRDUterus::Construct(const G4String& volumeName,G4VPhysicalVolume* mother,  
						const G4String& colourName, G4bool wireFrame, G4bool sensitivity)
{
 
  G4cout << "Construct " << volumeName << G4endl;
 
 HumanPhantomMaterial* material = new HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;

 MirdMapConstants* mmc = new MirdMapConstants();  
  mmc->initTrunk();
  std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  G4double ct = scaleZ*trk[ageGroup].Ct;//our zcenter is MIRD 0.5*ct z coordinate
  mmc->initUterus();
  std::map<std::string,uterus_struct> utr = mmc->GetUterusMap();  
  G4double a = scaleXY*utr[ageGroup].a;
  G4double b = scaleXY*utr[ageGroup].b;
  G4double c = scaleZ*utr[ageGroup].c;
  G4double y0 = scaleXY*utr[ageGroup].y0;  
  G4double z0 = scaleZ*utr[ageGroup].z0;  
  G4double y1 = scaleXY*utr[ageGroup].y1;  
  delete mmc;

  //must rotate to use zcuts!!!
  //remember. y0,y1 are negative!!!!

 G4double ax= a;//scaleXY*2.5*cm; //a
 G4double by= c;//scaleXY*1.5*cm; //c
 G4double cz= b;//scaleZ*5.*cm; //b

 G4double zcut1= -b;//-scaleZ*5.* cm; //-b
 G4double zcut2= -(y1-y0);//scaleZ*2.5*cm; //y1-y0 nop -(y1-y0) like from -b. OK volume checked!!
 //ok..see pancreas comments!
 G4Ellipsoid* uterus = new G4Ellipsoid("Uterus",
				       ax, by, cz,
				       zcut1, zcut2);

  G4LogicalVolume* logicUterus = new G4LogicalVolume(uterus,
						     soft,
						     "logical" + volumeName,
						     0, 0, 0);


  G4RotationMatrix* rm = new G4RotationMatrix();
  rm -> rotateX(90.* degree); 
  
  // Define rotation and position here!
  G4VPhysicalVolume* physUterus = new G4PVPlacement(rm,
				//G4ThreeVector(0. *cm, scaleZ*2*cm,-scaleXY*21 *cm),
				G4ThreeVector(0. *cm, y0,0.5*(0.0-ct)+z0),//-scaleXY*21 *cm),
						    "physicalUterus", //y0
  			       logicUterus,
			       mother,
			       false,
			       0, checkOverlaps);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicUterus->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  //  G4VisAttributes* UterusVisAtt = new G4VisAttributes(G4Colour(0.85,0.44,0.84));

  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* UterusVisAtt = new G4VisAttributes(colour);  
  UterusVisAtt->SetForceSolid(wireFrame);
  logicUterus->SetVisAttributes(UterusVisAtt);

  G4cout << "Uterus created !!!!!!" << G4endl;

  // Testing Uterus Volume
  G4double UterusVol = logicUterus->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Uterus = " << UterusVol/cm3 << " cm^3" << G4endl;
  
  // Testing Uterus Material
  G4String UterusMat = logicUterus->GetMaterial()->GetName();
  G4cout << "Material of Uterus = " << UterusMat << G4endl;
  
  // Testing Density
  G4double UterusDensity = logicUterus->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << UterusDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double UterusMass = (UterusVol)*UterusDensity;
  G4cout << "Mass of Uterus = " << UterusMass/gram << " g" << G4endl;
  
  VOL=UterusVol;RHO=UterusDensity;
  return physUterus;
}
