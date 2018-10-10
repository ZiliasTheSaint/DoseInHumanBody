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
#include "MIRDLeftBreast.hh"
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
#include "G4EllipticalTube.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "HumanPhantomColour.hh"

#include "MirdMapConstants.hh"

MIRDLeftBreast::MIRDLeftBreast()
{
}

MIRDLeftBreast::~MIRDLeftBreast()
{
}

G4VPhysicalVolume* MIRDLeftBreast::Construct(const G4String& volumeName,G4VPhysicalVolume* mother,
						    const G4String& colourName, G4bool wireFrame, G4bool sensitivity)
{
  G4cout << "Construct " << volumeName << G4endl;
 
 HumanPhantomMaterial* material = new HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;

 MirdMapConstants* mmc = new MirdMapConstants();
 mmc->initTrunk();
 std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
 G4double at = scaleXY*trk[ageGroup].At;
 G4double bt = scaleXY*trk[ageGroup].Bt;
 G4double ct = scaleXY*trk[ageGroup].Ct;
 mmc->initBreasts();
 std::map<std::string,breasts_struct> brs = mmc->GetBreastsMap();  
 G4double a = scaleXY*brs[ageGroup].a;
 G4double b = scaleXY*brs[ageGroup].b;
 G4double c = scaleXY*brs[ageGroup].c;
 G4double x0 = scaleXY*brs[ageGroup].x0;
 G4double z0 = scaleZ*brs[ageGroup].z0;
 delete mmc;
 G4double y0 = -bt*sqrt(1.0-x0*x0/at/at);//~-8.66cm
 //G4cout << "Construct " << volumeName << G4endl;
 G4double ax= a;//scaleXY*4.95* cm;
 G4double by= b;//scaleXY*4.35* cm;
 G4double cz= c;//scaleZ*4.15*cm;

 G4Ellipsoid* oneLeftBreast = new G4Ellipsoid("OneLeftBreast",
				      ax, by, cz);

 G4double dx= at;//scaleXY*20.* cm;
 G4double dy= bt;//scaleXY*10.* cm;
 G4double dz= 0.5*ct;//scaleZ*35.*cm;

 G4EllipticalTube* Trunk = new G4EllipticalTube("Trunk",dx, dy, dz );

 //Note: breats is on mother volume...unrottated. Our real Trunk is rottated! THIS TRUNK IS NOT!!!!!!!!!!
 G4RotationMatrix* rm_relative = new G4RotationMatrix();
 rm_relative -> rotateX(90. * degree);

 G4SubtractionSolid* breast = new G4SubtractionSolid("LeftBreast",
						      oneLeftBreast,
						      Trunk,
						      0,//rm_relative,//ROTATE LATER!!!!!!!!!!!!!!!!!!!!!!!!!!
                              //G4ThreeVector(-scaleXY*10.*cm, 0.0*cm, -scaleXY*8.66*cm));							  
							  //G4ThreeVector(-scaleXY*10.*cm, 0.0*cm, -scaleZ*8.66*cm));

							  //G4ThreeVector(-x0, 0.0, y0));
							  G4ThreeVector(-x0, -y0, 0.0));
 //substract 1-2=>breast-trunk=>trunk center goes to these breast coord and subtract
 //-x0 means left breast -y0 means in front!
 //NOTE THAT GEANT4 COORDINATE SYSTEM=MIRD5 COORDINATE SYSTEM!!! 
 //(We know: xbreast=xtrunk;ybreast=-ztrunk;zbreast=ytrunk <=> CONCLUSION: xMIRD=x;yMIRD=-z;zMIRD=y;)
 
 G4LogicalVolume* logicLeftBreast = new G4LogicalVolume(breast, soft,
							 "logical" + volumeName, 0, 0,0);
     
  // Define rotation and position here!
  G4VPhysicalVolume* physLeftBreast = 
	  //new G4PVPlacement(0,G4ThreeVector(scaleXY*10*cm, scaleZ*52.* cm, scaleXY*8.66 *cm),
	  
	  ////our phantom versus original MIRD phantom:xMIRD=x;yMIRD=-z;zMIRD=y and therefore:
	  ////////////////////////////////////x=xMIRD;y=zMIRD;z=-yMIRD
	  
	  //best way: rotate here along with the trunk!!!!!
	  new G4PVPlacement(rm_relative//0
	     ,G4ThreeVector(x0, z0, -y0),
      			       "physicalLeftBreast",
  			       logicLeftBreast,
			       mother,
			       false,
			       0, checkOverlaps);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicLeftBreast->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  
  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* LeftBreastVisAtt = new G4VisAttributes(colour);
  LeftBreastVisAtt->SetForceSolid(wireFrame);
 
  logicLeftBreast->SetVisAttributes(LeftBreastVisAtt);

  G4cout << "LeftBreast created !!!!!!" << G4endl;
  
  // Testing LeftBreast Volume
  G4double LeftBreastVol = logicLeftBreast->GetSolid()->GetCubicVolume();
  G4cout << "Volume of LeftBreast = " << LeftBreastVol/cm3 << " cm^3" << G4endl;
  
  // Testing LeftBreast Material
  G4String LeftBreastMat = logicLeftBreast->GetMaterial()->GetName();
  G4cout << "Material of LeftBreast = " << LeftBreastMat << G4endl;
  
  // Testing Density
  G4double LeftBreastDensity = logicLeftBreast->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << LeftBreastDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double LeftBreastMass = (LeftBreastVol)*LeftBreastDensity;
  G4cout << "Mass of LeftBreast = " << LeftBreastMass/gram << " g" << G4endl;

  VOL=LeftBreastVol;RHO=LeftBreastDensity;
  return physLeftBreast;
}
