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
#include "MIRDHeart.hh"
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
#include "G4Sphere.hh"
#include "G4UnionSolid.hh"

#include "MirdMapConstants.hh"

MIRDHeart::MIRDHeart()
{
	VOL=740*cm3;
	G4double d = 0.9869 *g/cm3;//softTissue
	RHO=d;
}

MIRDHeart::~MIRDHeart()
{

}

G4VPhysicalVolume* MIRDHeart::Construct(const G4String&,G4VPhysicalVolume*,
				    const G4String&, G4bool, G4bool)
{
  
 /*G4cout << " MIRD Heart is not available yet !!!! " << G4endl;
  
 G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;

 G4double ax= 4.00* cm;
 G4double by= 4.00 *cm;
 G4double cz= 7.00 *cm;
 G4double zcut1= -7.00 *cm;
 G4double zcut2= 0.0 *cm;

 G4Ellipsoid* heart1 =  new G4Ellipsoid("Heart1",ax, by, cz, zcut1, zcut2);

 G4double rmin =0.*cm;
 G4double rmax = 3.99*cm;
 G4double startphi = 0. * degree;
 G4double deltaphi = 360. * degree;
 G4double starttheta = 0. * degree;
 G4double deltatheta = 90. * degree;
 
 G4Sphere* heart2 = new G4Sphere("Heart2", rmin,rmax,
                                           startphi,   deltaphi,
                                           starttheta, deltatheta);

 G4UnionSolid* heart = new G4UnionSolid("Heart", heart1, heart2);

 G4LogicalVolume* logicHeart = new G4LogicalVolume(heart, soft,
						   "HeartVolume",
						   0, 0, 0);

  G4RotationMatrix* matrix = new G4RotationMatrix();
  matrix -> rotateY(25. * degree); 
  
  // Define rotation and position here!
  G4VPhysicalVolume* physHeart = new G4PVPlacement(matrix,G4ThreeVector(0.0,-3.0*cm, 15.32 *cm),
      			       "physicalHeart",
  			       logicHeart,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicHeart->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* HeartVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  HeartVisAtt->SetForceSolid(true);
  logicHeart->SetVisAttributes(HeartVisAtt);

  G4cout << "Heart created !!!!!!" << G4endl;

  // Testing Heart Volume
  G4double HeartVol = logicHeart->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Heart = " << HeartVol/cm3 << " cm^3" << G4endl;
  
  // Testing Heart Material
  G4String HeartMat = logicHeart->GetMaterial()->GetName();
  G4cout << "Material of Heart = " << HeartMat << G4endl;
  
  // Testing Density
  G4double HeartDensity = logicHeart->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << HeartDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double HeartMass = (HeartVol)*HeartDensity;
  G4cout << "Mass of Heart = " << HeartMass/gram << " g" << G4endl;

  VOL=HeartVol;RHO=HeartDensity;
  return physHeart;
  */
	//VOL = VOL*scaleXY*scaleXY*scaleZ;//scalexy^2*scaleZ=M/M0=V/V0
	MirdMapConstants* mmc = new MirdMapConstants();
  mmc->initHeart();
  std::map<std::string,heart_struct> esp = mmc->GetHeartMap();    
  VOL = scaleXY*scaleXY*scaleZ*esp[ageGroup].volume;  
  return 0;
}
