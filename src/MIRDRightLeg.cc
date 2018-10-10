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
#include "MIRDRightLeg.hh"

#include "globals.hh"

#include "G4SDManager.hh"
#include "G4Cons.hh"

#include "G4VisAttributes.hh"
#include "HumanPhantomMaterial.hh"
#include "G4EllipticalTube.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UnionSolid.hh"
#include "HumanPhantomColour.hh"

#include "MirdMapConstants.hh"

MIRDRightLeg::MIRDRightLeg()
{
}

MIRDRightLeg::~MIRDRightLeg()
{
}

G4VPhysicalVolume* MIRDRightLeg::Construct(const G4String& volumeName,G4VPhysicalVolume* mother, 
						  const G4String& colourName, G4bool wireFrame, G4bool sensitivity)
{
 
  G4cout << "Construct"<< volumeName << G4endl;

  HumanPhantomMaterial* material = new HumanPhantomMaterial();
  G4Material* soft = material -> GetMaterial("soft_tissue");
 
  MirdMapConstants* mmc = new MirdMapConstants();
  mmc->initTrunk();
  std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  G4double at = scaleXY*trk[ageGroup].At;
  at=0.5*at;//divide in two parts=>two legs
  mmc->initLegs();
  std::map<std::string,legs_struct> lgs = mmc->GetLegsMap(); 
  G4double cl=scaleZ*lgs[ageGroup].Cl;
  G4double cpl=scaleZ*lgs[ageGroup].Cpl;
  delete mmc;

  G4double R2=(at/4.0)*(cpl-cl)/cpl;//adult=1 to adjust match leg bone

  G4double rmin1 = 0.* cm;
  G4double rmin2 = 0.* cm;
  G4double dz= cl;//scaleZ*80.0 * cm; 
  G4double rmax1= 2.0*R2;//scaleXY*2.0 * cm;
  G4double rmax2= at;//scaleXY*10. * cm;
  G4double startphi= 0.* degree;
  G4double deltaphi= 360. * degree;

  G4Cons* leg1 = new G4Cons("Leg1",  
			   rmin1, rmax1, 
			   rmin2, rmax2, dz/2., 
			   startphi, deltaphi);
  
  //  G4UnionSolid* RightLeg = new G4UnionSolid("RightLeg", leg1, leg1,0,
  //                                    G4ThreeVector(20. *cm, 0., 0.));

  G4LogicalVolume* logicRightLeg = new G4LogicalVolume(leg1,
						   soft,
						   "logical" + volumeName,
						    0, 0, 0);
						   
  G4RotationMatrix* rm = new G4RotationMatrix();
  rm -> rotateX(90.* degree);

  G4VPhysicalVolume* physRightLeg = new G4PVPlacement(rm,
                                
				//G4ThreeVector(-scaleXY*10. * cm, -scaleZ*40. * cm, 0. *cm),
				G4ThreeVector(-at, -0.5*dz, 0.0),
      			       "physicalRightLeg",
  			       logicRightLeg,
			       mother,
			       false,
			       0, checkOverlaps);

  // Sensitive Body Part
  if (sensitivity == true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicRightLeg->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  //G4VisAttributes* RightLegVisAtt = new G4VisAttributes(G4Colour(0.94,0.5,0.5));
  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* RightLegVisAtt = new G4VisAttributes(colour);
  RightLegVisAtt->SetForceSolid(wireFrame);
  logicRightLeg->SetVisAttributes(RightLegVisAtt);

  G4cout << "RightLeg created !!!!!!" << G4endl;

  // Testing RightLeg Volume
  G4double RightLegVol = logicRightLeg->GetSolid()->GetCubicVolume();
  G4cout << "Volume of RightLeg = " << RightLegVol/cm3 << " cm^3" << G4endl;
  
  // Testing RightLeg Material
  G4String RightLegMat = logicRightLeg->GetMaterial()->GetName();
  G4cout << "Material of RightLeg = " << RightLegMat << G4endl;
  
  // Testing Density
  G4double RightLegDensity = logicRightLeg->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << RightLegDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double RightLegMass = (RightLegVol)*RightLegDensity;
  G4cout << "Mass of RightLeg = " << RightLegMass/gram << " g" << G4endl;

  VOL=RightLegVol;RHO=RightLegDensity;
  return physRightLeg;
}
