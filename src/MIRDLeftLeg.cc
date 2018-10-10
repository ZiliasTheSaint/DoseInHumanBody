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
#include "MIRDLeftLeg.hh"

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

MIRDLeftLeg::MIRDLeftLeg()
{
}

MIRDLeftLeg::~MIRDLeftLeg()
{
}

G4VPhysicalVolume* MIRDLeftLeg::Construct(const G4String& volumeName, G4VPhysicalVolume* mother, 
						 const G4String& colourName, G4bool wireFrame,G4bool sensitivity)
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
  //R2 leg=>2R2 (its diameter) for leg

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
  
  G4LogicalVolume* logicLeftLeg = new G4LogicalVolume(leg1,
						   soft,
						   "logical" + volumeName,
						    0, 0, 0);
						   
  G4RotationMatrix* rm = new G4RotationMatrix();
  rm -> rotateX(90.* degree);

  G4VPhysicalVolume* physLeftLeg = new G4PVPlacement(rm,
                                
				//G4ThreeVector(scaleXY*10. * cm, -scaleZ*40. * cm, 0. *cm),
				G4ThreeVector(at, -0.5*dz, 0.0),
      			       "physicalLeftLeg",
  			       logicLeftLeg,
			       mother,
			       false,
			       0, checkOverlaps);

  // Sensitive Body Part
  if (sensitivity == true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicLeftLeg->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  
  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* LeftLegVisAtt = new G4VisAttributes(colour);
  LeftLegVisAtt->SetForceSolid(wireFrame);
  logicLeftLeg->SetVisAttributes(LeftLegVisAtt);

  G4cout << "LeftLeg created !!!!!!" << G4endl;

  // Testing LeftLeg Volume
  G4double LeftLegVol = logicLeftLeg->GetSolid()->GetCubicVolume();
  G4cout << "Volume of LeftLeg = " << LeftLegVol/cm3 << " cm^3" << G4endl;
  
  // Testing LeftLeg Material
  G4String LeftLegMat = logicLeftLeg->GetMaterial()->GetName();
  G4cout << "Material of LeftLeg = " << LeftLegMat << G4endl;
  
  // Testing Density
  G4double LeftLegDensity = logicLeftLeg->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << LeftLegDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double LeftLegMass = (LeftLegVol)*LeftLegDensity;
  G4cout << "Mass of LeftLeg = " << LeftLegMass/gram << " g" << G4endl;

  VOL=LeftLegVol;RHO=LeftLegDensity;
  return physLeftLeg;
}
