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
#include "MIRDRibCage.hh"

#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "HumanPhantomMaterial.hh"
#include "G4SubtractionSolid.hh"
#include "G4EllipticalTube.hh"
#include "G4PVReplica.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "HumanPhantomColour.hh"

#include "MirdMapConstants.hh"

MIRDRibCage::MIRDRibCage()
{
}

MIRDRibCage::~MIRDRibCage()
{
  }

G4VPhysicalVolume* MIRDRibCage::Construct(const G4String& volumeName, G4VPhysicalVolume* mother,
						 const G4String& colourName, G4bool wireFrame, G4bool sensitivity)
{
  HumanPhantomMaterial* material = new HumanPhantomMaterial();   
  G4cout << "Construct "<< volumeName  <<G4endl;   
  G4Material* skeleton = material -> GetMaterial("skeleton");
  G4Material* soft = material -> GetMaterial("soft_tissue"); 
  delete material; 

  MirdMapConstants* mmc = new MirdMapConstants();  
  mmc->initTrunk();
  std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  G4double ct = scaleZ*trk[ageGroup].Ct;//our zcenter is MIRD 0.5*ct z coordinate
  mmc->initRibCage();
  std::map<std::string,ribCage_struct> rbc = mmc->GetRibCageMap();  
  G4double a = scaleXY*rbc[ageGroup].a;
  G4double b = scaleXY*rbc[ageGroup].b;
  G4double d = scaleXY*rbc[ageGroup].d;
  G4double z1 = scaleZ*rbc[ageGroup].z1;  
  G4double z2 = scaleZ*rbc[ageGroup].z2;  
  G4double c = scaleZ*rbc[ageGroup].c;    
  delete mmc;
  G4double zcenter = 0.5*(z1+z2);
  G4double dx= a;//scaleXY*17. *cm; // a2
  G4double dy= b;//scaleXY*9.8 * cm; //b2
  G4double thickness= z2-z1;//scaleZ*32.4 * cm; // z2/2 of cage!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  G4EllipticalTube* outCage = new G4EllipticalTube("outCage",dx, dy, thickness/2.);

  dx = a-d;//scaleXY*16.4 * cm; // a1
  dy = b-d;//scaleXY*9.2 * cm; // b1
  G4double dz = thickness;//scaleZ*34. *cm; // z2/2

  G4EllipticalTube* inCage = new G4EllipticalTube("inCage",dx, dy, dz/2.);

  G4SubtractionSolid* cage = new G4SubtractionSolid("Cage",
						     outCage,
						     inCage, 0, G4ThreeVector(0.*cm, 0.*cm, 0. * cm)); 

 
  G4LogicalVolume* logicRibCage = new G4LogicalVolume(cage, soft, "LogicalCage", 0, 0, 0);

  G4VPhysicalVolume* physRibCage = new G4PVPlacement(0,
	  //G4ThreeVector(0.0, 0.0, thickness/2. + scaleZ*0.1 * cm),
	  G4ThreeVector(0.0, 0.0, 0.5*(0.0-ct)+zcenter),//thickness/2.),
						      // with respect to the trunk!!!!!!!!!!!!!!!!!!!!!!!!!!!
						      "physicalRibCage",
						      logicRibCage,
						      mother,
						      false,
						      0, checkOverlaps);
	
  
  G4double xx = a;//scaleXY*17.*cm;
  G4double yy = b;//scaleXY*9.8*cm;
  G4double ribThickness = c;//scaleZ*1.4*cm;
  G4EllipticalTube* rib_out = new G4EllipticalTube("rib_out",xx, yy, ribThickness/2.);	
  
  xx = a-d;//scaleXY*16.5 *cm;
  yy = b-d;//scaleXY*9.3 * cm;
  G4double zz = ribThickness;//scaleZ*1.5 * cm;  
  G4EllipticalTube* rib_in = new G4EllipticalTube("rib_in",xx, yy, zz/2.);
  G4SubtractionSolid* rib = new G4SubtractionSolid("rib",rib_out, rib_in);

  G4LogicalVolume* logicRib= new G4LogicalVolume(rib, skeleton, "logical" + volumeName, 0, 0, 0);

  //------------
 // int nofpairs = floor(thickness/(2.0*ribThickness));//always 12!!
  //---------------
  G4double rcenter = -0.5*thickness+ribThickness*0.5;// with respect to the RIB CAGE!!!!!!!!!!!!!!!!!!!!!
  //ribs=equally spaced
  physRib1 = new G4PVPlacement(0,
	  //G4ThreeVector(0.0, 0.0, scaleZ*(- 32.2*cm/2. + 0.8 *cm)),
	  G4ThreeVector(0.0, 0.0, rcenter),//scaleZ*(- 32.2*cm/2. + 0.8 *cm)),
			// with respect to the RIB CAGE!!!!!!!!!!!!!!!!!!!!!
						      "physicalRib",
						      logicRib,
						      physRibCage,
						      false,
						      0, checkOverlaps);

   rcenter = rcenter +2.0*ribThickness;
   physRib2 = new G4PVPlacement(0,//G4ThreeVector(0.0, 0.0,  scaleZ*( - 32.2*cm/2. + 0.8 *cm + 2.8 *cm)),//+2.8cm!!!
	   G4ThreeVector(0.0, 0.0,  rcenter),//+2.8cm!!!
						      "physicalRib",
						      logicRib,
						      physRibCage,
						      false,
						      0, checkOverlaps);
   rcenter = rcenter +2.0*ribThickness;
   physRib3 = new G4PVPlacement(0,//G4ThreeVector(0.0, 0.0, (-thickness/2. + scaleZ*0.8 * cm + scaleZ*5.6 *cm)),
						  G4ThreeVector(0.0, 0.0,  rcenter),
						      "physicalRib",
						      logicRib,
						      physRibCage,
						      false,
						      0, checkOverlaps);
   rcenter = rcenter +2.0*ribThickness;
  physRib4 = new G4PVPlacement(0,//G4ThreeVector(0.0, 0.0, (-thickness/2. +scaleZ* 0.8 * cm + scaleZ*8.4 *cm)),
						       G4ThreeVector(0.0, 0.0,  rcenter),
						      "physicalRib",
						      logicRib,
						      physRibCage,
						      false,
						      0, checkOverlaps);
  rcenter = rcenter +2.0*ribThickness;
 physRib5 = new G4PVPlacement(0,//G4ThreeVector(0.0, 0.0, (-thickness/2. + scaleZ*0.8 * cm + scaleZ*11.2 *cm)),
						      G4ThreeVector(0.0, 0.0,  rcenter),
						      "physicalRib",
						      logicRib,
						      physRibCage,
						      false,
						      0, checkOverlaps);
 rcenter = rcenter +2.0*ribThickness;
 physRib6 = new G4PVPlacement(0,//G4ThreeVector(0.0, 0.0, (-thickness/2. +  scaleZ*0.8 * cm + scaleZ*14. *cm)),
						      G4ThreeVector(0.0, 0.0,  rcenter),
						      "physicalRib",
						      logicRib,
						      physRibCage,
						      false,
						      0, checkOverlaps);
 rcenter = rcenter +2.0*ribThickness;
 physRib7 = new G4PVPlacement(0,//G4ThreeVector(0.0, 0.0, (-thickness/2. + scaleZ*0.8 *cm + scaleZ*16.8 *cm)),
						      G4ThreeVector(0.0, 0.0,  rcenter),
						      "physicalRib",
						      logicRib,
						      physRibCage,
						      false,
						      0, checkOverlaps);
 rcenter = rcenter +2.0*ribThickness;
 physRib8 = new G4PVPlacement(0,//G4ThreeVector(0.0, 0.0, (-thickness/2. + scaleZ*0.8 *cm + scaleZ*19.6 *cm)),
						      G4ThreeVector(0.0, 0.0,  rcenter),
						      "physicalRib",
						      logicRib,
						      physRibCage,
						      false,
						      0, checkOverlaps);
 rcenter = rcenter +2.0*ribThickness;
 physRib9 = new G4PVPlacement(0,//G4ThreeVector(0.0, 0.0, (-thickness/2. +scaleZ* 0.8*cm + scaleZ*22.4 *cm)),
						      G4ThreeVector(0.0, 0.0,  rcenter),
						      "physicalRib",
						      logicRib,
						      physRibCage,
						      false,
						      0, checkOverlaps);
 rcenter = rcenter +2.0*ribThickness;
 physRib10 = new G4PVPlacement(0,//G4ThreeVector(0.0, 0.0, (-thickness/2. +scaleZ* 0.8*cm + scaleZ*25.2 *cm)),
						     G4ThreeVector(0.0, 0.0,  rcenter),
						      "physicalRib",
						      logicRib,
						      physRibCage,
						      false,
						      0, checkOverlaps);
 rcenter = rcenter +2.0*ribThickness;
 physRib11 = new G4PVPlacement(0,//G4ThreeVector(0.0, 0.0, (-thickness/2. + scaleZ*0.8*cm + scaleZ*28. *cm)),
						     G4ThreeVector(0.0, 0.0,  rcenter),
						      "physicalRib",
						      logicRib,
						      physRibCage,
						      false,
						      0, checkOverlaps);
 rcenter = rcenter +2.0*ribThickness;
  physRib12 = new G4PVPlacement(0,//G4ThreeVector(0.0, 0.0, (-thickness/2. + scaleZ*0.8*cm + scaleZ*30.8 *cm)),
						      G4ThreeVector(0.0, 0.0,  rcenter),
						      "physicalRib",
						      logicRib,
						      physRibCage,
						      false,
						      0, checkOverlaps);
  
  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    //    logicRibCage->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
    logicRib->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  //logicRibCage -> SetVisAttributes(G4VisAttributes::Invisible);//THIS CAUSES ERROR!!!

  //G4VisAttributes* RibCageVisAtt = new G4VisAttributes(G4Colour(0.46,0.53,0.6));

  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* RibCageVisAtt = new G4VisAttributes(colour); 
  RibCageVisAtt->SetForceSolid(wireFrame);
  logicRib->SetVisAttributes(RibCageVisAtt);

  G4cout << "RibCage created !!!!!!" << G4endl;
  // Testing Pelvis Volume
  G4double RibCageVol = logicRib->GetSolid()->GetCubicVolume();
  G4cout << "Volume of RibCage = " << ((RibCageVol)*12.)/cm3 << " cm^3" << G4endl;
  
  // Testing RibCage Material
  G4String RibCageMat = logicRib->GetMaterial()->GetName();
  G4cout << "Material of RibCage = " << RibCageMat << G4endl;
  
  // Testing Density
  G4double RibCageDensity = logicRib->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << RibCageDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double RibCageMass = (RibCageVol)* RibCageDensity * 12;// 12 is the total number of ribs;
  G4cout << "Mass of RibCage = " << (RibCageMass)/gram << " g" << G4endl;

  
  VOL=RibCageVol*12.0;RHO=RibCageDensity;
  return physRibCage;
}
