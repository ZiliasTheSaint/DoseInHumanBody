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
#include "MIRDLowerLargeIntestine.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4EllipticalTube.hh"
#include "G4UnionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Torus.hh"
#include "HumanPhantomMaterial.hh"
#include "HumanPhantomColour.hh"

#include "MirdMapConstants.hh"

MIRDLowerLargeIntestine::MIRDLowerLargeIntestine()
{
}

MIRDLowerLargeIntestine::~MIRDLowerLargeIntestine()
{

}

G4VPhysicalVolume* MIRDLowerLargeIntestine::Construct(const G4String& volumeName,
							     G4VPhysicalVolume* mother,
							     const G4String& colourName, G4bool wireFrame,G4bool sensitivity)
{
  G4cout << "Construct "<< volumeName <<G4endl;
 
 HumanPhantomMaterial* material = new HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;
 
 MirdMapConstants* mmc = new MirdMapConstants();  
  mmc->initTrunk();
  std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  G4double ct = scaleZ*trk[ageGroup].Ct;//our zcenter is MIRD 0.5*ct z coordinate
  mmc->initUpperLargeIntestine();
  std::map<std::string,upperLargeIntestine_struct> uli = mmc->GetUpperLargeIntestineMap();  
  G4double y0 = scaleXY*uli[ageGroup].y0;
  mmc->initLowerLargeIntestine();
  std::map<std::string,lowerLargeIntestine_struct> lli = mmc->GetLowerLargeIntestineMap();  
  G4double a = scaleXY*lli[ageGroup].a;
  G4double b = scaleXY*lli[ageGroup].b;
  G4double d = scaleXY*lli[ageGroup].d;
  G4double x1 = scaleXY*lli[ageGroup].x1;  
  G4double mx = scaleXY*lli[ageGroup].mx;  
  G4double my = scaleXY*lli[ageGroup].my;  
  G4double z1 = scaleZ*lli[ageGroup].z1;    
  G4double z2 = scaleZ*lli[ageGroup].z2;
  mmc->initSigmoidColon();
  std::map<std::string,sigmoidColon_struct> sgc = mmc->GetSigmoidColonMap();  
  G4double aC = scaleXY*sgc[ageGroup].a;
  G4double bC = scaleXY*sgc[ageGroup].b;  
  G4double dC = scaleXY*sgc[ageGroup].d;
  G4double x0C = scaleXY*sgc[ageGroup].x0;  
  G4double z0C = scaleZ*sgc[ageGroup].z0;  
  G4double r1C = scaleZ*sgc[ageGroup].R1;     //the key!!! 
  G4double r2C = scaleZ*sgc[ageGroup].R2;      
  delete mmc;
  G4double zcenter = 0.5*(z1+z2);
 G4double dx = a;//scaleXY*1.88 * cm; //a
 G4double dy = b;//scaleXY*2.13 *cm; //b
 G4double dz = 0.5*(z2-z1);//scaleZ*7.64 *cm; //(z1-z2)/2

 G4EllipticalTube* DescendingColonLowerLargeIntestine = new G4EllipticalTube("DiscendingColon",dx, dy, dz);


  G4double rmin= 0.0 *cm;
  G4double rmax = aC;//a;//scaleXY*1.88 * cm;//a
  G4double rtor= r1C;//scaleXY*5.72*cm; //R1
  G4double startphi= 0. * degree;
  G4double deltaphi= 90. * degree;

  G4Torus* SigmoidColonUpLowerLargeIntestine = new G4Torus("SigmoidColonUpLowerLargeIntestine",
							    rmin, rmax,rtor,
							    startphi, deltaphi);

  rtor = r2C;//scaleXY*3. * cm;//R2
  G4VSolid* SigmoidColonDownLowerLargeIntestine = new G4Torus("SigmoidColonDownLowerLargeIntestine",
							      rmin, rmax,
							      rtor,startphi,deltaphi);

  G4RotationMatrix* relative_rm =  new G4RotationMatrix();
  relative_rm -> rotateY(180. * degree);
  relative_rm -> rotateZ(90. * degree);

  G4UnionSolid*  SigmoidColonLowerLargeIntestine = new G4UnionSolid( "SigmoidColonLowerLargeIntestine",
								      SigmoidColonUpLowerLargeIntestine,
								      SigmoidColonDownLowerLargeIntestine,
								      relative_rm,
								      //G4ThreeVector(0.0,scaleXY*8.72*cm,0.0));// R1 + R2
									  G4ThreeVector(0.0,r1C+r2C,0.0));  
 
  G4RotationMatrix* relative_rm_2 =  new G4RotationMatrix();
  relative_rm_2 -> rotateX(90. * degree);

  G4UnionSolid* LowerLargeIntestine = new G4UnionSolid( "LowerLargeIntestine",
						       DescendingColonLowerLargeIntestine,
							SigmoidColonLowerLargeIntestine,
							relative_rm_2,
							//G4ThreeVector(-scaleXY*5.72*cm,0.0*cm, -scaleXY*7.64*cm)
							G4ThreeVector(-r1C,0.0*cm, -0.5*(z2-z1))
							); // -rtor,0, -dz


  G4LogicalVolume* logicLowerLargeIntestine = new G4LogicalVolume( LowerLargeIntestine, soft,
								   "logical" + volumeName,
								   0, 0, 0);
  
  G4VPhysicalVolume* physLowerLargeIntestine = new G4PVPlacement(0,           // R1+ R2, -2.36 (y0), z0 
	  //modified angle =0 =>same y0 as upperLarge or transversalColon!!!!
								 //G4ThreeVector(scaleXY*8.72*cm, -scaleZ*2.36*cm,-scaleXY*18.64 *cm),
								 //G4ThreeVector(r1C+r2C, -0.5*my,0.5*(0.0-ct)+zcenter),//-scaleZ*18.64 *cm),
					G4ThreeVector(r1C+r2C, y0,0.5*(0.0-ct)+zcenter),//-scaleZ*18.64 *cm),
								 "physicalLowerLargeIntestine",
								 logicLowerLargeIntestine,
								 mother,
								 false,
								 0, checkOverlaps);
  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicLowerLargeIntestine->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  //G4VisAttributes* LowerLargeIntestineVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* LowerLargeIntestineVisAtt = new G4VisAttributes(colour);
  LowerLargeIntestineVisAtt->SetForceSolid(wireFrame);
  logicLowerLargeIntestine->SetVisAttributes(LowerLargeIntestineVisAtt);

  G4cout << "LowerLargeIntestine created !!!!!!" << G4endl;

  // Testing LowerLargeIntestine Volume
  G4double LowerLargeIntestineVol = logicLowerLargeIntestine->GetSolid()->GetCubicVolume();
  G4cout << "Volume of LowerLargeIntestine = " << LowerLargeIntestineVol/cm3 << " cm^3" << G4endl;
  
  // Testing LowerLargeIntestine Material
  G4String LowerLargeIntestineMat = logicLowerLargeIntestine->GetMaterial()->GetName();
  G4cout << "Material of LowerLargeIntestine = " << LowerLargeIntestineMat << G4endl;
  
  // Testing Density
  G4double LowerLargeIntestineDensity = logicLowerLargeIntestine->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << LowerLargeIntestineDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double LowerLargeIntestineMass = (LowerLargeIntestineVol)*LowerLargeIntestineDensity;
  G4cout << "Mass of LowerLargeIntestine = " << LowerLargeIntestineMass/gram << " g" << G4endl;

  VOL=LowerLargeIntestineVol;RHO=LowerLargeIntestineDensity;
  return physLowerLargeIntestine;
}
