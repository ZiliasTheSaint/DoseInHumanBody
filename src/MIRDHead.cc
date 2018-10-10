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
#include "MIRDHead.hh"
#include "globals.hh"
#include "HumanPhantomMaterial.hh"
#include "G4SDManager.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4UnionSolid.hh"
#include "G4Ellipsoid.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include "G4EllipticalTube.hh"
#include "HumanPhantomColour.hh"

#include "G4Tubs.hh"
#include "MirdMapConstants.hh"

MIRDHead::MIRDHead()
{  
  material = new HumanPhantomMaterial();
}

MIRDHead::~MIRDHead()
{
  delete material;
}

G4VPhysicalVolume* MIRDHead::Construct(const G4String& volumeName,G4VPhysicalVolume* mother,
					      const G4String& colourName, G4bool wireFrame, G4bool sensitivity)
{
  G4cout << "Construct " << volumeName <<G4endl;

  //G4cout << "SCALE XY = " << scaleXY << "SCALE Z="<<scaleZ <<G4endl;//WORKS!!

  G4Material* soft = material -> GetMaterial("soft_tissue");  
  
  MirdMapConstants* mmc = new MirdMapConstants();
  mmc->initTrunk();
  std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  G4double bt = scaleXY*trk[ageGroup].Bt;
  G4double ct = scaleZ*trk[ageGroup].Ct;
  mmc->initSpine();
  std::map<std::string,spine_struct> spn = mmc->GetSpineMap();  
  //G4double a = scaleXY*spn[ageGroup].a;
  G4double b = scaleXY*spn[ageGroup].b;//max
  G4double y0 = scaleXY*spn[ageGroup].y0;//max
  mmc->initHead();
  std::map<std::string,head_struct> hed = mmc->GetHeadMap();  
  G4double rh = scaleXY*hed[ageGroup].Rh;
  G4double ch0 = scaleZ*hed[ageGroup].Ch0;
  G4double ah = scaleXY*hed[ageGroup].Ah;
  G4double bh = scaleXY*hed[ageGroup].Bh;
  G4double ch1 = scaleZ*hed[ageGroup].Ch1;
  G4double ch2 = scaleZ*hed[ageGroup].Ch2;
  delete mmc;
  //ensure it encompass the upper spine even if phantom neck does not look well!:
  G4double deltaRh2=bt-rh;//move neck in the back instead....TO WORK!!!
  //G4double deltaRh=y0+2.0*b-rh;
  //if (deltaRh>0)
//	rh=rh+deltaRh;//2.0*b;
  G4Tubs* neck = new G4Tubs("neck", 0.0, rh, 0.5*ch0, 0.0, 2.0*pi);//right circular cylinder
  G4Ellipsoid* head10 = new G4Ellipsoid("Head10", ah, bh, ch2, 0.0, ch2);//top half ellipsoid..zcut1=0, zcut2=zsemiaxis
  G4EllipticalTube* head20 = new G4EllipticalTube("Head20", ah, bh, 0.5*ch1);//right elliptical tube
  G4UnionSolid* headWithoutNeck = new G4UnionSolid("HeadWN",head20,head10,
  				0, // Rotation 
  				G4ThreeVector(0.0, 0.0, 0.5*ch1)//relative to first solid here, which is head20...at his top
				//ellipsoid is already a half!!!General 0.5z1size+0.5z2size. Centers are translated
				);
  
  G4UnionSolid* headWithNeck = new G4UnionSolid("HeadN",neck,headWithoutNeck
  				,0, // Rotation 
  				//G4ThreeVector(0.0, 0.0, 0.5*(ch0+ch1))//relative to first solid here, which is neck...at his top
				G4ThreeVector(0.0, -deltaRh2, 0.5*(ch0+ch1))//relative to first solid here, which is neck...at his top
				//hwithoutn 0.5zsize = 0.5(ch1+ch2)
				);
  // Ellipsoid
 /* G4double ax = scaleXY*7.0 * cm;
  G4double by = scaleXY*10.0 * cm;
  G4double cz = scaleZ*8.50 * cm;
  G4double zcut1 = 0.0 * cm;
  G4double zcut2 = scaleZ*8.5 * cm;

  G4Ellipsoid* head1 = new G4Ellipsoid("Head1", ax, by, cz, zcut1, zcut2);

   G4double dx = scaleXY*7.0 * cm;
   G4double dy = scaleXY*10.0 * cm;
   G4double dz = scaleZ*7.75 * cm;
 

  G4EllipticalTube* head2 = new G4EllipticalTube("Head2", dx, dy, dz);

  G4UnionSolid* head = new G4UnionSolid("Head",head2,head1,
  				0, // Rotation 
  				G4ThreeVector(0.* cm, 0.*cm, scaleZ*7.7500 * cm) );*/

  G4LogicalVolume* logicHead = //new G4LogicalVolume(head, soft,"logical" + volumeName,
	  //new G4LogicalVolume(neck, soft,"logical" + volumeName,
	  new G4LogicalVolume(headWithNeck, soft,"logical" + volumeName,
						   0, 0,0);
  G4RotationMatrix* rm = new G4RotationMatrix();
  rm -> rotateX(90.* degree);
  
  // Define rotation and position here! TAKEN FROM TUTORIALS:
  G4VPhysicalVolume* physHead = new G4PVPlacement(rm,//ROTATION OF MOTHER FRAME!!!!!!!!!!!!!!!!!!!!!!!
						  //G4ThreeVector(0.* cm,scaleZ*77.75 *cm, 0.*cm),
						  //G4ThreeVector(0.0,ct+0.5*ch0, 0.0),
						  G4ThreeVector(0.0,ct+0.5*ch0, -deltaRh2),//POSITION IN ROTATED FRAME!!!!!!!!!!!!!
						  "physicalHead",
						  logicHead,
						  mother,
						  false,
						  0, checkOverlaps);
  //placement with G4Transform3D(rm,G4ThreeVector(...)) instead of above means ROT OF DAUGHTER FRAME and POSITION IN MOTHER FRAME!!

  // Sensitive Body Part
 
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicHead->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
    G4cout <<SDman->FindSensitiveDetector("BodyPartSD")->GetName()<< G4endl;
    SDman->FindSensitiveDetector("BodyPartSD")->SetVerboseLevel(1);

  }

  // Visualization Attributes

  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* HeadVisAtt = new G4VisAttributes(colour);

 HeadVisAtt->SetForceSolid(wireFrame);
 HeadVisAtt->SetVisibility(true);
 // HeadVisAtt->SetLineWidth(0.7* mm);
  //HeadVisAtt-> SetForceAuxEdgeVisible(true);
  logicHead->SetVisAttributes(HeadVisAtt);

  // Testing Head Volume
  G4double HeadVol = logicHead->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Head = " << HeadVol/cm3 << " cm^3" << G4endl;
  //G4cout << "Volume of Head = " << HeadVol << " ???? " << G4endl;
  
  // Testing Head Material
  G4String HeadMat = logicHead->GetMaterial()->GetName();
  G4cout << "Material of Head = " << HeadMat << G4endl;
  
  // Testing Density
  G4double HeadDensity = logicHead->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << HeadDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double HeadMass = (HeadVol)*HeadDensity;
  G4cout << "Mass of Head = " << HeadMass/gram << " g" << G4endl;

  //HeadMass = logicHead->GetMass();
  //G4cout << "Mass of Head function = " << HeadMass/gram << " g" << G4endl;
  
  VOL=HeadVol;RHO=HeadDensity;
  return physHead;
}
