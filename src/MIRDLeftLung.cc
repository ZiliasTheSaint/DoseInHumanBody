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
#include "MIRDLeftLung.hh"
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
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Box.hh"
#include "HumanPhantomColour.hh"

#include "MirdMapConstants.hh"

MIRDLeftLung::MIRDLeftLung()
{
}

MIRDLeftLung::~MIRDLeftLung()
{

}

G4VPhysicalVolume* MIRDLeftLung::Construct(const G4String& volumeName,G4VPhysicalVolume* mother, 
						  const G4String& colourName, G4bool wireFrame,G4bool sensitivity)
{

  G4cout << "Construct " << volumeName << G4endl;
 
 HumanPhantomMaterial* material = new HumanPhantomMaterial();
 G4Material* lung_material = material -> GetMaterial("lung_material");
 delete material;

 MirdMapConstants* mmc = new MirdMapConstants();  
  mmc->initTrunk();
  std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  G4double ct = scaleZ*trk[ageGroup].Ct;//our zcenter is MIRD 0.5*ct z coordinate
  mmc->initLungs();
  std::map<std::string,lungs_struct> lgn = mmc->GetLungsMap();  
  G4double a = scaleXY*lgn[ageGroup].a;
  G4double b = scaleXY*lgn[ageGroup].b;
  G4double c = scaleZ*lgn[ageGroup].c;
  G4double x0 = scaleXY*lgn[ageGroup].x0;  
  G4double z0 = scaleZ*lgn[ageGroup].z0;
  //G4double x1r = scaleXY*lgn[ageGroup].x1r;
  G4double x1l = scaleXY*lgn[ageGroup].x1l;
  //G4double y1r = scaleXY*lgn[ageGroup].y1r;=>y1l
  G4double y1l = scaleXY*lgn[ageGroup].y1l;
  //G4double z1r = scaleZ*lgn[ageGroup].z1r;=>z0
  //G4double z2r = scaleZ*lgn[ageGroup].z2r;=>z2l
  G4double z2l = scaleZ*lgn[ageGroup].z2l;
  delete mmc;

 G4double ax = a;//scaleXY*5. *cm; //a
 G4double by = b;//scaleXY*7.5 *cm; //b
 G4double cz = c;//scaleZ*24.*cm; //c
 G4double zcut1 = 0.0 *cm; 
 G4double zcut2=cz;//scaleZ*24. *cm;
 
 G4Ellipsoid* oneLung = new G4Ellipsoid("OneLung",ax, by, cz, zcut1,zcut2);
 
 ax= a-(-x0+abs(x1l));//see right lung comments
 by= b+y1l;
 cz= z2l-z0;
 //G4Ellipsoid* subtrLung = new G4Ellipsoid("subtrLung",ax, by, cz);
 G4Box* subtrLung = new G4Box("subtrLung",ax, by, cz);//must be a box!!!!!!!!!!!!!!!!!!!!!!
 G4double zshift=z0-z0;
 
 G4SubtractionSolid* lung1 =  new G4SubtractionSolid("Lung1", oneLung,subtrLung,
		0, G4ThreeVector(-a,-b,zshift));//-a here

 //THESE SUBTRACTIONS ARE NOT GOOD..see http://h-pylori-symptoms.blogspot.ro/2012/09/healthy-lungs-tips-surely-will-help-you_21.html
 /*
 ax= scaleXY*5.*cm; 
 by= scaleXY*7.5*cm; 
 cz= scaleZ*24.*cm;
 G4Ellipsoid* subtrLung = new G4Ellipsoid("subtrLung",ax, by, cz);
 // y<0
 G4double dx = scaleXY*5.5* cm;
 G4double dy = scaleXY*8.5 * cm;
 G4double dz = scaleZ*24. * cm;
 G4Box* box = new G4Box("Box", dx, dy, dz);
 //G4SubtractionSolid* section = new G4SubtractionSolid("BoxSub", subtrLung, box, 0, G4ThreeVector(0.*cm, 8.5* cm, 0.*cm)); 
G4SubtractionSolid* section2 = new G4SubtractionSolid("BoxSub2", subtrLung, box, 0, G4ThreeVector(0.*cm, scaleXY*8.5* cm, 0.*cm)); 
//G4SubtractionSolid* lung1 =  new G4SubtractionSolid("Lung1", oneLung,
//				       section,
//				       0, G4ThreeVector(6.*cm,0*cm,0.0*cm)); 
 G4SubtractionSolid* lung2 =  new G4SubtractionSolid("Lung2", oneLung,
					       section2,
					      0, G4ThreeVector(-scaleXY*6.*cm,0*cm,0.0*cm));
 // G4RotationMatrix* matrix = new G4RotationMatrix();  
 // matrix->rotateX(180. * degree);
 //matrix ->rotateZ(180.*degree);
 //matrix -> rotateY(180.* degree);

 //G4UnionSolid* lungs = new G4UnionSolid("Lungs", lung1, lung2, matrix, G4ThreeVector(17*cm, 0., 0.));

 */
 G4LogicalVolume* logicLeftLung = new G4LogicalVolume(lung1//oneLung//lung2
	 ,lung_material,
						  "logical" + volumeName, 0, 0, 0); 
  

  G4VPhysicalVolume* physLeftLung = new G4PVPlacement(0,
	  //G4ThreeVector(scaleXY*8.50 *cm, 0.0*cm, scaleZ*8.5*cm),
	  G4ThreeVector(x0, 0.0, 0.5*(0.0-ct)+z0),
						  "physicalLeftLung",                    
  			       logicLeftLung,
			       mother,
			       false,
			       0, checkOverlaps);


  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicLeftLung->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
 }

  // Visualization Attributes
  //G4VisAttributes* LeftLungVisAtt = new G4VisAttributes(G4Colour(0.25,0.41,0.88));
  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* LeftLungVisAtt = new G4VisAttributes(colour);
  LeftLungVisAtt->SetForceSolid(wireFrame);
  logicLeftLung->SetVisAttributes(LeftLungVisAtt); 

  G4cout << "LeftLung created !!!!!!" << G4endl;

  // Testing LeftLung Volume
  G4double LeftLungVol = logicLeftLung->GetSolid()->GetCubicVolume();
 
 G4cout << "Volume of LeftLung = " << (LeftLungVol)/cm3 << " cm^3" << G4endl;
  
  // Testing LeftLung Material
  G4String LeftLungMat = logicLeftLung->GetMaterial()->GetName();
  G4cout << "Material of LeftLung = " << LeftLungMat << G4endl;
  
  // Testing Density
  G4double LeftLungDensity = logicLeftLung->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << LeftLungDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double LeftLungMass = (LeftLungVol)*LeftLungDensity;
  G4cout << "Mass of LeftLung = " << LeftLungMass/gram << " g" << G4endl;
  
  VOL=LeftLungVol;RHO=LeftLungDensity;
  return physLeftLung;
}
