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
#include "MIRDRightLung.hh"
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

MIRDRightLung::MIRDRightLung()
{
}

MIRDRightLung::~MIRDRightLung()
{

}

G4VPhysicalVolume* MIRDRightLung::Construct(const G4String& volumeName,G4VPhysicalVolume* mother, 
						   const G4String& colourName, G4bool wireFrame,G4bool sensitivity)
{

  G4cout << "Construct " << volumeName <<G4endl;
 
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
  G4double x1r = scaleXY*lgn[ageGroup].x1r;
  G4double y1r = scaleXY*lgn[ageGroup].y1r;
  G4double z1r = scaleZ*lgn[ageGroup].z1r;
  G4double z2r = scaleZ*lgn[ageGroup].z2r;
  delete mmc;

 G4double ax = a;//scaleXY*5. *cm; //a
 G4double by = b;//scaleXY*7.5 *cm; //b
 G4double cz = c;//scaleZ*24.*cm; //c
 G4double zcut1 = 0.0 *cm; 
 G4double zcut2=cz;//scaleZ*24. *cm;
  
 G4Ellipsoid* oneLung = new G4Ellipsoid("OneLung",ax, by, cz, zcut1,zcut2);//half

 ax= a-(x0-abs(x1r));//x>x1r exclusion when below conditions are met!!x1r<0,x0=8.50=>a-(x0-|x1r|) left gap!!
 //in fact: xmine=xreal+x0; xreal>x1; x1<0.=>xmine>x1+x0=x0-|x1|
 //for left: xmine=xreal-x0;x<x1l exclusion=>xmine<x1l-x0
 by= b+y1r;//a shift by b follows!...y<y1r!!!!
 cz= z2r-z1r;//z1r<z<z2r
 //G4Ellipsoid* subtrLung = new G4Ellipsoid("subtrLung",ax, by, cz);//must be a box!!!!!!!!!!!!!!!!!!!!!!
 //box size y larger then onelung..ok!
 G4Box* subtrLung = new G4Box("subtrLung",ax, by, cz);//must be a box!!!!!!!!!!!!!!!!!!!!!!
 G4double zshift=z1r-z0;
 
 G4SubtractionSolid* lung1 =  new G4SubtractionSolid("Lung1", oneLung,subtrLung,
		0, G4ThreeVector(a,-b,zshift));

 //THESE SUBTRACTIONS ARE NOT GOOD..see http://h-pylori-symptoms.blogspot.ro/2012/09/healthy-lungs-tips-surely-will-help-you_21.html
 /*ax= scaleXY*5.*cm; 
 by= scaleXY*7.5*cm; 
 cz= scaleZ*24.*cm;
 G4Ellipsoid* subtrLung = new G4Ellipsoid("subtrLung",ax, by, cz);
  // y<0
 G4double dx = scaleXY*5.5* cm;
 G4double dy = scaleXY*8.5 * cm;
 G4double dz = scaleZ*24. * cm;
 G4Box* box = new G4Box("Box", dx, dy, dz);
 
 G4SubtractionSolid* section = new G4SubtractionSolid("BoxSub", subtrLung, box, 0, G4ThreeVector(0.*cm, scaleXY*8.5* cm, 0.*cm)); 
 //G4SubtractionSolid* section2 = new G4SubtractionSolid("BoxSub2", subtrLung, box, 0, G4ThreeVector(0.*cm, -8.5* cm, 0.*cm)); 
  G4SubtractionSolid* lung1 =  new G4SubtractionSolid("Lung1", oneLung,
					       section,
  0, G4ThreeVector(scaleXY*6.*cm,0*cm,0.0*cm));
*/
 G4LogicalVolume* logicRightLung = new G4LogicalVolume(lung1//oneLung//lung1//
	 ,lung_material,"logical" + volumeName, 0, 0, 0); 
  
 G4RotationMatrix* matrix = new G4RotationMatrix();
 matrix -> rotateZ(-360. * deg);
  G4VPhysicalVolume* physRightLung = new G4PVPlacement(0,
	  //G4ThreeVector(-scaleXY*8.50 *cm, 0.0*cm, scaleZ*8.5*cm),
	  G4ThreeVector(-x0, 0.0*cm, 0.5*(0.0-ct)+z0),
						  "physicalRightLung",                    
  			       logicRightLung,
			       mother,
			       false,
			       0, checkOverlaps);


  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicRightLung->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
 }

  // Visualization Attributes
  // G4VisAttributes* RightLungVisAtt = new G4VisAttributes(G4Colour(0.25,0.41,0.88));
  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* RightLungVisAtt = new G4VisAttributes(colour);
   RightLungVisAtt->SetForceSolid(wireFrame);
  logicRightLung->SetVisAttributes(RightLungVisAtt); 

  G4cout << "RightLung created !!!!!!" << G4endl;

  // Testing RightLung Volume
  G4double RightLungVol = logicRightLung->GetSolid()->GetCubicVolume();
 
 G4cout << "Volume of RightLung = " << (RightLungVol)/cm3 << " cm^3" << G4endl;
  
  // Testing RightLung Material
  G4String RightLungMat = logicRightLung->GetMaterial()->GetName();
  G4cout << "Material of RightLung = " << RightLungMat << G4endl;
  
  // Testing Density
  G4double RightLungDensity = logicRightLung->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << RightLungDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double RightLungMass = (RightLungVol)*RightLungDensity;
  G4cout << "Mass of RightLung = " << RightLungMass/gram << " g" << G4endl;
  
  VOL=RightLungVol;RHO=RightLungDensity;
  return physRightLung;
}
