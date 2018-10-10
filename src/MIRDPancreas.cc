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
#include "MIRDPancreas.hh"

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
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "HumanPhantomColour.hh"

#include "MirdMapConstants.hh"

MIRDPancreas::MIRDPancreas()
{
}

MIRDPancreas::~MIRDPancreas()
{

}

G4VPhysicalVolume* MIRDPancreas::Construct(const G4String& volumeName,G4VPhysicalVolume* mother,
						  const G4String& colourName
						  ,G4bool wireFrame, G4bool sensitivity)
{

  G4cout << "Construct "<< volumeName << G4endl;
 
 HumanPhantomMaterial* material = new HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;

 MirdMapConstants* mmc = new MirdMapConstants();  
  mmc->initTrunk();
  std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  G4double ct = scaleZ*trk[ageGroup].Ct;//our zcenter is MIRD 0.5*ct z coordinate
  mmc->initPancreas();
  std::map<std::string,pancreas_struct> pcr = mmc->GetPancreasMap();  
  G4double a = scaleXY*pcr[ageGroup].a;
  G4double b = scaleXY*pcr[ageGroup].b;
  G4double c = scaleZ*pcr[ageGroup].c;
  G4double x0 = scaleXY*pcr[ageGroup].x0;  
  G4double x1 = scaleXY*pcr[ageGroup].x1;  
  G4double z0 = scaleZ*pcr[ageGroup].z0;    
  delete mmc;
  //must put x on z axis to cut, so rotate about Y to match mird, so we have left z up x inside y!!!
  G4double ax= c;//scaleXY*3.*cm; //c
  G4double by= b;//scaleXY*1.*cm;//b
  G4double cz= a;//scaleZ*15.*cm;//a
  G4double zcut1= 0;//-a;//-scaleZ*15. *cm;// -a
  G4double zcut2= a;//0.0 *cm; 
  //TAKE INTO ACCOUNT ONLY PART BETWEEN ZCUTS!!!!!!!!!!!!!!!!!!NOT REMOVED BETWEEN CUTS!!!!!!!!!!!!!!!
  G4Ellipsoid* pancreasFirst =  new G4Ellipsoid("PancreasFirst",ax, by, cz,
						zcut1, zcut2);

  G4double xx = 2.0*c;//scaleXY*6. * cm;// 2*c
  G4double yy = 2.0*b;//scaleXY*2. * cm;// 2*b
  G4double zz = 2.0*(a-(x1-x0));//scaleZ*12. * cm; // cz - x1 = 3 cm
  G4Box* subtrPancreas = new G4Box("SubtrPancreas",xx/2., yy/2., zz/2.);
  G4double placeBoxATz=a;//good!
  G4SubtractionSolid* pancreas = new G4SubtractionSolid("pancreas",
							pancreasFirst,
							subtrPancreas,
							0, 
							//G4ThreeVector(-scaleXY*3 * cm,0.0,-scaleZ*9.*cm));
							G4ThreeVector(c,0.0,placeBoxATz));
  //whole y plane/place box CENTER at this coordinate of first solid!!
  
  G4LogicalVolume* logicPancreas = new G4LogicalVolume(pancreas, soft,
						       "logical" + volumeName,
						       0, 0, 0);
  G4RotationMatrix* rotation = new G4RotationMatrix();
  //rotation ->rotateY(90. * degree);
  rotation ->rotateY(-90. * degree);//along +y clockwise!
  //NOTE GEANT4 COORDINATE SYSTEM = MIRD5 COORDINATE SYSTEM!!!!!+x at right,+y inside and +z up!!!!

  G4VPhysicalVolume* physPancreas = new G4PVPlacement(rotation,
						      //G4ThreeVector(-0. *cm, 0.0, scaleXY*2*cm),//x0, 0, 2 cm
							  G4ThreeVector(x0, 0.0, 0.5*(0.0-ct)+z0),////x0, 0, 2 cm
      			       "physicalPancreas",
  			       logicPancreas,
			       mother,
			       false,
			       0, checkOverlaps);
  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicPancreas->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  // G4VisAttributes* PancreasVisAtt = new G4VisAttributes(G4Colour(0.28,0.82,0.8));
  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* PancreasVisAtt = new G4VisAttributes(colour);
  PancreasVisAtt->SetForceSolid(wireFrame);
  logicPancreas->SetVisAttributes(PancreasVisAtt);

  G4cout << "Pancreas created !!!!!!" << G4endl;

  // Testing Pancreas Volume
  G4double PancreasVol = logicPancreas->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Pancreas = " << PancreasVol/cm3 << " cm^3" << G4endl;
  
  // Testing Pancreas Material
  G4String PancreasMat = logicPancreas->GetMaterial()->GetName();
  G4cout << "Material of Pancreas = " << PancreasMat << G4endl;
  
  // Testing Density
  G4double PancreasDensity = logicPancreas->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << PancreasDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double PancreasMass = (PancreasVol)*PancreasDensity;
  G4cout << "Mass of Pancreas = " << PancreasMass/gram << " g" << G4endl;

  VOL=PancreasVol;RHO=PancreasDensity;
  return physPancreas;
}
