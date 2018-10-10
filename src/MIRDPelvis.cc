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
#include "MIRDPelvis.hh"

#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "HumanPhantomMaterial.hh"
#include "G4EllipticalTube.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4Box.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "HumanPhantomColour.hh"

#include "MirdMapConstants.hh"

MIRDPelvis::MIRDPelvis()
{
}

MIRDPelvis::~MIRDPelvis()
{

}

G4VPhysicalVolume* MIRDPelvis::Construct(const G4String& volumeName,G4VPhysicalVolume* mother, 
						const G4String& colourName, G4bool wireFrame,G4bool sensitivity)
{
   HumanPhantomMaterial* material = new HumanPhantomMaterial();   
   G4cout << "Construct " << volumeName <<G4endl;   
  G4Material* skeleton = material -> GetMaterial("skeleton"); 
  delete material;
  
   MirdMapConstants* mmc = new MirdMapConstants();  
  mmc->initTrunk();
  std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  G4double ct = scaleZ*trk[ageGroup].Ct;//our zcenter is MIRD 0.5*ct z coordinate
  mmc->initPelvis();
  std::map<std::string,pelvis_struct> plv = mmc->GetPelvisMap();  
  G4double a1 = scaleXY*plv[ageGroup].a1;
  G4double b1 = scaleXY*plv[ageGroup].b1;
  G4double a2 = scaleXY*plv[ageGroup].a2;
  G4double b2 = scaleXY*plv[ageGroup].b2;  
  G4double y01 = scaleXY*plv[ageGroup].y01;  
  G4double y02 = scaleXY*plv[ageGroup].y02;  
  G4double y1 = scaleXY*plv[ageGroup].y1;  
  G4double z1 = scaleZ*plv[ageGroup].z1;    
  G4double z2 = scaleZ*plv[ageGroup].z2;    
  delete mmc;
  //0<z<z2
  double zcenter = 0.5*z2;
  G4double dx= a2;//scaleXY*12. *cm; // a2
  G4double dy= b2;//scaleXY*12. * cm; //b2
  G4double dz= 0.5*z2;//scaleZ*11. * cm; // z2/2
  G4VSolid* outPelvis = new G4EllipticalTube("OutPelvis",dx, dy, dz);
  
  dx = a1;//scaleXY*11.3 * cm; // a1
  dy = b1;//scaleXY*11.3* cm; // b1
  dz = 0.5*z2;//scaleZ*12.0 *cm; // z2/2!!!!!!!!!!!!!+1?????????? 
  G4VSolid* inPelvis = new G4EllipticalTube("InPelvis",dx, dy, dz);

  G4double placeInRelativeToOut_y=y01-y02;//-0.8; see first and second pelvis equation for tubes!!
  G4SubtractionSolid* firstPelvis = new G4SubtractionSolid("FirstPelvis", outPelvis, inPelvis, 0,
	  //G4ThreeVector(0.*cm, -scaleXY*0.8 *cm, 0. * cm)); 
	  G4ThreeVector(0.*cm, placeInRelativeToOut_y, 0. * cm)); 

  G4double x = 2.0*a2;//scaleXY*28. * cm; // a2 * 2//2.0*a2;//
  G4double y = 2.0*b2;//scaleXY*28. * cm; //b2*2//2.0*b2;//
  G4double z = z2;//scaleZ*24. *cm; // z2//z2;//
  G4VSolid* subPelvis = new G4Box("SubtrPelvis", x/2., y/2., z/2.);	   
    
  G4SubtractionSolid* secondPelvis = new G4SubtractionSolid("SecondPelvis",
							    firstPelvis,
							    subPelvis, 0, 
							    //G4ThreeVector(0.0,-scaleXY*14. * cm, 0.*cm));
								G4ThreeVector(0.0,-b2, 0.*cm));
                             // half of the y size of the box OK!!!!! y>y02 means y>0 
  /////////////////////////////////ok///////////////////////////////////////////////////////////////
  x = 2.0*a2;//scaleXY*28. * cm; 
  y = 2.0*(b2-(y1-y02));//scaleXY*28. * cm;
  z = 2.0*(z1);//scaleZ*24. *cm;==========z0 which is zmin=0!!!!
  subPelvis = new G4Box("SubtrPelvis", x/2., y/2., z/2.);	   
  G4double placeY=b2;
  G4double placeZ=-0.5*z2;
  //z<z1=>y>y1 removal
  //G4double placeRelativeTo_z = -(0.5*z2-(z1-zcenter));//(z2-z1);//enter with z1 inside!!//(0.5*z2-(z1-zcenter));//z<z1=>this+z2..ok!!
  //G4double placeRelativeTo_y = +(y-(y1-y02));//y1);//enter with y1 inside!!//(y1-y02));//abs(y02)));//y>y1 removal
  //note y1 has nothing to do with y0!!!nop...ca la lung si ca la pancreas!!!!!!
  //yes, but here y>y02 is already taken into account!!!!!!!!!!!!!!!!!!nop it has nothing to do wiyh it!
  G4SubtractionSolid* pelvis = new G4SubtractionSolid("Pelvis", secondPelvis, subPelvis,
						      0, 
						      //G4ThreeVector(0.0, scaleXY*22. * cm, -scaleZ*9. *cm)); 
							  //G4ThreeVector(0.0, placeRelativeTo_y, placeRelativeTo_z));

							  G4ThreeVector(0.0, placeY, placeZ)); 

 
  G4LogicalVolume* logicPelvis = new G4LogicalVolume(pelvis, skeleton,
						     "logical" + volumeName, 0, 0, 0);
  
 
  G4VPhysicalVolume* physPelvis = new G4PVPlacement
	  //(0,G4ThreeVector(0.0, -scaleXY*3. * cm,-scaleZ*24. * cm),// 0, y02, z position
	  (0,G4ThreeVector(0.0, y02,0.5*(0.0-ct)+zcenter),//-scaleZ*24. * cm),// 0, y02, z position
						    // with respect to the trunk 
      			       "physicalPelvis",
  			       logicPelvis,
			       mother,
			       false,
			       0, checkOverlaps);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicPelvis->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  //  G4VisAttributes* PelvisVisAtt = new G4VisAttributes(G4Colour(0.46,0.53,0.6));
 
  HumanPhantomColour* colourPointer = new HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* PelvisVisAtt = new G4VisAttributes(colour);
  PelvisVisAtt->SetForceSolid(wireFrame);
  logicPelvis->SetVisAttributes(PelvisVisAtt);

  G4cout << "Pelvis created !!!!!!" << G4endl;

  // Testing Pelvis Volume
  G4double PelvisVol = logicPelvis->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Pelvis = " << PelvisVol/cm3 << " cm^3" << G4endl;
  
  // Testing Pelvis Material
  G4String PelvisMat = logicPelvis->GetMaterial()->GetName();
  G4cout << "Material of Pelvis = " << PelvisMat << G4endl;
  
  // Testing Density
  G4double PelvisDensity = logicPelvis->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << PelvisDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double PelvisMass = (PelvisVol)*PelvisDensity;
  G4cout << "Mass of Pelvis = " << PelvisMass/gram << " g" << G4endl;

  VOL=PelvisVol;RHO=PelvisDensity;
  return physPelvis;
}
