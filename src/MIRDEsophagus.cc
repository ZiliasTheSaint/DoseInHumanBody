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
#include "MIRDEsophagus.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "HumanPhantomMaterial.hh"
#include "G4SDManager.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4Ellipsoid.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include "G4EllipticalTube.hh"
#include "G4Box.hh"
#include <cmath>

#include "MirdMapConstants.hh"

MIRDEsophagus::MIRDEsophagus()
{
	VOL=45.*cm3;//??

	G4double d = 0.9869 *g/cm3;//softTissue
	RHO=d;
}

MIRDEsophagus::~MIRDEsophagus()
{
}

G4VPhysicalVolume* MIRDEsophagus::Construct(const G4String&,G4VPhysicalVolume*,
                                          const G4String&, G4bool, G4bool)
{
  //VOL = VOL*scaleXY*scaleXY*scaleZ;//scalexy^2*scaleZ=M/M0=V/V0
  
  MirdMapConstants* mmc = new MirdMapConstants();
	mmc->initEsophagus();
  std::map<std::string,esophagus_struct> esp = mmc->GetEsophagusMap();    
  G4double a = scaleXY*esp[ageGroup].a;  
  G4double b = scaleXY*esp[ageGroup].b;
  G4double d = scaleXY*esp[ageGroup].d;    
  G4double z2 = scaleZ*esp[ageGroup].z2;
  G4double z3 = scaleZ*esp[ageGroup].z3;
  G4double r = scaleXY*esp[ageGroup].r;
  G4double xp1 =scaleXY*esp[ageGroup].x1p;
  G4double xp2 =scaleXY*esp[ageGroup].x2p;
  delete mmc;

  double pi =3.14159265359;
  G4double vol1 = pi*(z3-z2)*(a*b-(a-d)*(b-d));
  G4double vol2 = pi*r*r*(xp2-xp1);

  VOL=vol1+vol2;
  //G4cout<<"pi= "<<pi<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<G4endl;

  return 0;
}
