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
//
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
//D. Fulea, National Institute of Public Health Bucharest, Cluj-Napoca Regional Center, Romania
//
#ifndef MIRDBodyFactory_h
#define MIRDBodyFactory_h 1

#include "VBodyFactory.hh"
#include "VOrgan.hh"
#include <map>
class VBodyFactory;
class G4VPhysicalVolume;
class VOrgan;

class MIRDBodyFactory: public VBodyFactory
{
public:
  MIRDBodyFactory();
 ~MIRDBodyFactory();
  G4VPhysicalVolume* CreateOrgan(const G4String&, G4VPhysicalVolume*, 
				 const G4String&,G4bool,G4bool);

  G4double GetBodyVOL(G4String body){
	  return bodyVOL[body];
  };
  G4double GetBodyRHO(G4String body){
	  return bodyRHO[body];
  };

  void SetScaleXY(G4double val);
  void SetScaleZ(G4double val);
  void SetAgeGroup(G4String val);
protected:
  G4double scaleXY;//scale on phantom transversal plane, x-y plane in MIRD5 description
  G4double scaleZ;//scale on phantom height, z axis in MIRD5 description.
  G4String ageGroup;
private:
  std::map<std::string,VOrgan*> organ; 
  //-----------------------------------------
  std::map<std::string,G4double> bodyVOL;
  std::map<std::string,G4double> bodyRHO;
};
#endif
