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
#ifndef MIRDLeftOvary_h
#define MIRDLeftOvary_h 1

#include "G4VPhysicalVolume.hh"
#include "VOrgan.hh"

class G4VPhysicalVolume;

class MIRDLeftOvary: public VOrgan
{
public:

  MIRDLeftOvary();
  ~MIRDLeftOvary();
  G4VPhysicalVolume* Construct(const G4String&, G4VPhysicalVolume*,
				    const G4String&, G4bool, G4bool);
  G4double getCubicVolume(){return VOL;};
  G4double getDensity(){return RHO;};

private:
  G4double VOL;
  G4double RHO;

};
#endif
