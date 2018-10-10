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
#ifndef HumanPhantomMessenger_h
#define HumanPhantomMessenger_h 1

class HumanPhantomConstruction;

class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADoubleAndUnit;

#include "G4UImessenger.hh"
#include "globals.hh"
#include <iostream>

class HumanPhantomMessenger: public G4UImessenger
{
public:
  HumanPhantomMessenger(HumanPhantomConstruction* myUsrPhtm);
  ~HumanPhantomMessenger();
    
  void SetNewValue(G4UIcommand* command, G4String newValue);

  //void AddBodyPart(G4String);	      // Set Body Parts Sensitivity

private:
  HumanPhantomConstruction*           myUserPhantom;

  G4UIdirectory*                 phantomDir;
  //G4UIdirectory*                 bpDir;

  G4UIcmdWithAString*            modelCmd; 
  G4UIcmdWithAString*            sexCmd;  
  G4UIcmdWithAString*            ageCmd;
  G4UIcmdWithADoubleAndUnit* heightCmd;
  G4UIcmdWithADoubleAndUnit* massCmd;
  G4UIcmdWithAString*            ageGroupCmd;
  //G4UIcmdWithAString*            bodypartCmd;
  G4UIcmdWithoutParameter*       endCmd;

  G4String                       bodypart;
  G4bool                         bps;

  G4UIcmdWithADoubleAndUnit* ftotal_diameterCmd;
  G4UIcmdWithADoubleAndUnit* ftotal_heightCmd;
  G4UIcmdWithADoubleAndUnit* fsourceTopCmd;
};

#endif

