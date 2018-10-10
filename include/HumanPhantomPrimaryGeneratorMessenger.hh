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

#ifndef HumanPhantomPrimaryGeneratorMessenger_h
#define HumanPhantomPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class HumanPhantomPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;

class HumanPhantomPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
  HumanPhantomPrimaryGeneratorMessenger(HumanPhantomPrimaryGeneratorAction*);
  ~HumanPhantomPrimaryGeneratorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  HumanPhantomPrimaryGeneratorAction*    primary; 
  G4UIcmdWithAString*               beamCmd;

  G4UIdirectory*          xfieldDir;
  G4UIdirectory*        eventDir; 
  G4UIdirectory*        sourceDir;

  G4UIcmdWithADoubleAndUnit* fieldCenterYCmd;
  G4UIcmdWithADoubleAndUnit* fieldCenterXCmd;
  G4UIcmdWithADoubleAndUnit* fieldHeightCmd;
  G4UIcmdWithADoubleAndUnit* fieldWidthCmd;
  G4UIcmdWithADoubleAndUnit* azimuthalRotationAngleCmd;
  G4UIcmdWithADoubleAndUnit* polarRotationAngleCmd;

  G4UIcmdWithADoubleAndUnit* sliceThicknessCmd;
  G4UIcmdWithADoubleAndUnit* angleIncrementCmd;
  G4UIcmdWithAString*        isFanBeamCmd;
  G4UIcmdWithAString*        isHelicalScanCmd;
  G4UIcmdWithAString*        isHalfFieldCmd;
  G4UIcmdWithAString*        isDentalPanoramicCmd;
  G4UIcmdWithAString*        isSpectrumCmd;
  //===============
  G4UIcmdWithAString*        kvCmd;
  G4UIcmdWithAString*        filtrationCmd;
  G4UIcmdWithAString*        anodAngleCmd;
  G4UIcmdWithAString*        rippleCmd;
  G4UIcmdWithAString*        anodMaterialCmd;
  //=================
  G4UIcmdWithAString*        DAP_uGymm2Cmd;
  G4UIcmdWithAString*        CTDI_uGyCmd;
  G4UIcmdWithADoubleAndUnit* FCACmd;

  G4UIcmdWithAString*        MASCmd;
  G4UIcmdWithAString*        useMASCmd;
  G4UIcmdWithAString*        pitchCmd;
  G4UIcmdWithAnInteger* PrintCmd;

  G4UIcmdWithAnInteger* quantaCmd;
  G4UIcmdWithAString*        activityCmd;
  G4UIcmdWithAString*        exposureTimeCmd;
  G4UIcmdWithAString*        radiationYieldCmd;
  G4UIcmdWithAString*        isIsotropicSourceCmd;
  G4UIcmdWithAString*        emmissionAreaCmd;
  G4UIcmdWithAString*        kermaCmd;
  G4UIcmdWithAString*        mentrCmd;
  G4UIcmdWithAString*        particleCounterCmd;
};

#endif

