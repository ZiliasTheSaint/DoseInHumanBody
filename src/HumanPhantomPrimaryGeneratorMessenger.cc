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
// D. Fulea, National Institute of Public Health Bucharest, Cluj-Napoca Regional Center, Romania
// 

#include "HumanPhantomPrimaryGeneratorMessenger.hh"
#include "HumanPhantomPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"

HumanPhantomPrimaryGeneratorMessenger::HumanPhantomPrimaryGeneratorMessenger(HumanPhantomPrimaryGeneratorAction* primaryGun)
  :primary(primaryGun)
{ 
  beamCmd = new G4UIcmdWithAString("/gun/setBeam",this);
  beamCmd->SetGuidance("Choose the type of beam");
  beamCmd->SetGuidance("Choice : beamAlongX, beamAlongY, beamAlongZ, isotropicFlux, rectangleField, CTScan ");
  beamCmd->SetParameterName("choice",true);
  beamCmd->SetCandidates("beamAlongX beamAlongY beamAlongZ isotropicFlux rectangleField CTScan");
  beamCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  eventDir = new G4UIdirectory("/event/");
  eventDir->SetGuidance("event control");

  sourceDir = new G4UIdirectory("/source/");
  sourceDir->SetGuidance("radiation source control");

  PrintCmd = new G4UIcmdWithAnInteger("/event/printModulo",this);
  PrintCmd->SetGuidance("Print events modulo n");
  PrintCmd->SetParameterName("EventNb",false);
  PrintCmd->SetRange("EventNb>0");

  xfieldDir = new G4UIdirectory("/xfield/");
  xfieldDir->SetGuidance("XField control");

  fieldCenterYCmd = new G4UIcmdWithADoubleAndUnit("/xfield/centerY",this);
  fieldCenterYCmd->SetGuidance("Define field center Y position.");
  fieldCenterYCmd->SetParameterName("Size",false);
  //fieldCenterYCmd->SetRange("Size>0.");
  fieldCenterYCmd->SetUnitCategory("Length");
  fieldCenterYCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fieldCenterXCmd = new G4UIcmdWithADoubleAndUnit("/xfield/centerX",this);
  fieldCenterXCmd->SetGuidance("Define field center X position.");
  fieldCenterXCmd->SetParameterName("Size",false);
  //fieldCenterXCmd->SetRange("Size>0.");
  fieldCenterXCmd->SetUnitCategory("Length");
  fieldCenterXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fieldWidthCmd = new G4UIcmdWithADoubleAndUnit("/xfield/width",this);
  fieldWidthCmd->SetGuidance("Define field width.");
  fieldWidthCmd->SetParameterName("Size",false);
  fieldWidthCmd->SetRange("Size>0.");
  fieldWidthCmd->SetUnitCategory("Length");
  fieldWidthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fieldHeightCmd = new G4UIcmdWithADoubleAndUnit("/xfield/height",this);
  fieldHeightCmd->SetGuidance("Define field height.");
  fieldHeightCmd->SetParameterName("Size",false);
  fieldHeightCmd->SetRange("Size>0.");
  fieldHeightCmd->SetUnitCategory("Length");
  fieldHeightCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  azimuthalRotationAngleCmd = new G4UIcmdWithADoubleAndUnit("/xfield/azimuthalAngle",this);
  azimuthalRotationAngleCmd->SetGuidance("Define azimuthal angle of field central axis relative to phantom symmetry axis.");
  azimuthalRotationAngleCmd->SetParameterName("Angle",false);  
  azimuthalRotationAngleCmd->SetUnitCategory("Angle");
  azimuthalRotationAngleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  polarRotationAngleCmd = new G4UIcmdWithADoubleAndUnit("/xfield/polarAngle",this);
  polarRotationAngleCmd->SetGuidance("Define polar angle of field central axis relative to phantom symmetry axis.");
  polarRotationAngleCmd->SetParameterName("Angle",false);  
  polarRotationAngleCmd->SetUnitCategory("Angle");
  polarRotationAngleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  sliceThicknessCmd = new G4UIcmdWithADoubleAndUnit("/xfield/sliceThickness",this);
  sliceThicknessCmd->SetGuidance("Define slice thickness for CT scan.");
  sliceThicknessCmd->SetParameterName("Size",false); 
  sliceThicknessCmd->SetRange("Size>0.");
  sliceThicknessCmd->SetUnitCategory("Length");
  sliceThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  angleIncrementCmd = new G4UIcmdWithADoubleAndUnit("/xfield/angleIncrement",this);
  angleIncrementCmd->SetGuidance("Define angle increment for CT scan.");
  angleIncrementCmd->SetParameterName("Size",false); 
  angleIncrementCmd->SetRange("Size>0.");
  angleIncrementCmd->SetUnitCategory("Angle");
  angleIncrementCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  isFanBeamCmd = new G4UIcmdWithAString("/xfield/SetFanBeam",this);
  isFanBeamCmd->SetGuidance("Set on if CTScan is fan beam or off if pencil beam.");
  isFanBeamCmd->SetGuidance("  Choice : on(default), off");
  isFanBeamCmd->SetParameterName("choice",true);
  isFanBeamCmd->SetDefaultValue("on");
  isFanBeamCmd->SetCandidates("on off");
  isFanBeamCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  isHelicalScanCmd = new G4UIcmdWithAString("/xfield/SetHelicalScan",this);
  isHelicalScanCmd->SetGuidance("Set on if CTScan is helical or off if CTScan consists only of one full rotation followed by a table-feed translation and repeat.");
  isHelicalScanCmd->SetGuidance("  Choice : on(default), off");
  isHelicalScanCmd->SetParameterName("choice",true);
  isHelicalScanCmd->SetDefaultValue("on");
  isHelicalScanCmd->SetCandidates("on off");
  isHelicalScanCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  isHalfFieldCmd = new G4UIcmdWithAString("/xfield/SetHalfField",this);
  isHalfFieldCmd->SetGuidance("Set on if CTScan is asymetric meaning half field or off if normal scan.");
  isHalfFieldCmd->SetGuidance("  Choice : on(default), off");
  isHalfFieldCmd->SetParameterName("choice",true);
  isHalfFieldCmd->SetDefaultValue("on");
  isHalfFieldCmd->SetCandidates("on off");
  isHalfFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  isDentalPanoramicCmd = new G4UIcmdWithAString("/xfield/SetDentalPanoramic",this);
  isDentalPanoramicCmd->SetGuidance("Set on if dental panoramic with 180 rotation or off if normal CTScan.");
  isDentalPanoramicCmd->SetGuidance("  Choice : on(default), off");
  isDentalPanoramicCmd->SetParameterName("choice",true);
  isDentalPanoramicCmd->SetDefaultValue("on");
  isDentalPanoramicCmd->SetCandidates("on off");
  isDentalPanoramicCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  isSpectrumCmd = new G4UIcmdWithAString("/xfield/isSpectrum?",this);
  isSpectrumCmd->SetGuidance("Set yes if we have an energy spectrum, and no otherwise.");
  isSpectrumCmd->SetGuidance("  Choice : yes(default), no");
  isSpectrumCmd->SetParameterName("choice",true);
  isSpectrumCmd->SetDefaultValue("yes");
  isSpectrumCmd->SetCandidates("yes no");
  isSpectrumCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  kvCmd = new G4UIcmdWithAString("/xfield/kVp",this);
  kvCmd->SetGuidance("Set Xray tube kilovoltage.");  
  kvCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  filtrationCmd = new G4UIcmdWithAString("/xfield/filtration",this);
  filtrationCmd->SetGuidance("Set Xray tube filtration [mmAl].");  
  filtrationCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  anodMaterialCmd = new G4UIcmdWithAString("/xfield/anodMaterial",this);
  anodMaterialCmd->SetGuidance("Set Xray tube anodMaterial.");  
  anodMaterialCmd->SetCandidates("W MO RH");
  anodMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  anodAngleCmd = new G4UIcmdWithAString("/xfield/anodAngle",this);
  anodAngleCmd->SetGuidance("Set Xray tube anode angle [deg].");  
  anodAngleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  rippleCmd = new G4UIcmdWithAString("/xfield/ripple",this);
  rippleCmd->SetGuidance("Set Xray tube waveform ripple [%].");  
  rippleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  FCACmd = new G4UIcmdWithADoubleAndUnit("/xfield/focusToPhantomSymmetryAxisDistance",this);
  FCACmd->SetGuidance("Define focus to phantom symmetry axis distance.");
  FCACmd->SetParameterName("Size",false); 
  FCACmd->SetRange("Size>0.");
  FCACmd->SetUnitCategory("Length");
  FCACmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DAP_uGymm2Cmd = new G4UIcmdWithAString("/xfield/DAP[uGymm2]",this);
  DAP_uGymm2Cmd->SetGuidance("Set free in air dose area product.");  
  DAP_uGymm2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  CTDI_uGyCmd = new G4UIcmdWithAString("/xfield/CTDI[uGy]",this);
  CTDI_uGyCmd->SetGuidance("Set CTDI (or dose/rotation free in air).");  
  CTDI_uGyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  MASCmd = new G4UIcmdWithAString("/xfield/mAs",this);
  MASCmd->SetGuidance("Set mAs.");  
  MASCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  useMASCmd = new G4UIcmdWithAString("/xfield/use_mAs_for_DAP_calculation?",this);
  useMASCmd->SetGuidance("Set yes if want to use mAs (and Kerma from spectrum) for DAP calculation.");
  useMASCmd->SetGuidance("  Choice : yes(default), no");
  useMASCmd->SetParameterName("choice",true);
  useMASCmd->SetDefaultValue("yes");
  useMASCmd->SetCandidates("yes no");
  useMASCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  pitchCmd = new G4UIcmdWithAString("/xfield/pitch",this);
  pitchCmd->SetGuidance("Set pitch factor.");  
  pitchCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  quantaCmd = new G4UIcmdWithAnInteger("/source/quantaMode",this);
  quantaCmd->SetGuidance("Quanta emmission mode requires knowledge of source activity or kerma measurements or number of quantas");
  quantaCmd->SetParameterName("quantaMode",false);
  quantaCmd->SetRange("quantaMode>=0");

  activityCmd = new G4UIcmdWithAString("/source/activity",this);
  activityCmd->SetGuidance("Set source activity in Bq.");  
  activityCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  exposureTimeCmd = new G4UIcmdWithAString("/source/exposureTime",this);
  exposureTimeCmd->SetGuidance("Set target exposure time in seconds.");  
  exposureTimeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  isIsotropicSourceCmd = new G4UIcmdWithAString("/source/isotropicSource",this);
  isIsotropicSourceCmd->SetGuidance("Set on if isotropic radiation source.");
  isIsotropicSourceCmd->SetGuidance("  Choice : on(default), off");
  isIsotropicSourceCmd->SetParameterName("choice",true);
  isIsotropicSourceCmd->SetDefaultValue("on");
  isIsotropicSourceCmd->SetCandidates("on off");
  isIsotropicSourceCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  radiationYieldCmd = new G4UIcmdWithAString("/source/radiationYield",this);
  radiationYieldCmd->SetGuidance("Set radiation yield.");  
  radiationYieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  emmissionAreaCmd = new G4UIcmdWithAString("/source/emmissionArea",this);
  emmissionAreaCmd->SetGuidance("Set emmission area if not isotropic source ...in mm2.");  
  emmissionAreaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  kermaCmd = new G4UIcmdWithAString("/source/kerma",this);
  kermaCmd->SetGuidance("Set kerma in Gy at entry surface.");  
  kermaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  mentrCmd = new G4UIcmdWithAString("/source/mentr",this);
  mentrCmd->SetGuidance("Set mass-energy transfer coefficient for photons at incident energy in medium ...cm2/g.");  
  mentrCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  particleCounterCmd = new G4UIcmdWithAString("/source/particleCounter",this);
  particleCounterCmd->SetGuidance("Set the number of quantas hitting the target.");  
  particleCounterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

HumanPhantomPrimaryGeneratorMessenger::~HumanPhantomPrimaryGeneratorMessenger()
{
  delete beamCmd;
  delete fieldCenterYCmd;
  delete fieldCenterXCmd;
  delete fieldWidthCmd;
  delete fieldHeightCmd;
  delete xfieldDir;
  delete PrintCmd;
  delete eventDir;
  delete sourceDir;
  delete quantaCmd;
  delete activityCmd;
  delete exposureTimeCmd;
  delete radiationYieldCmd;
  delete isIsotropicSourceCmd;
  delete emmissionAreaCmd;
  delete kermaCmd;
  delete mentrCmd;
  delete particleCounterCmd;
  delete azimuthalRotationAngleCmd;
  delete polarRotationAngleCmd;
  delete sliceThicknessCmd;
  delete angleIncrementCmd;
  delete isFanBeamCmd;
  delete isHelicalScanCmd;
  delete isHalfFieldCmd;
  delete isDentalPanoramicCmd;
  delete isSpectrumCmd;
  delete kvCmd;
  delete filtrationCmd;
  delete anodAngleCmd;
  delete rippleCmd;
  delete anodMaterialCmd;
  delete DAP_uGymm2Cmd;
  delete CTDI_uGyCmd;
  delete FCACmd;
  delete pitchCmd;
}

void HumanPhantomPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if( command == beamCmd )
   { primary->SetBeam(newValue);
   }

  if(command == PrintCmd)
    {primary->SetPrintModulo(PrintCmd->GetNewIntValue(newValue));}

  if( command == fieldCenterYCmd ) {
    primary
      ->SetFieldCenterY(fieldCenterYCmd->GetNewDoubleValue(newValue));
  }

  if( command == fieldCenterXCmd ) {
    primary
      ->SetFieldCenterX(fieldCenterXCmd->GetNewDoubleValue(newValue));
  }

  if( command == fieldWidthCmd ) {
    primary
      ->SetFieldWidth(fieldWidthCmd->GetNewDoubleValue(newValue));
  }

  if( command == fieldHeightCmd ) {
    primary
      ->SetFieldHeight(fieldHeightCmd->GetNewDoubleValue(newValue));
  }

  if( command == azimuthalRotationAngleCmd ) {
    primary
      ->SetAzimuthalRotationAngle(azimuthalRotationAngleCmd->GetNewDoubleValue(newValue));
  }

  if( command == polarRotationAngleCmd ) {
    primary
      ->SetPolarRotationAngle(polarRotationAngleCmd->GetNewDoubleValue(newValue));
  }

  if( command == isFanBeamCmd )
   { primary->SetFanBeam(newValue);
   }

  if( command == isHelicalScanCmd )
   { primary->SetHelicalScan(newValue);
   }

  if( command == isHalfFieldCmd )
   { primary->SetHalfField(newValue);
   }

  if( command == isDentalPanoramicCmd )
   { primary->SetDentalPanoramic(newValue);
   }

  if( command == sliceThicknessCmd ) {
    primary
      ->SetSliceThickness(sliceThicknessCmd->GetNewDoubleValue(newValue));
  }

    if( command == angleIncrementCmd ) {
    primary
      ->SetAngleIncrement(angleIncrementCmd->GetNewDoubleValue(newValue));
  }

	if( command == isSpectrumCmd )
   { primary->SetIfSpectrum(newValue);
   }

	if( command == kvCmd )
   { primary->SetKv(newValue);
   }

	if( command == filtrationCmd )
   { primary->SetFiltration(newValue);
   }

	if( command == anodMaterialCmd )
   { primary->SetAnodMaterial(newValue);
   }

	if( command == anodAngleCmd )
   { primary->SetAnodAngle(newValue);
   }

	if( command == rippleCmd )
   { primary->SetRipple(newValue);
   }

	if( command == FCACmd ) {
    primary
      ->SetFocusToPhantomCentralAxisDistance(FCACmd->GetNewDoubleValue(newValue));
  }

	if( command == DAP_uGymm2Cmd )
   { primary->SetDAP_uGymm2(newValue);
   }

	if( command == CTDI_uGyCmd )
   { primary->SetCTDI(newValue);
   }

	if( command == MASCmd )
   { primary->SetMAS(newValue);
   }

	if( command == useMASCmd )
   { primary->SetIfUseMASforDAP(newValue);
   }

	if( command == pitchCmd )
   { primary->SetPitch(newValue);
   }

	//==============
	if(command == quantaCmd)
    {primary->SetIQUANTA(quantaCmd->GetNewIntValue(newValue));}

	if( command == activityCmd )
   { primary->SetActivity_Bq(newValue);
   }

	if( command == exposureTimeCmd )
   { primary->SetExposureTime_s(newValue);
   }

	if( command == radiationYieldCmd )
   { primary->SetRadiationYield(newValue);
   }

	if( command == isIsotropicSourceCmd )
   { primary->SetIsotropicSource(newValue);
   }

	if( command == emmissionAreaCmd )
   { primary->SetEmmissionArea_mm2(newValue);
   }

	if( command == kermaCmd )
   { primary->SetKerma_Gy(newValue);
   }

	if( command == mentrCmd )
   { primary->SetMassEnergyTransferCoefficient_cm2Perg(newValue);
   }

	if( command == particleCounterCmd )
   { primary->SetParticleCounter(newValue);
   }
}


