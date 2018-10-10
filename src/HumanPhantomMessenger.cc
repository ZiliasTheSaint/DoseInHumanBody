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
#include "HumanPhantomMessenger.hh"
#include "HumanPhantomConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include "globals.hh"

#include "G4RunManager.hh"

HumanPhantomMessenger::HumanPhantomMessenger(HumanPhantomConstruction* myUsrPhtm)
  :myUserPhantom(myUsrPhtm),bps(false)
{ 
  phantomDir = new G4UIdirectory("/phantom/");
  phantomDir->SetGuidance("Set Your Phantom.");
  
  //bpDir = new G4UIdirectory("/bodypart/");
  //bpDir->SetGuidance("Add Body Part to Phantom");

  modelCmd = new G4UIcmdWithAString("/phantom/setPhantomModel",this);
  modelCmd->SetGuidance("Set Phantom: MIRD, MAMMO");//ORNLFemale, ORNLMale, MIX, MIRDHead, ORNLHead.");
  modelCmd->SetParameterName("phantomModel",true);
  modelCmd->SetDefaultValue("MIRDHead");
  modelCmd->SetCandidates("MIRD MIRDHead MAMMO");//ORNLFemale ORNLMale MIX MIRDHead ORNLHead");
  modelCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  sexCmd = new G4UIcmdWithAString("/phantom/setPhantomSex",this);
  sexCmd->SetGuidance("Set sex of Phantom: Male or Female.");
  sexCmd->SetParameterName("phantomSex",true);
  sexCmd->SetDefaultValue("Female");
  sexCmd->SetCandidates("Male Female");
  sexCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  ageCmd = new G4UIcmdWithAString("/phantom/setPhantomAge",this);
  ageCmd->SetGuidance("Set Age of Phantom.");
  ageCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  //bodypartCmd = new G4UIcmdWithAString("/bodypart/addBodyPart",this);
  //bodypartCmd->SetGuidance("Add a Body Part to Phantom");
  //bodypartCmd->SetParameterName("bpName",true);
  //bodypartCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  endCmd = new G4UIcmdWithoutParameter("/phantom/buildNewPhantom",this);
  endCmd->SetGuidance("Build your Phantom.");
  endCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fsourceTopCmd = new G4UIcmdWithADoubleAndUnit("/phantom/focusToBreastDistance",this);
  fsourceTopCmd->SetGuidance("Define focus to breast distance.");
  fsourceTopCmd->SetParameterName("Size",false);
  fsourceTopCmd->SetRange("Size>0.");
  fsourceTopCmd->SetUnitCategory("Length");
  fsourceTopCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ftotal_diameterCmd = new G4UIcmdWithADoubleAndUnit("/phantom/breast_total_diameter",this);
  ftotal_diameterCmd->SetGuidance("Define total diameter.");
  ftotal_diameterCmd->SetParameterName("Size",false);
  ftotal_diameterCmd->SetRange("Size>0.");
  ftotal_diameterCmd->SetUnitCategory("Length");
  ftotal_diameterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
  ftotal_heightCmd = new G4UIcmdWithADoubleAndUnit("/phantom/breast_total_height",this);
  ftotal_heightCmd->SetGuidance("Define total height.");
  ftotal_heightCmd->SetParameterName("Size",false);
  ftotal_heightCmd->SetRange("Size>0.");
  ftotal_heightCmd->SetUnitCategory("Length");
  ftotal_heightCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  heightCmd = new G4UIcmdWithADoubleAndUnit("/phantom/phantom_height",this);
  heightCmd->SetGuidance("Define phantom height.");
  heightCmd->SetParameterName("Size",false);
  heightCmd->SetRange("Size>0.");
  heightCmd->SetUnitCategory("Length");
  heightCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  massCmd = new G4UIcmdWithADoubleAndUnit("/phantom/phantom_mass",this);
  massCmd->SetGuidance("Define phantom weight.");
  massCmd->SetUnitCategory("Mass");
  massCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ageGroupCmd = new G4UIcmdWithAString("/phantom/setAgeGroup",this);
  ageGroupCmd->SetGuidance("Set age group.");
  ageGroupCmd->SetParameterName("phantomAgeGroup",true);
  ageGroupCmd->SetDefaultValue("adult");
  ageGroupCmd->SetCandidates("newborn age_1 age_5 age_10 age_15 adult");
  ageGroupCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
}

HumanPhantomMessenger::~HumanPhantomMessenger()
{
  delete  modelCmd;
  delete  sexCmd;
  delete  ageCmd;
  delete  massCmd;
  delete  heightCmd;
  delete ageGroupCmd;
  //delete  bodypartCmd;
  delete  endCmd;
  delete  phantomDir;
  //delete  bpDir;
  delete ftotal_diameterCmd;
  delete ftotal_heightCmd;
  delete fsourceTopCmd;
}

void HumanPhantomMessenger::SetNewValue(G4UIcommand* command,G4String newValue){ 

  if( command == modelCmd )
    { 
      myUserPhantom->SetPhantomModel(newValue); 
    }       
  if( command == sexCmd )
    { 
      myUserPhantom->SetPhantomSex(newValue); 
    }
  if( command == ageCmd )
    { 
      myUserPhantom->SetPhantomAge(newValue); 
    }
  //if( command == bodypartCmd )
  //  {
  //    AddBodyPart(newValue);
  //  }
  if( command == endCmd )
    { 
      G4cout << 
	" ****************>>>> NEW PHANTOM CONSTRUCTION <<<<***************** " 
	     << G4endl;
    }

  if( command == ftotal_diameterCmd ) {
    myUserPhantom
      ->SetTotal_diameter(ftotal_diameterCmd->GetNewDoubleValue(newValue));
  }
    
  if( command == ftotal_heightCmd ) {
    myUserPhantom
      ->SetTotal_height(ftotal_heightCmd->GetNewDoubleValue(newValue));
  }

  if( command == fsourceTopCmd ) {
    myUserPhantom
      ->SetPointSourceOrFrontalBeamToDetectorDistance(fsourceTopCmd->GetNewDoubleValue(newValue));
  }

  if( command == heightCmd ) {
    myUserPhantom
      ->SetPhantomHeightForScaling(heightCmd->GetNewDoubleValue(newValue));
  }
  if( command == massCmd ) {
    myUserPhantom
      ->SetPhantomMassForScaling(massCmd->GetNewDoubleValue(newValue));
  }

  if( command == ageGroupCmd )
    { 
      myUserPhantom->SetAgeGroup(newValue); 
    }
}

/*void  HumanPhantomMessenger::AddBodyPart(G4String newBodyPartSensitivity)
{

  char* str = new char[newBodyPartSensitivity.length()+1];

  strcpy(str, newBodyPartSensitivity.c_str()); 
  
  std::string bodypart = strtok(str," ");

  std::string sensitivity = strtok(NULL," ");

  if(sensitivity=="yes"){
    bps=true;
  }else{
    bps=false;
  }

  G4cout << " >>> Body Part = " << bodypart << "\n"
	 << " >>> Sensitivity = " << sensitivity << G4endl;

  myUserPhantom->SetBodyPartSensitivity(bodypart,bps);
}*/

