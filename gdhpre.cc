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
// $Id: exampleN03.cc,v 1.39 2010-12-01 05:56:17 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

#include "QGSP_BIC_EMY.hh"

#include "HumanPhantomConstruction.hh"
#include "HumanPhantomPrimaryGeneratorAction.hh"
#include "HumanPhantomSteppingAction.hh"
#include "HumanPhantomEventAction.hh"
#include "HumanPhantomRunAction.hh"
#include "SteppingVerbose.hh"
#include "HumanPhantomPhysicsList.hh"
#include "MyModularPhysicsList.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc,char** argv)
{	
  // Choose the Random engine
  //
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
       
  // User Verbose output class
  //
  G4VSteppingVerbose::SetInstance(new SteppingVerbose);

  // Construct the default run manager
  //
  G4RunManager * runManager = new G4RunManager;

  // Set mandatory initialization classes
  //
  HumanPhantomConstruction* userPhantom = new HumanPhantomConstruction();
  runManager->SetUserInitialization(userPhantom);
  
  // Physics list
  //G4VModularPhysicsList* physicsList = new QGSP_BIC_EMY;
  //physicsList->SetVerboseLevel(1);
  //runManager->SetUserInitialization(physicsList);
  //runManager->SetUserInitialization(new HumanPhantomPhysicsList);   
  runManager->SetUserInitialization(new MyModularPhysicsList());   

  // Set user action classes  
  runManager->SetUserAction(new HumanPhantomPrimaryGeneratorAction);

#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  
  visManager->Initialize();  
#endif

  runManager->SetUserAction(new HumanPhantomRunAction);  
  HumanPhantomEventAction* eventAction = new HumanPhantomEventAction();
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(new HumanPhantomSteppingAction()); 
  
  // Initialize G4 kernel
  //runManager->Initialize();
  //INITIALIZATION IS CALLED FROM run.mac AFTER PHANTOM IS SET!!!
  //The vis.mac for vizualization is called from run.mac also!!!

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

//#ifdef G4VIS_USE
  //    UImanager->ApplyCommand("/control/execute vis.mac");
  //The vis.mac for vizualization is called from run.mac also!!!

  if (argc!=1)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
	  //std::cin.get();//hold to press enter!!
    }
  else
    {  // interactive mode : define UI session
		
	  #ifdef G4UI_USE
       G4cout << " UI session starts ..." << G4endl;
       G4UIExecutive* ui = new G4UIExecutive(argc, argv);
       UImanager->ApplyCommand("/control/execute run.mac");     
       ui->SessionStart();
       delete ui;
	 #endif
	   
  }  		
  
//#endif      
 
  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

