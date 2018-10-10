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
#include "HumanPhantomEventAction.hh"
#include "HumanPhantomHit.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "HumanPhantomRunAction.hh"
#include "G4RunManager.hh"
#include "HumanPhantomPrimaryGeneratorAction.hh"

#include "HumanPhantomConstruction.hh"

HumanPhantomEventAction::HumanPhantomEventAction():
hitCollectionID(-1)
{ 
 runAct = (HumanPhantomRunAction*)G4RunManager::GetRunManager()->GetUserRunAction(); 
}
 
HumanPhantomEventAction::~HumanPhantomEventAction()
{
}

void HumanPhantomEventAction::BeginOfEventAction(const G4Event* evt)
{
	G4int pmodulo=HumanPhantomPrimaryGeneratorAction::printModulo;
	G4int evtNb = evt->GetEventID();
    if (evtNb%pmodulo == 0) { 
     G4cout << "\n---> Begin of event: " << evtNb 
		 <<" ;Incident energy: "<<G4BestUnit(HumanPhantomPrimaryGeneratorAction::incidentEnergy,"Energy") << G4endl;
	}
	if(!HumanPhantomConstruction::isMammo){

		energyTotal["logicalHead"]=0.;
		energyTotal["logicalSkull"]=0.;
		energyTotal["logicalBrain"]=0.;
		energyTotal["logicalTrunk"]=0.;
		energyTotal["logicalLeftLeg"]=0.;
		energyTotal["logicalRightLeg"]=0.; 
		energyTotal["logicalLeftArmBone"]=0.;
		energyTotal["logicalRightArmBone"]=0.;
		energyTotal["logicalLeftLegBone"]=0.;
		energyTotal["logicalRightLegBone"]=0.;

		energyTotal["logicalUpperSpine"]=0.;
		energyTotal["logicalLeftScapula"]=0.; 
		energyTotal["logicalRightScapula"]=0.;
		energyTotal["logicalLeftAdrenal"]=0.; 
		energyTotal["logicalRightAdrenal"]=0.;
		energyTotal["logicalMiddleLowerSpine"]=0.;
		energyTotal["logicalPelvis"]=0.;
		energyTotal["logicalStomach"]=0.;
		energyTotal["logicalUpperLargeIntestine"]=0.;
		energyTotal["logicalLowerLargeIntestine"]=0.;

		energyTotal["logicalRibCage"]=0.;
		energyTotal["logicalSpleen"]=0.;
		energyTotal["logicalPancreas"]=0.; 
		energyTotal["logicalLeftKidney"]=0.;
		energyTotal["logicalRightKidney"]=0.;
		energyTotal["logicalUrinaryBladder"]=0.;
		energyTotal["logicalLeftLung"]=0.;
		energyTotal["logicalRightLung"]=0.;	 
		energyTotal["logicalLeftOvary"]=0.;
		energyTotal["logicalRightOvary"]=0.;

		energyTotal["logicalUterus"]=0.;
		energyTotal["logicalLeftBreast"]=0.;
		energyTotal["logicalRightBreast"]=0.; 
 
		energyTotal["logicalThyroid"]=0.0;//!!!!!!!!!!
        energyTotal["logicalLiver"]=0.0;//!!!!!!!!!
		energyTotal["logicalHeart"]=0.0;//!!!!!!!!!!!!!
		
		energyTotal["logicalLeftTesticle"]=0.0;//!!!!!!!!!!!!!
		energyTotal["logicalRightTesticle"]=0.0;//!!!!!!!!!!!!!
		energyTotal["logicalFacialBones"]=0.;//////////////

		energyTotal["logicalLeftClavicle"]=0.0;
		energyTotal["logicalRightClavicle"]=0.0;
		energyTotal["logicalGallBladder"]=0.0;
		energyTotal["logicalEsophagus"]=0.0;
		energyTotal["logicalSmallIntestine"]=0.0;
		energyTotal["logicalThymus"]=0.0;

		energyTotal["logicalmaleGenitalia"]=0.0;//

		G4SDManager * SDman = G4SDManager::GetSDMpointer();  

		if (hitCollectionID==-1) {
			hitCollectionID = SDman->GetCollectionID("HumanPhantomCollection");
		}

		/*thyroidEDEP=0.;liverEDEP=0.;heartEDEP=0.;leftClavicleEDEP=0.;rightClavicleEDEP=0.;
		gallBladderEDEP=0.0;esophagusEDEP=0.;smallIntestineEDEP=0.;leftTesticleEDEP=0.;rightTesticleEDEP=0.;
		thymusEDEP=0.;*/
	}
	else{//if is mamo
		phener=0.;
	}
}
 
void HumanPhantomEventAction::EndOfEventAction(const G4Event* evt)
{  
	if(!HumanPhantomConstruction::isMammo){
		//==================this will be removed
		/*doCollectionBasedScoring =true;
		if (thyroidEDEP>0){
			Fill("logicalThyroid",thyroidEDEP);
			doCollectionBasedScoring =false;
		}
		if (liverEDEP>0){
			Fill("logicalLiver",liverEDEP);
			doCollectionBasedScoring =false;
		}
		if (heartEDEP>0){
			Fill("logicalHeart",heartEDEP);
			doCollectionBasedScoring =false;
		}
		if (leftClavicleEDEP>0){
			Fill("logicalLeftClavicle",leftClavicleEDEP);
			doCollectionBasedScoring =false;
		}
		if (rightClavicleEDEP>0){
			Fill("logicalRightClavicle",rightClavicleEDEP);
			doCollectionBasedScoring =false;
		}
		if (gallBladderEDEP>0){
			Fill("logicalGallBladder",gallBladderEDEP);
			doCollectionBasedScoring =false;
		}
		if (esophagusEDEP>0){
			Fill("logicalEsophagus",esophagusEDEP);//NA
			doCollectionBasedScoring =false;
		}
		if (smallIntestineEDEP>0){
			Fill("logicalSmallIntestine",smallIntestineEDEP);
			doCollectionBasedScoring =false;
		}
		if (leftTesticleEDEP>0){
			Fill("logicalLeftTesticle",leftTesticleEDEP);
			doCollectionBasedScoring =false;
		}
		if (rightTesticleEDEP>0){
			Fill("logicalRightTesticle",rightTesticleEDEP);
			doCollectionBasedScoring =false;
		}
		if (thymusEDEP>0){
			Fill("logicalThymus",thymusEDEP);
			doCollectionBasedScoring =false;
		}*/
		//==================================
		G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
 
		HumanPhantomHitsCollection* HC = 0;

		if (HCE)
			HC = (HumanPhantomHitsCollection*)(HCE->GetHC(hitCollectionID));

		if (HC)
		{
			G4int hitNumber = HC->entries();
			G4double edep =0;
			G4String bodyPart;
			G4double mass=1.;
			for (G4int i=0;i<hitNumber;i++) 
			{
				edep = (*HC)[i]->GetEdep();
				bodyPart = (*HC)[i]->GetBodyPartID();
				//if(doCollectionBasedScoring)//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NO..we can have dose in thyroid and spline..miss spline!		  
					Fill(bodyPart, edep);
	     
	        }
	   }
    	  
		

		//================================
       totalEventEnergyDeposit();
	}
	else{//if is mamo
		runAct->fillPerEvent(phener);
	}
}

void HumanPhantomEventAction:: Fill(G4String bodypartName, 
				      G4double energyDeposit)

{
    energyTotal[bodypartName] += energyDeposit;
}

void HumanPhantomEventAction::totalEventEnergyDeposit() 
{

    G4RunManager* runManager = G4RunManager::GetRunManager();
    HumanPhantomRunAction* pointerRun = (HumanPhantomRunAction*)(runManager->GetUserRunAction());

    std::map<std::string,G4double>::iterator i = energyTotal.begin();
    std::map<std::string,G4double>::iterator end = energyTotal.end();

    while(i!=end)
    {

       G4String bodypart = i->first;
       G4double energyDep = i->second;
      
       if(energyDep != 0.)
	   {
	       pointerRun->Fill(bodypart, energyDep);
	   }
       
	   i++;
    }
  
}
