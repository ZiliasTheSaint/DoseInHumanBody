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
#ifndef HumanPhantomEventAction_h
#define HumanPhantomEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <map>

class G4Event;

class HumanPhantomRunAction;
class HumanPhantomEventAction : public G4UserEventAction
{
public:
  HumanPhantomEventAction();
  ~HumanPhantomEventAction();

public:
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
  void AddRealEnergy(G4double de){phener +=de;};
  //--------------------
 /* void AddThyroidEnergy(G4double de)
  {thyroidEDEP +=de;};

  void AddLiverEnergy(G4double de)
  {liverEDEP +=de;};

  void AddHeartEnergy(G4double de)
  {heartEDEP +=de;};

  void AddLeftClavicleEnergy(G4double de)
  {leftClavicleEDEP +=de;};

  void AddRightClavicleEnergy(G4double de)
  {rightClavicleEDEP +=de;};

  void AddGallBladderEnergy(G4double de)
  {gallBladderEDEP +=de;};

  void AddEsophagusEnergy(G4double de)
  {esophagusEDEP +=de;};

  void AddSmallIntestineEnergy(G4double de)
  {smallIntestineEDEP +=de;};

  void AddLeftTesticleEnergy(G4double de)
  {leftTesticleEDEP +=de;};

  void AddRightTesticleEnergy(G4double de)
  {rightTesticleEDEP +=de;};

  void AddThymusEnergy(G4double de)
  {thymusEDEP +=de;};
  *///--------------------
private:
  void Fill(G4String bodypartName, G4double energyDeposit);
  void totalEventEnergyDeposit();
 
  G4int hitCollectionID; 
  std::map<std::string,G4double> energyTotal;
  
  G4double  phener;
  HumanPhantomRunAction*  runAct;

  /*G4double thyroidEDEP;G4double liverEDEP; 
  G4double heartEDEP;G4double leftClavicleEDEP;G4double rightClavicleEDEP;G4double gallBladderEDEP;
  G4double esophagusEDEP;//NA
  G4double smallIntestineEDEP;G4double leftTesticleEDEP;G4double rightTesticleEDEP;G4double thymusEDEP;

  G4bool doCollectionBasedScoring;//false if in thyroid or other non-geometrical organs...do not score in remainder volume i.e. head, trunk.
  */
};
#endif

    
