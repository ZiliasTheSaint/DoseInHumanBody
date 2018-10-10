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

#ifndef FemaleBuilder_h
#define FemaleBuilder_h 1

#include "PhantomBuilder.hh"

class FemaleBuilder: public PhantomBuilder
{
public:
  FemaleBuilder();
  ~FemaleBuilder();

  void BuildLeftBreast(const G4String&, G4bool, G4bool);
  G4double getMassOfLeftBreast(){return volOfLeftBreast*rhoOfLeftBreast;};
  void BuildRightBreast(const G4String&, G4bool, G4bool);
  G4double getMassOfRightBreast(){return volOfRightBreast*rhoOfRightBreast;};
  void BuildLeftOvary(const G4String&, G4bool, G4bool);
  G4double getMassOfLeftOvary(){return volOfLeftOvary*rhoOfLeftOvary;};
  void BuildRightOvary(const G4String&, G4bool, G4bool);
  G4double getMassOfRightOvary(){return volOfRightOvary*rhoOfRightOvary;};
  void BuildUterus(const G4String&, G4bool, G4bool);
  G4double getMassOfUterus(){return volOfUterus*rhoOfUterus;};

  G4double getMassOfTrunk(){
	  return rhoOfTrunk*(volOfTrunk-
		  volOfLeftArmBone-volOfRightArmBone-volOfLeftScapula-volOfRightScapula-
	      volOfMiddleLowerSpine-volOfPelvis-volOfUrinaryBladder-
	volOfHeart-volOfLeftLung-volOfRightLung-volOfLiver-
	      volOfStomach-volOfRibCage-volOfSpleen-volOfUpperLargeIntestine-volOfLowerLargeIntestine-
	      volOfLeftKidney-volOfRightKidney-volOfLeftAdrenal-volOfRightAdrenal-volOfPancreas-
		  volOfLeftOvary-volOfRightOvary-volOfUterus);
  };//@@@@@@@@@@@@
  
  G4double getMassOfLeftTesticle(){return 1.0;};
  G4double getMassOfRightTesticle(){return 1.0;};

  protected:
  G4double volOfLeftBreast;G4double rhoOfLeftBreast;
  G4double volOfRightBreast;G4double rhoOfRightBreast;

  G4double volOfLeftOvary;G4double rhoOfLeftOvary;
  G4double volOfRightOvary;G4double rhoOfRightOvary;
  G4double volOfUterus;G4double rhoOfUterus;

};
#endif
