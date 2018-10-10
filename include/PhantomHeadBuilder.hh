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

#ifndef PhantomHeadBuilder_h
#define PhantomHeadBuilder_h 1

#include "G4VPhysicalVolume.hh"
#include "BasePhantomBuilder.hh"
#include "globals.hh"

class BasePhantomBuilder;
class G4VPhysicalVolume;
class VBodyFactory;
class PhantomHeadBuilder: public BasePhantomBuilder
{
public:
  PhantomHeadBuilder();
  ~PhantomHeadBuilder();

  void BuildHead(const G4String&,G4bool,G4bool);G4double getMassOfHead();
  void BuildSkull(const G4String&,G4bool,G4bool);G4double getMassOfSkull();
  void BuildBrain(const G4String&,G4bool,G4bool);G4double getMassOfBrain();
  void BuildFacial(const G4String&,G4bool,G4bool);

  void SetModel(G4String);
  void SetMotherVolume(G4VPhysicalVolume*);

 
  G4VPhysicalVolume* GetPhantom();
  //------------------
  G4double getMassOfFacial(){return volOfFacial*rhoOfFacial;};

  G4double getMassOfTrunk(){return 1.0;};
  G4double getMassOfLeftLeg(){return 1.0;};
  G4double getMassOfRightLeg(){return 1.0;};
  G4double getMassOfUpperSpine(){return 1.0;};
  G4double getMassOfMiddleLowerSpine(){return 1.0;};
  G4double getMassOfLeftLegBone(){return 1.0;};
  G4double getMassOfRightLegBone(){return 1.0;};
  G4double getMassOfLeftArmBone(){return 1.0;};
  G4double getMassOfRightArmBone()  {return 1.0;};
  G4double getMassOfRibCage(){return 1.0;};
  G4double getMassOfPelvis(){return 1.0;};
  G4double getMassOfLeftScapula(){return 1.0;};
  G4double getMassOfRightScapula(){return 1.0;};
  G4double getMassOfLeftAdrenal(){return 1.0;};
  G4double getMassOfRightAdrenal(){return 1.0;};
  G4double getMassOfHeart(){return 1.0;};
  G4double getMassOfLeftLung(){return 1.0;};
  G4double getMassOfRightLung(){return 1.0;};
  G4double getMassOfStomach(){return 1.0;};
  G4double getMassOfUpperLargeIntestine(){return 1.0;};
  G4double getMassOfLowerLargeIntestine(){return 1.0;};
  G4double getMassOfLeftKidney(){return 1.0;};
  G4double getMassOfRightKidney(){return 1.0;};
  G4double getMassOfLiver(){return 1.0;};
  G4double getMassOfPancreas(){return 1.0;};
  G4double getMassOfSpleen(){return 1.0;};
  G4double getMassOfUrinaryBladder(){return 1.0;};
  G4double getMassOfThyroid(){return 1.0;};

  G4double getMassOfLeftBreast(){return 1.0;};
  G4double getMassOfRightBreast(){return 1.0;};
  G4double getMassOfLeftOvary(){return 1.0;};
  G4double getMassOfRightOvary(){return 1.0;};
  G4double getMassOfUterus(){return 1.0;};

  G4double getMassOfEsophagus(){return 1.0;};
  G4double getMassOfLeftClavicle(){return 1.0;};
  G4double getMassOfRightClavicle(){return 1.0;};
  G4double getMassOfGallBladder(){return 1.0;};
  G4double getMassOfSmallIntestine(){return 1.0;};
  G4double getMassOfThymus(){return 1.0;};
  G4double getMassOfLeftTesticle(){return 1.0;};
  G4double getMassOfRightTesticle(){return 1.0;};

protected: 
  VBodyFactory* body;

  G4String model;

  G4VPhysicalVolume* motherVolume;
  G4VPhysicalVolume* headVolume;

  G4double volOfHead;G4double rhoOfHead;
  G4double volOfSkull;G4double rhoOfSkull;
  G4double volOfBrain;G4double rhoOfBrain;
  G4double volOfFacial;G4double rhoOfFacial;
};
#endif
