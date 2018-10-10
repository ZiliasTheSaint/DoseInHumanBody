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

#ifndef PhantomBuilder_h
#define PhantomBuilder_h 1

#include "G4VPhysicalVolume.hh"
#include "BasePhantomBuilder.hh"
#include "globals.hh"

class BasePhantomBuilder;
class G4VPhysicalVolume;
class VBodyFactory;
class PhantomBuilder: public BasePhantomBuilder
{
public:
  PhantomBuilder();
  ~PhantomBuilder();

  void BuildHead(const G4String&,G4bool,G4bool);
  void BuildTrunk(const G4String&,G4bool,G4bool);
  void BuildLeftLeg(const G4String&,G4bool,G4bool);
  void BuildRightLeg(const G4String&,G4bool,G4bool);

  void BuildUpperSpine(const G4String&,G4bool,G4bool);
    G4double getMassOfUpperSpine(){return volOfUpperSpine*rhoOfUpperSpine;};
  void BuildMiddleLowerSpine(const G4String&,G4bool,G4bool);
    G4double getMassOfMiddleLowerSpine(){return volOfMiddleLowerSpine*rhoOfMiddleLowerSpine;};
  void BuildLeftLegBone(const G4String&,G4bool,G4bool);
	G4double getMassOfLeftLegBone(){return volOfLeftLegBone*rhoOfLeftLegBone;};
  void BuildRightLegBone(const G4String&,G4bool,G4bool);
    G4double getMassOfRightLegBone(){return volOfRightLegBone*rhoOfRightLegBone;};
  void BuildLeftArmBone(const G4String&,G4bool,G4bool);
    G4double getMassOfLeftArmBone(){return volOfLeftArmBone*rhoOfLeftArmBone;};
  void BuildRightArmBone(const G4String&,G4bool,G4bool);
    G4double getMassOfRightArmBone(){return volOfRightArmBone*rhoOfRightArmBone;};
 
  void BuildSkull(const G4String&,G4bool,G4bool);
    G4double getMassOfSkull(){return volOfSkull*rhoOfSkull;};
  void BuildFacial(const G4String&,G4bool,G4bool);
    G4double getMassOfFacial(){return volOfFacial*rhoOfFacial;};

  void BuildRibCage(const G4String&,G4bool,G4bool);
    G4double getMassOfRibCage(){return volOfRibCage*rhoOfRibCage;};
  void BuildPelvis(const G4String&,G4bool,G4bool);
    G4double getMassOfPelvis(){return volOfPelvis*rhoOfPelvis;};
  void BuildLeftScapula(const G4String&,G4bool,G4bool);
    G4double getMassOfLeftScapula(){return volOfLeftScapula*rhoOfLeftScapula;};
  void BuildRightScapula(const G4String&,G4bool,G4bool);
    G4double getMassOfRightScapula(){return volOfRightScapula*rhoOfRightScapula;};
  void BuildLeftAdrenal(const G4String&,G4bool,G4bool);
    G4double getMassOfLeftAdrenal(){return volOfLeftAdrenal*rhoOfLeftAdrenal;};
  void BuildRightAdrenal(const G4String&,G4bool,G4bool);
    G4double getMassOfRightAdrenal(){return volOfRightAdrenal*rhoOfRightAdrenal;};  

  void BuildBrain(const G4String&,G4bool,G4bool);
     G4double getMassOfBrain(){return volOfBrain*rhoOfBrain;};

  void BuildHeart(const G4String&,G4bool,G4bool);
     G4double getMassOfHeart(){return volOfHeart*rhoOfHeart;};

  void BuildLeftLung(const G4String&,G4bool,G4bool);
     G4double getMassOfLeftLung(){return volOfLeftLung*rhoOfLeftLung;};
  void BuildRightLung(const G4String&,G4bool,G4bool);
     G4double getMassOfRightLung(){		 
		 return volOfRightLung*rhoOfRightLung;
	 };

  void BuildStomach(const G4String&,G4bool,G4bool);
     G4double getMassOfStomach(){return volOfStomach*rhoOfStomach;};
  void BuildUpperLargeIntestine(const G4String&,G4bool,G4bool);
     G4double getMassOfUpperLargeIntestine(){return volOfUpperLargeIntestine*rhoOfUpperLargeIntestine;};
  void BuildLowerLargeIntestine(const G4String&,G4bool,G4bool);
     G4double getMassOfLowerLargeIntestine(){return volOfLowerLargeIntestine*rhoOfLowerLargeIntestine;};
 
  void BuildEsophagus(const G4String&,G4bool,G4bool);
   G4double getMassOfEsophagus(){return volOfEsophagus*rhoOfEsophagus;};
  void BuildLeftClavicle(const G4String&,G4bool,G4bool);
   G4double getMassOfLeftClavicle(){return volOfLeftClavicle*rhoOfLeftClavicle;};
  void BuildRightClavicle(const G4String&,G4bool,G4bool);
   G4double getMassOfRightClavicle(){return volOfRightClavicle*rhoOfRightClavicle;};
  void BuildGallBladder(const G4String&,G4bool,G4bool);
   G4double getMassOfGallBladder(){return volOfGallBladder*rhoOfGallBladder;};
  void BuildSmallIntestine(const G4String&,G4bool,G4bool);
   G4double getMassOfSmallIntestine(){return volOfSmallIntestine*rhoOfSmallIntestine;};
  void BuildThymus(const G4String&,G4bool,G4bool);
   G4double getMassOfThymus(){return volOfThymus*rhoOfThymus;};

  void BuildLeftKidney(const G4String&,G4bool,G4bool);
     G4double getMassOfLeftKidney(){return volOfLeftKidney*rhoOfLeftKidney;};
  void BuildRightKidney(const G4String&,G4bool,G4bool);
     G4double getMassOfRightKidney(){return volOfRightKidney*rhoOfRightKidney;};

  void BuildLiver(const G4String&,G4bool,G4bool);
     G4double getMassOfLiver(){return volOfLiver*rhoOfLiver;};
  void BuildPancreas(const G4String&,G4bool,G4bool);
     G4double getMassOfPancreas(){return volOfPancreas*rhoOfPancreas;};
  void BuildSpleen(const G4String&,G4bool,G4bool);
     G4double getMassOfSpleen(){return volOfSpleen*rhoOfSpleen;};
  void BuildUrinaryBladder(const G4String&,G4bool,G4bool);
     G4double getMassOfUrinaryBladder(){return volOfUrinaryBladder*rhoOfUrinaryBladder;};

  void BuildThyroid(const G4String&,G4bool,G4bool);
     G4double getMassOfThyroid(){return volOfThyroid*rhoOfThyroid;};

  //======================================================================
  G4double getMassOfLeftLeg(){
	  return rhoOfLeftLeg*(volOfLeftLeg-volOfLeftLegBone);
  };

  G4double getMassOfRightLeg(){
	  return rhoOfRightLeg*(volOfRightLeg-volOfRightLegBone);
  };

  G4double getMassOfHead(){return rhoOfHead*(volOfHead-
	volOfThyroid-
       -volOfSkull-volOfBrain-volOfUpperSpine
	   -volOfFacial
	   );
  };

  G4double getMassOfTrunk(){
	   return rhoOfTrunk*(volOfTrunk-
       volOfLeftArmBone-volOfRightArmBone-volOfLeftScapula-volOfRightScapula-
	   volOfMiddleLowerSpine-volOfPelvis-volOfUrinaryBladder-
	volOfHeart-volOfLeftLung-volOfRightLung-volOfLiver-
	volOfLeftClavicle-volOfRightClavicle-volOfGallBladder-
volOfEsophagus-
    volOfSmallIntestine-volOfThymus-
	   volOfStomach-volOfRibCage-volOfSpleen-volOfUpperLargeIntestine-volOfLowerLargeIntestine-
	   volOfLeftKidney-volOfRightKidney-volOfLeftAdrenal-volOfRightAdrenal-volOfPancreas
	   );
  };
  //============================================================
  void SetModel(G4String);
  void SetMotherVolume(G4VPhysicalVolume*);
 
  G4VPhysicalVolume* GetPhantom();

  void SetScaleXY(G4double val);
  void SetScaleZ(G4double val);
  void SetAgeGroup(G4String val);
protected: 

  G4double scaleXY;//scale on phantom transversal plane, x-y plane in MIRD5 description
  G4double scaleZ;//scale on phantom height, z axis in MIRD5 description.
  G4String ageGroup;

  VBodyFactory* body;

  G4String model;

  G4VPhysicalVolume* motherVolume;
  G4VPhysicalVolume* headVolume;
  G4VPhysicalVolume* trunkVolume;
  G4VPhysicalVolume* leftLegVolume;
  G4VPhysicalVolume* rightLegVolume;  
  G4VPhysicalVolume* maleGenitaliaVolume;

  G4double volOfHead;G4double rhoOfHead;
  G4double volOfSkull;G4double rhoOfSkull;
  G4double volOfFacial;G4double rhoOfFacial;
  G4double volOfBrain;G4double rhoOfBrain;
  G4double volOfTrunk;G4double rhoOfTrunk;
  G4double volOfLeftLeg;G4double rhoOfLeftLeg;
  G4double volOfRightLeg;G4double rhoOfRightLeg;
  G4double volOfLeftLegBone;G4double rhoOfLeftLegBone;
  G4double volOfRightLegBone;G4double rhoOfRightLegBone;
  G4double volOfLeftArmBone;G4double rhoOfLeftArmBone;
  G4double volOfRightArmBone;G4double rhoOfRightArmBone;
  G4double volOfLeftScapula;G4double rhoOfLeftScapula;
  G4double volOfRightScapula;G4double rhoOfRightScapula;
  G4double volOfUpperSpine;G4double rhoOfUpperSpine;
  G4double volOfMiddleLowerSpine;G4double rhoOfMiddleLowerSpine;
  G4double volOfPelvis;G4double rhoOfPelvis;
  G4double volOfHeart;G4double rhoOfHeart;//
  G4double volOfLeftLung;G4double rhoOfLeftLung;
  G4double volOfRightLung;G4double rhoOfRightLung;
  G4double volOfStomach;G4double rhoOfStomach;
  G4double volOfRibCage;G4double rhoOfRibCage;
  G4double volOfSpleen;G4double rhoOfSpleen;
  G4double volOfUpperLargeIntestine;G4double rhoOfUpperLargeIntestine;
  G4double volOfLowerLargeIntestine;G4double rhoOfLowerLargeIntestine;
  G4double volOfLeftKidney;G4double rhoOfLeftKidney;
  G4double volOfRightKidney;G4double rhoOfRightKidney;
  G4double volOfLeftAdrenal;G4double rhoOfLeftAdrenal;
  G4double volOfRightAdrenal;G4double rhoOfRightAdrenal;
  G4double volOfLiver;G4double rhoOfLiver;//
  G4double volOfPancreas;G4double rhoOfPancreas;
  G4double volOfUrinaryBladder;G4double rhoOfUrinaryBladder;
  G4double volOfThyroid;G4double rhoOfThyroid;//

  G4double volOfLeftClavicle;G4double rhoOfLeftClavicle;//
  G4double volOfRightClavicle;G4double rhoOfRightClavicle;
  G4double volOfGallBladder;G4double rhoOfGallBladder;
  G4double volOfEsophagus;G4double rhoOfEsophagus;
  G4double volOfSmallIntestine;G4double rhoOfSmallIntestine;
  G4double volOfThymus;G4double rhoOfThymus;
};
#endif
