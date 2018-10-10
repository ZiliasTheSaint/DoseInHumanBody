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

#ifndef BasePhantomBuilder_h
#define BasePhantomBuilder_h 1

#include "G4VPhysicalVolume.hh"
class G4VPhysicalVolume;
class BasePhantomBuilder
{
public:

  BasePhantomBuilder();
  virtual ~BasePhantomBuilder();
 
  virtual void BuildHead(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfHead()=0;
  virtual void BuildTrunk(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfTrunk()=0;
  
  virtual void BuildUpperSpine(const G4String&,G4bool,G4bool) {return ;}virtual G4double getMassOfUpperSpine()=0;
  virtual void BuildMiddleLowerSpine(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfMiddleLowerSpine()=0;
  virtual void BuildLeftLeg(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfLeftLeg()=0;
  virtual void BuildRightLeg(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfRightLeg()=0;
  virtual void BuildLeftLegBone(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfLeftLegBone()=0;
  virtual void BuildRightLegBone(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfRightLegBone()=0;
  virtual void BuildLeftArmBone(const G4String&,G4bool,G4bool) {return ;}virtual G4double getMassOfLeftArmBone()=0;
  virtual void BuildRightArmBone(const G4String&,G4bool,G4bool) {return ;}virtual G4double getMassOfRightArmBone()=0;
  virtual void BuildSkull(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfSkull()=0;

  virtual void BuildFacial(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfFacial()=0;

  virtual void BuildRibCage(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfRibCage()=0;
  virtual void BuildPelvis(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfPelvis()=0;
 
  virtual void BuildBrain(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfBrain()=0;
  virtual void BuildHeart(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfHeart()=0;
  virtual void BuildLeftLung(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfLeftLung()=0;
  virtual void BuildRightLung(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfRightLung()=0;
  virtual void BuildStomach(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfStomach()=0;
  virtual void BuildUpperLargeIntestine(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfUpperLargeIntestine()=0;
  virtual void BuildLowerLargeIntestine(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfLowerLargeIntestine()=0;
  
  virtual void BuildEsophagus(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfEsophagus()=0;
  virtual void BuildLeftClavicle(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfLeftClavicle()=0;
  virtual void BuildRightClavicle(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfRightClavicle()=0;
  virtual void BuildLeftTesticle(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfLeftTesticle()=0;
  virtual void BuildRightTesticle(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfRightTesticle()=0;
  virtual void BuildGallBladder(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfGallBladder()=0;
  virtual void BuildSmallIntestine(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfSmallIntestine()=0;
  virtual void BuildThymus(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfThymus()=0;

  virtual void BuildLeftKidney(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfLeftKidney()=0;
  virtual void BuildRightKidney(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfRightKidney()=0;
  virtual void BuildLeftAdrenal(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfLeftAdrenal()=0;
  virtual void BuildRightAdrenal(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfRightAdrenal()=0;
  virtual void BuildLiver(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfLiver()=0;
  virtual void BuildPancreas(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfPancreas()=0;
  virtual void BuildSpleen(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfSpleen()=0;
  virtual void BuildUrinaryBladder(const G4String& ,G4bool,G4bool) {return ;};virtual G4double getMassOfUrinaryBladder()=0;
  virtual void BuildThyroid(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfThyroid()=0;
  virtual void BuildLeftScapula(const G4String&,G4bool,G4bool){return;};virtual G4double getMassOfLeftScapula()=0;
  virtual void BuildRightScapula(const G4String&,G4bool,G4bool){return;};virtual G4double getMassOfRightScapula()=0;
  
  virtual void SetModel(G4String) {return ;};
  virtual void SetMotherVolume(G4VPhysicalVolume*) {return;};
  virtual G4VPhysicalVolume* GetPhantom() {return 0;};
  virtual void SetScaleXY(G4double) {return ;};
  virtual void SetScaleZ(G4double) {return ;};
  virtual void SetAgeGroup(G4String) {return ;};

  virtual void BuildLeftOvary(const G4String&,G4bool,G4bool ) {return ;};virtual G4double getMassOfLeftOvary()=0;  
  virtual void BuildRightOvary(const G4String&,G4bool,G4bool) {return ;};virtual G4double getMassOfRightOvary()=0;
  virtual void BuildUterus(const G4String&,G4bool,G4bool){return;};virtual G4double getMassOfUterus()=0;
  virtual void BuildLeftBreast(const G4String&,G4bool,G4bool){return;};virtual G4double getMassOfLeftBreast()=0;
  virtual void BuildRightBreast(const G4String&,G4bool,G4bool){return;};virtual G4double getMassOfRightBreast()=0;
  virtual void BuildVoxelLeftBreast(const G4String&,G4bool,G4bool){return;};//
  virtual void BuildVoxelRightBreast(const G4String&,G4bool,G4bool){return;};//

};
#endif
