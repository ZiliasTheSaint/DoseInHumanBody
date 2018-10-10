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
//
#ifndef HumanPhantomConstruction_H
#define HumanPhantomConstruction_H 1

#include "G4VUserDetectorConstruction.hh"

#include "HumanPhantomMessenger.hh"

#include "globals.hh"
#include <map>

class G4VPhysicalVolume;
class HumanPhantomMaterial;
class G4LogicalVolume;
//class HumanPhantomPrimaryGeneratorAction;
class HumanPhantomConstruction : public G4VUserDetectorConstruction
{
  public:
     HumanPhantomConstruction();
    ~HumanPhantomConstruction();
     G4VPhysicalVolume* Construct();

	 static G4bool isMammo;
	 const G4VPhysicalVolume* GetDetector()  {return physicalDetector;};
	 G4LogicalVolume* GetDetector_logical(){return logicDetector;};//used for computing mass and dose

	 void SetPointSourceOrFrontalBeamToDetectorDistance(G4double val){fUpperAir_height=val;};
	 G4double GetPointSourceOrFrontalBeamToDetectorDistance()           {return fUpperAir_height;};

	 void SetTotal_diameter (G4double val){ftotal_diameter=val;};
	 void SetTotal_height (G4double val){ftotal_height=val;};
	 G4double GetWorldSizeZ()           {return ftotal_height+fUpperAir_height;};
	 G4double GetWorldSizeRadius()           {return 0.5*ftotal_diameter;};

  void SetBodyPartSensitivity(G4String, G4bool);

  void SetPhantomSex(G4String);
  void SetPhantomModel(G4String);

  G4String GetPhantomSex(){return sex;};

  void SetPhantomAge(std::string val){std::istringstream buffer(val); buffer>>age;};
  double GetPhantomAge(){return age;};

  //G4VPhysicalVolume* GetMotherVolume(){return mother;};
  G4double GetBodypartMass(G4String bodypart){
	  return bodypartMass[bodypart];
  };
  double GetBodypartWeight(G4String bodypart){
	  return bodypartWeight[bodypart];
  };

  G4double GetPhantomTotalMass(){
		return phantomTotalMass;
  };
 
  G4double GetWorldSize(){
		return fWorldSize;
  };

  G4double GetScaleXY(){
	  return scaleXY;
  };
  G4double GetScaleZ(){
	  return scaleZ;
  };
  G4String GetAgeGroup(){
	  return ageGroup;
  };

  void SetPhantomHeightForScaling(G4double val){
	  phantomHeight=val;
  };
  void SetPhantomMassForScaling(G4double val){
	  phantomMass=val;
  };
  void SetAgeGroup(G4String val){
	  ageGroup=val;
  }
 private:
  HumanPhantomMaterial* material;
  HumanPhantomMessenger* messenger;

  G4VPhysicalVolume* motherWorld;
  G4VPhysicalVolume* ConstructWorld();
  G4VPhysicalVolume* BuildMaleGenitalia();
  G4VPhysicalVolume* ConstructMAMMO();
  G4VPhysicalVolume* physicalDetector; //pointer to the physical active detector
  G4LogicalVolume* logicDetector;
  G4double fUpperAir_height;
  G4double ftotal_diameter;//=14.0*cm;//cm
  G4double ftotal_height;//=5.0*cm;//cm

  G4String                 model;
  G4String                 sex;
  double				   age;
  std::map<std::string,G4bool> sensitivities;
  std::map<std::string,G4double> bodypartMass;//storing REAL mass. Constructed organs are a mess...see RIB CAGE!
  std::map<std::string,double> bodypartWeight;//storing weight for effective dose calculation!
  //===========
  G4double fWorldSize;
  G4double phantomTotalMass;
  void SetPhantomTotalMass(G4double mass){
		phantomTotalMass=mass;
  }
  //=======================
  G4double phantomHeight;//user defined
  G4double phantomMass;//user defined
  G4double actualPhantomHeight;//this phantom height = 174*cm, see LEG, TRUNK, HEAD file
  G4double actualFemalePhantomMass;//this female mass = 73.3201*kg, run a simulation to see this!
  G4double actualMalePhantomMass;//this male mass = 71.7475*kg, run a simulation to see this!
  G4double scaleXY;//scale on phantom transversal plane, x-y plane in MIRD5 description
  G4double scaleZ;//scale on phantom height, z axis in MIRD5 description.
  G4String ageGroup;//newborn, age 1,....,adult
};

#endif

