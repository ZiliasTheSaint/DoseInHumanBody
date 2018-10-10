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
// Authors: D. Fulea, National Institute of Public Health Bucharest, Cluj-Napoca Regional Center, Romania
// 

#ifndef MirdMapConstants_h
#define MirdMapConstants_h 1

#include "globals.hh"
#include <map>

struct trunk_struct{
	G4double At;
	G4double Bt;
	G4double Ct;
	G4double volume;
	G4double mass;
};
struct head_struct{
	G4double Rh;
	G4double Ah;
	G4double Bh;
	G4double Ch0;
	G4double Ch1;
	G4double Ch2;
	G4double volume;
	G4double mass;
};
struct legs_struct{
	G4double Cl;
	G4double Cpl;//C'L	
	G4double volume;
	G4double mass;
};
struct skin_struct{
	G4double S;
	G4double total_volume;
};
struct spine_struct{
	G4double a;
	G4double b;
	G4double y0;
	G4double z1;
	G4double z2;
	G4double z3;
	G4double z4;
	G4double volume;
};
struct brain_struct{
	G4double a;
	G4double b;
	G4double c;
	G4double volume;
};
struct skull_struct{
	G4double d;
	G4double volume;
};
struct armbone_struct{
	G4double a;
	G4double b;
	G4double x0;
	G4double z2;
	G4double volume;
};
struct adrenals_struct{
	G4double a;
	G4double b;
	G4double c;
	G4double x0;
	G4double y0;
	G4double z0;
	G4double theta;
	G4double volume;
};
struct breasts_struct{
	G4double a;
	G4double b;
	G4double c;
	G4double x0;	
	G4double z0;	
	G4double volume;
};
struct kidneys_struct{
	G4double a;
	G4double b;
	G4double c;
	G4double x0;	
	G4double y0;	
	G4double z0;	
	G4double x1;	
	G4double volume;
};
struct lungs_struct{
	G4double a;
	G4double b;
	G4double c;
	G4double x0;	
	G4double z0;	
	G4double x1r;	
	G4double y1r;
	G4double z1r;	
	G4double z2r;	
	G4double x1l;
	G4double y1l;	
	G4double z2l;	
	G4double volume;
};
struct ovaries_struct{
	G4double a;
	G4double b;
	G4double c;
	G4double x0;	
	G4double z0;	
	G4double volume;
};
struct testes_struct{
	G4double a;
	G4double b;
	G4double c;
	G4double y0;		
	G4double volume;
};
struct uterus_struct{
	G4double a;
	G4double b;
	G4double c;
	G4double y0;	
	G4double z0;	
	G4double y1;	
	G4double volume;
};
struct scapula_struct{
	G4double a1;
	G4double a2;
	G4double b;
	G4double z1;	
	G4double z2;	
	G4double m1;//tan!!	
	G4double m2;//tan!!	
	G4double volume;
};
struct urinaryBladder_struct{
	G4double a;
	G4double b;
	G4double c;
	G4double d;	
	G4double y0;	
	G4double z0;	
	G4double volume;
};
struct spleen_struct{
	G4double a;
	G4double b;
	G4double c;
	G4double x0;	
	G4double y0;	
	G4double z0;	
	G4double volume;
};
struct stomach_struct{
	G4double a;
	G4double b;
	G4double c;
	G4double d;
	G4double x0;	
	G4double y0;	
	G4double z0;	
	G4double volume;
};
struct pancreas_struct{
	G4double a;
	G4double b;
	G4double c;	
	G4double x0;	
	G4double z0;	
	G4double x1;	
	G4double volume;
};
struct pelvis_struct{
	G4double a1;
	G4double b1;
	G4double a2;	
	G4double b2;	
	G4double y01;	
	G4double y02;	
	G4double y1;	
	G4double z1;	
	G4double z2;	
	G4double volume;
};
struct ribCage_struct{
	G4double a;
	G4double b;
	G4double d;		
	G4double z1;	
	G4double z2;
	G4double c;
	G4double volume;
};
struct upperLargeIntestine_struct{
	G4double a;
	G4double b;
	G4double d;		
	G4double x0;	
	G4double y0;
	G4double z1;
	G4double z2;
	G4double volume;
};
struct transverseColon_struct{
	G4double b;
	G4double c;
	G4double d;		
	G4double y0;	
	G4double z0;
	G4double x1;	
	G4double volume;
};
struct lowerLargeIntestine_struct{
	G4double a;
	G4double b;
	G4double d;		
	G4double x1;	
	G4double mx;
	G4double my;
	G4double z1;
	G4double z2;
	G4double volume;
};
struct sigmoidColon_struct{
	G4double a;
	G4double b;
	G4double d;		
	G4double x0;	
	G4double z0;
	G4double R1;	
	G4double R2;	
	G4double volume;
};
struct thyroid_struct{
	G4double R;
	G4double r;
	G4double c;		
	G4double y0;	
	G4double volume;
};
struct thymus_struct{
	G4double a;
	G4double b;
	G4double c;		
	G4double y0;
	G4double z0;
	G4double volume;
};
struct liver_struct{
	G4double a;
	G4double b;
	G4double xm;		
	G4double ym;
	G4double zm;
	G4double z1;
	G4double z2;
	G4double volume;
};
struct heart_struct{
	double alpha1;
	double beta1;
	double gamma1;		
	double alpha2;
	double beta2;
	double gamma2;
	double alpha3;
	double beta3;
	double gamma3;
	G4double vx;
	G4double avy;
	G4double lavz;
	G4double ravz;
	G4double ax;
	G4double tlvw;
	G4double trvw;
	G4double tax;
	G4double x0;
	G4double y0;
	G4double z0;
	G4double volume;
};
struct smallIntestine_struct{
	G4double a;
	G4double b;
	G4double y0;		
	G4double y1;
	G4double y2;
	G4double z1;
	G4double z2;
	G4double volume;
};
struct esophagus_struct{
	G4double a;
	G4double b;
	G4double d;
	G4double y0;
	G4double z2;
	G4double z3;
	G4double r;
	G4double x1p;
	G4double x2p;
	G4double z1;	
	double alpha1;
	double beta1;
	double gamma1;		
	double alpha2;
	double beta2;
	double gamma2;
	double alpha3;
	double beta3;
	double gamma3;		
};
struct gallBladder_struct{	
	double alpha1;
	double beta1;
	double gamma1;		
	double alpha2;
	double beta2;
	double gamma2;
	double alpha3;
	double beta3;
	double gamma3;		
	G4double r1;
	G4double r2;
	double s;
	G4double h;
	G4double x0;
	G4double y0;
	G4double z0;
	G4double volume;
};
struct clavicles_struct{
	G4double y0;
	G4double z1;
	G4double R;		
	G4double r;
	double cotan1;
	double cotan2;
	G4double volume;
};
struct facial_struct{
	G4double a1;
	G4double b1;
	G4double d;		
	G4double z1;
	G4double z5;//z2
	G4double volume;
};

class MirdMapConstants
{
public:
  MirdMapConstants();
  virtual ~MirdMapConstants();
  
  void initTrunk();
  void initHead();
  void initLegs();
  void initSkin();
  void initSpine();
  void initBrain();
  void initSkull();
  void initArmBone();
  void initAdrenals();
  void initBreasts();
  void initKidneys();
  void initLungs();
  void initOvaries();
  void initTestes();
  void initUterus();
  void initScapula();
  void initUrinaryBladder();
  void initSpleen();
  void initStomach();
  void initPancreas();
  void initPelvis();
  void initRibCage();
  void initUpperLargeIntestine();//with colon
  void initLowerLargeIntestine();
  void initTransverseColon();
  void initSigmoidColon();
  void initThyroid();
  void initThymus();
  void initLiver();
  void initHeart();
  void initSmallIntestine();
  void initEsophagus();
  void initGallBladder();
  void initClavicles();
  void initFacial();

  std::map<std::string,trunk_struct> GetTrunkMap(){
	  return trunk;
  };
  std::map<std::string,head_struct> GetHeadMap(){
	  return head;
  };
  std::map<std::string,legs_struct> GetLegsMap(){
	  return legs;
  };
  std::map<std::string,skin_struct> GetSkinMap(){
	  return skin;
  };
  std::map<std::string,spine_struct> GetSpineMap(){
	  return spine;
  };
  std::map<std::string,brain_struct> GetBrainMap(){
	  return brain;
  };
  std::map<std::string,skull_struct> GetSkullMap(){
	  return skull;
  };
  std::map<std::string,armbone_struct> GetArmBoneMap(){
	  return armbone;
  };
  std::map<std::string,adrenals_struct> GetAdrenalsMap(){
	  return adrenals;
  };
  std::map<std::string,breasts_struct> GetBreastsMap(){
	  return breasts;
  };
  std::map<std::string,kidneys_struct> GetKidneysMap(){
	  return kidneys;
  };
  std::map<std::string,lungs_struct> GetLungsMap(){
	  return lungs;
  };
  std::map<std::string,ovaries_struct> GetOvariesMap(){
	  return ovaries;
  };
  std::map<std::string,testes_struct> GetTestesMap(){
	  return testes;
  };
  std::map<std::string,uterus_struct> GetUterusMap(){
	  return uterus;
  };
  std::map<std::string,scapula_struct> GetScapulaMap(){
	  return scapula;
  };
  std::map<std::string,urinaryBladder_struct> GetUrinaryBladderMap(){
	  return urinaryBladder;
  };
  std::map<std::string,spleen_struct> GetSpleenMap(){
	  return spleen;
  };
  std::map<std::string,stomach_struct> GetStomachMap(){
	  return stomach;
  };
  std::map<std::string,pancreas_struct> GetPancreasMap(){
	  return pancreas;
  };
  std::map<std::string,pelvis_struct> GetPelvisMap(){
	  return pelvis;
  };
  std::map<std::string,ribCage_struct> GetRibCageMap(){
	  return ribCage;
  };
  std::map<std::string,upperLargeIntestine_struct> GetUpperLargeIntestineMap(){
	  return upperLargeIntestine;
  };
  std::map<std::string,lowerLargeIntestine_struct> GetLowerLargeIntestineMap(){
	  return lowerLargeIntestine;
  };
  std::map<std::string,transverseColon_struct> GetTransverseColonMap(){
	  return transverseColon;
  };
  std::map<std::string,sigmoidColon_struct> GetSigmoidColonMap(){
	  return sigmoidColon;
  };
  std::map<std::string,thyroid_struct> GetThyroidMap(){
	  return thyroid;
  };
  std::map<std::string,thymus_struct> GetThymusMap(){
	  return thymus;
  };
  std::map<std::string,liver_struct> GetLiverMap(){
	  return liver;
  };
  std::map<std::string,heart_struct> GetHeartMap(){
	  return heart;
  };
  std::map<std::string,smallIntestine_struct> GetSmallIntestineMap(){
	  return smallIntestine;
  };
  std::map<std::string,esophagus_struct> GetEsophagusMap(){
	  return esophagus;
  };
  std::map<std::string,gallBladder_struct> GetGallBladderMap(){
	  return gallBladder;
  };
  std::map<std::string,clavicles_struct> GetClaviclesMap(){
	  return clavicles;
  };
  std::map<std::string,facial_struct> GetFacialMap(){
	  return facial;
  };

private:	
	std::map<std::string,trunk_struct> trunk;
	std::map<std::string,head_struct> head;
	std::map<std::string,legs_struct> legs;
	std::map<std::string,skin_struct> skin;
	std::map<std::string,spine_struct> spine;
	std::map<std::string,brain_struct> brain;
	std::map<std::string,skull_struct> skull;
	std::map<std::string,armbone_struct> armbone;
	std::map<std::string,adrenals_struct> adrenals;
	std::map<std::string,breasts_struct> breasts;
	std::map<std::string,kidneys_struct> kidneys;
	std::map<std::string,lungs_struct> lungs;
	std::map<std::string,ovaries_struct> ovaries;
	std::map<std::string,testes_struct> testes;
	std::map<std::string,uterus_struct> uterus;
	std::map<std::string,scapula_struct> scapula;
	std::map<std::string,urinaryBladder_struct> urinaryBladder;
	std::map<std::string,spleen_struct> spleen;
	std::map<std::string,stomach_struct> stomach;
	std::map<std::string,pancreas_struct> pancreas;
	std::map<std::string,pelvis_struct> pelvis;
	std::map<std::string,ribCage_struct> ribCage;
	std::map<std::string,upperLargeIntestine_struct> upperLargeIntestine;
	std::map<std::string,lowerLargeIntestine_struct> lowerLargeIntestine;
	std::map<std::string,transverseColon_struct> transverseColon;
	std::map<std::string,sigmoidColon_struct> sigmoidColon;
	std::map<std::string,thyroid_struct> thyroid;
	std::map<std::string,thymus_struct> thymus;
	std::map<std::string,liver_struct> liver;
	std::map<std::string,heart_struct> heart;
	std::map<std::string,smallIntestine_struct> smallIntestine;
	std::map<std::string,esophagus_struct> esophagus;
	std::map<std::string,gallBladder_struct> gallBladder;
	std::map<std::string,clavicles_struct> clavicles;
	std::map<std::string,facial_struct> facial;
};
#endif
