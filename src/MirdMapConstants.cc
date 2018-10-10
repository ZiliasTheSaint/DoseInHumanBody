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
//source: http://ordose.ornl.gov/resources/phantom.html

#include "MirdMapConstants.hh"

MirdMapConstants::MirdMapConstants()
{ }

MirdMapConstants::~MirdMapConstants()
{
	//trunk.clear();
}

void MirdMapConstants::initTrunk(){
	trunk_struct strct;
	strct.At = 6.35 *cm;
	strct.Bt = 4.90 *cm;
	strct.Ct = 21.60 *cm;
	strct.volume = 2110 *cm3;
	strct.mass = 2.100 *kg;
	trunk["newborn"] = strct;

	strct.At = 8.80 *cm;
	strct.Bt = 6.50 *cm;
	strct.Ct = 30.70 *cm;
	strct.volume = 5520 *cm3;
	strct.mass = 5.530 *kg;
	trunk["age_1"] = strct;

	strct.At = 11.45 *cm;
	strct.Bt = 7.50 *cm;
	strct.Ct = 40.80 *cm;
	strct.volume = 11000 *cm3;
	strct.mass = 11.000 *kg;
	trunk["age_5"] = strct;

	strct.At = 13.90 *cm;
	strct.Bt = 8.40 *cm;
	strct.Ct = 50.80 *cm;
	strct.volume = 18600 *cm3;
	strct.mass = 18.700 *kg;
	trunk["age_10"] = strct;

	strct.At = 17.25 *cm;
	strct.Bt = 9.80 *cm;
	strct.Ct = 63.10 *cm;
	strct.volume = 33500 *cm3;
	strct.mass = 34.500 *kg;
	trunk["age_15"] = strct;

	strct.At = 20.00 *cm;
	strct.Bt = 10.00 *cm;
	strct.Ct = 70.00 *cm;
	strct.volume = 44000 *cm3;
	strct.mass = 44.800 *kg;
	trunk["adult"] = strct;
}

void MirdMapConstants::initHead(){
	head_struct strct;
	strct.Rh = 2.8 *cm;
	strct.Ah = 4.52 *cm;
	strct.Bh = 5.78 *cm;
	strct.Ch0 = 1.56 *cm;
	strct.Ch1 = 7.01 *cm;
	strct.Ch2 = 3.99 *cm;
	strct.volume = 965 *cm3;
	strct.mass = 1.020 *kg;
	head["newborn"] = strct;

	strct.Rh = 3.6 *cm;
	strct.Ah = 6.13 *cm;
	strct.Bh = 7.84 *cm;
	strct.Ch0 = 2.30 *cm;
	strct.Ch1 = 9.50 *cm;
	strct.Ch2 = 5.41 *cm;
	strct.volume = 2410 *cm3;
	strct.mass = 2.580 *kg;
	head["age_1"] = strct;
		
	strct.Rh = 3.8 *cm;
	strct.Ah = 7.13 *cm;
	strct.Bh = 9.05 *cm;
	strct.Ch0 = 3.30 *cm;
	strct.Ch1 = 10.70 *cm;
	strct.Ch2 = 6.31 *cm;
	strct.volume = 3670 *cm3;
	strct.mass = 4.000 *kg;
	head["age_5"] = strct;
		
	strct.Rh = 4.4 *cm;
	strct.Ah = 7.43 *cm;
	strct.Bh = 9.40 *cm;
	strct.Ch0 = 4.70 *cm;
	strct.Ch1 = 11.68 *cm;
	strct.Ch2 = 6.59 *cm;
	strct.volume = 4300 *cm3;
	strct.mass = 4.710 *kg;
	head["age_10"] = strct;
		
	strct.Rh = 5.2 *cm;
	strct.Ah = 7.77 *cm;
	strct.Bh = 9.76 *cm;
	strct.Ch0 = 7.70 *cm;
	strct.Ch1 = 12.35 *cm;
	strct.Ch2 = 6.92 *cm;
	strct.volume = 4900 *cm3;
	strct.mass = 5.410 *kg;
	head["age_15"] = strct;
		
	strct.Rh = 5.4 *cm;
	strct.Ah = 8.00 *cm;
	strct.Bh = 10.00 *cm;
	strct.Ch0 = 8.40 *cm;
	strct.Ch1 = 13.05 *cm;
	strct.Ch2 = 7.15 *cm;
	strct.volume = 5430 *cm3;
	strct.mass = 6.040 *kg;
	head["adult"] = strct;
}

void MirdMapConstants::initLegs(){
	legs_struct strct;
	strct.Cl = 16.8 *cm;
	strct.Cpl = 21.6 *cm;	
	strct.volume = 451 *cm3;
	strct.mass = 0.48 *kg;
	legs["newborn"] = strct;

	strct.Cl = 26.5 *cm;
	strct.Cpl = 37.1*cm;	
	strct.volume = 1470 *cm3;
	strct.mass = 1.6 *kg;
	legs["age_1"] = strct;

	strct.Cl = 48.0 *cm;
	strct.Cpl = 65.0 *cm;	
	strct.volume = 4380 *cm3;
	strct.mass = 4.78 *kg;
	legs["age_5"] = strct;

	strct.Cl = 66.0 *cm;
	strct.Cpl = 90.0 *cm;	
	strct.volume = 8930 *cm3;
	strct.mass = 9.74 *kg;
	legs["age_10"] = strct;

	strct.Cl = 78.0 *cm;
	strct.Cpl = 100.0 *cm;	
	strct.volume = 15400 *cm3;
	strct.mass = 16.800 *kg;
	legs["age_15"] = strct;

	strct.Cl = 80.0 *cm;
	strct.Cpl = 100.0 *cm;	
	strct.volume = 20800 *cm3;
	strct.mass = 22.660 *kg;
	legs["adult"] = strct;
}

void MirdMapConstants::initSkin(){
	skin_struct strct;
	strct.S = 0.07 *cm;	
	strct.total_volume=114 *cm3;
	skin["newborn"] = strct;

	strct.S = 0.08 *cm;
	strct.total_volume=261 *cm3;
	skin["age_1"] = strct;

	strct.S = 0.09 *cm;
	strct.total_volume=517 *cm3;
	skin["age_5"] = strct;

	strct.S = 0.10 *cm;
	strct.total_volume=854 *cm3;
	skin["age_10"] = strct;

	strct.S = 0.17 *cm;
	strct.total_volume=2050 *cm3;
	skin["age_15"] = strct;

	strct.S = 0.20 *cm;
	strct.total_volume=2890 *cm3;
	skin["adult"] = strct;
}

void MirdMapConstants::initSpine(){
	spine_struct strct;
	strct.a = 0.64 *cm;		
	strct.b = 1.23 *cm;		
	strct.y0 = 2.70 *cm;		
	strct.z1 = 6.79 *cm;		
	strct.z2 = 10.83 *cm;		
	strct.z3 = 21.60 *cm;		
	strct.z4 = 26.52 *cm;	//27.02 *cm;	///must adjust
	strct.volume=50 *cm3;
	spine["newborn"] = strct;

	strct.a = 0.88 *cm;		
	strct.b = 1.63 *cm;		
	strct.y0 = 3.58 *cm;		
	strct.z1 = 9.65 *cm;		
	strct.z2 = 15.39 *cm;		
	strct.z3 = 30.70 *cm;		
	strct.z4 = 37.51 *cm;//38.01 *cm;	///must adjust
	strct.volume=128 *cm3;
	spine["age_1"] = strct;

	strct.a = 1.15 *cm;		
	strct.b = 1.88 *cm;		
	strct.y0 = 4.13 *cm;		
	strct.z1 = 12.82 *cm;		
	strct.z2 = 20.46 *cm;		
	strct.z3 = 40.80 *cm;		
	strct.z4 = 48.83 *cm;	///must adjust..ok!
	strct.volume=245 *cm3;
	spine["age_5"] = strct;

	strct.a = 1.39 *cm;		
	strct.b = 2.10 *cm;		
	strct.y0 = 4.62 *cm;		
	strct.z1 = 15.97 *cm;		
	strct.z2 = 25.47 *cm;		
	strct.z3 = 50.80 *cm;		
	strct.z4 = 60.87 *cm;//59.89 *cm;	///must adjust
	strct.volume=403 *cm3;
	spine["age_10"] = strct;

	strct.a = 1.73 *cm;		
	strct.b = 2.45 *cm;		
	strct.y0 = 5.39 *cm;		
	strct.z1 = 19.83 *cm;		
	strct.z2 = 31.64 *cm;		
	strct.z3 = 63.10 *cm;		
	strct.z4 = 76.9*cm;//72.91 *cm;	///must adjust
	strct.volume=707 *cm3;
	spine["age_15"] = strct;

	strct.a = 2.00 *cm;		
	strct.b = 2.50 *cm;		
	strct.y0 = 5.50 *cm;		
	strct.z1 = 22.00 *cm;		
	strct.z2 = 35.10 *cm;		
	strct.z3 = 70.00 *cm;		
	strct.z4 = 85.0*cm;//80.54 *cm;	///must adjust
	strct.volume=920 *cm3;
	spine["adult"] = strct;
}

void MirdMapConstants::initBrain(){
	brain_struct strct;
	strct.a = 4.14 *cm;		
	strct.b = 5.40 *cm;		
	strct.c = 3.61 *cm;
	strct.volume=338 *cm3;
	brain["newborn"] = strct;

	strct.a = 5.63 *cm;		
	strct.b = 7.34 *cm;		
	strct.c = 4.91 *cm;		
	strct.volume=850 *cm3;
	brain["age_1"] = strct;

	strct.a = 6.34 *cm;		
	strct.b = 8.26 *cm;		
	strct.c = 5.52 *cm;		
	strct.volume=1210 *cm3;
	brain["age_5"] = strct;

	strct.a = 6.51 *cm;		
	strct.b = 8.48 *cm;		
	strct.c = 5.67 *cm;	
	strct.volume=1310 *cm3;
	brain["age_10"] = strct;

	strct.a = 6.58 *cm;		
	strct.b = 8.57 *cm;		
	strct.c = 5.73 *cm;		
	strct.volume=1350 *cm3;
	brain["age_15"] = strct;

	strct.a = 6.60 *cm;		
	strct.b = 8.60 *cm;		
	strct.c = 5.75 *cm;		
	strct.volume=1370 *cm3;
	brain["adult"] = strct;
}

void MirdMapConstants::initSkull(){
	skull_struct strct;
	strct.d = 0.20 *cm;		
	strct.volume=49.8 *cm3;
	skull["newborn"] = strct;

	strct.d = 0.30 *cm;		
	strct.volume=139 *cm3;
	skull["age_1"] = strct;

	strct.d = 0.56 *cm;		
	strct.volume=339 *cm3;
	skull["age_5"] = strct;

	strct.d = 0.67*cm;		
	strct.volume=434 *cm3;
	skull["age_10"] = strct;

	strct.d = 0.76 *cm;		
	strct.volume=508 *cm3;
	skull["age_15"] = strct;

	strct.d = 0.90 *cm;		
	strct.volume=618 *cm3;
	skull["adult"] = strct;
}
void MirdMapConstants::initArmBone(){
	armbone_struct strct;
	strct.a = 0.44 *cm;		
	strct.b = 1.32 *cm;		
	strct.x0 = 5.84 *cm;//+/-		
	strct.z2 = 21.29 *cm;		
	strct.volume=45.3 *cm3;
	armbone["newborn"] = strct;

	strct.a = 0.62 *cm;		
	strct.b = 1.76 *cm;		
	strct.x0 = 8.10 *cm;//+/-		
	strct.z2 = 30.26 *cm;		
	strct.volume=121 *cm3;
	armbone["age_1"] = strct;

	strct.a = 0.80 *cm;		
	strct.b = 2.03 *cm;		
	strct.x0 = 10.53 *cm;//+/-		
	strct.z2 = 40.22 *cm;		
	strct.volume=239 *cm3;
	armbone["age_5"] = strct;

	strct.a = 0.97 *cm;		
	strct.b = 2.27 *cm;		
	strct.x0 = 12.79 *cm;//+/-		
	strct.z2 = 50.07 *cm;		
	strct.volume=404 *cm3;
	armbone["age_10"] = strct;

	strct.a = 1.21 *cm;		
	strct.b = 2.65 *cm;		
	strct.x0 = 15.87 *cm;//+/-		
	strct.z2 = 62.20 *cm;		
	strct.volume=731 *cm3;
	armbone["age_15"] = strct;

	strct.a = 1.40 *cm;		
	strct.b = 2.70 *cm;		
	strct.x0 = 18.40 *cm;//+/-		
	strct.z2 = 69.00 *cm;		
	strct.volume=956 *cm3;
	armbone["adult"] = strct;
}
void MirdMapConstants::initAdrenals(){
	adrenals_struct strct;
	strct.a = 1.61 *cm;		
	strct.b = 0.54 *cm;		
	strct.c = 1.54 *cm;		
	strct.x0 = 1.41 *cm;//+/-
	strct.y0 = 2.45 *cm;
	strct.z0 = 11.73 *cm;		
	strct.theta = 63.3 *deg;//+/-		
	strct.volume=5.61 *cm3;
	adrenals["newborn"] = strct;

	strct.a = 1.05 *cm;		
	strct.b = 0.35 *cm;		
	strct.c = 2.20 *cm;		
	strct.x0 = 1.54 *cm;//+/-
	strct.y0 = 3.25 *cm;
	strct.z0 = 16.66 *cm;		
	strct.theta = 62.2 *deg;//+/-		
	strct.volume=3.39 *cm3;
	adrenals["age_1"] = strct;

	strct.a = 1.12 *cm;		
	strct.b = 0.37 *cm;		
	strct.c = 2.92 *cm;		
	strct.x0 = 2.00 *cm;//+/-
	strct.y0 = 3.75 *cm;
	strct.z0 = 22.14 *cm;		
	strct.theta = 59.3 *deg;//+/-		
	strct.volume=5.07 *cm3;
	adrenals["age_5"] = strct;

	strct.a = 1.17 *cm;		
	strct.b = 0.39 *cm;		
	strct.c = 3.63 *cm;		
	strct.x0 = 2.43 *cm;//+/-
	strct.y0 = 4.20 *cm;
	strct.z0 = 27.58 *cm;		
	strct.theta = 57.2 *deg;//+/-		
	strct.volume=6.94 *cm3;
	adrenals["age_10"] = strct;

	strct.a = 1.30 *cm;		
	strct.b = 0.43 *cm;		
	strct.c = 4.30 *cm;		
	strct.x0 = 3.02 *cm;//+/-
	strct.y0 = 4.90 *cm;
	strct.z0 = 34.26 *cm;		
	strct.theta = 55.6 *deg;//+/-		
	strct.volume=10.01 *cm3;
	adrenals["age_15"] = strct;

	strct.a = 1.50 *cm;		
	strct.b = 0.50 *cm;		
	strct.c = 5.00 *cm;		
	strct.x0 = 3.50 *cm;//+/-
	strct.y0 = 5.00 *cm;
	strct.z0 = 38.00 *cm;		
	strct.theta = 52.0 *deg;//+/-		
	strct.volume=15.7 *cm3;
	adrenals["adult"] = strct;
}
void MirdMapConstants::initBreasts(){
	breasts_struct strct;
	strct.a = 0.36 *cm;		
	strct.b = 0.36 *cm;		
	strct.c = 0.36 *cm;		
	strct.x0 = 3.18 *cm;//+/-
	strct.z0 = 16.05 *cm;		
	strct.volume=0.197 *cm3;
	breasts["newborn"] = strct;

	strct.a = 0.63 *cm;		
	strct.b = 0.63 *cm;		
	strct.c = 0.63 *cm;		
	strct.x0 = 4.40 *cm;//+/-
	strct.z0 = 22.81 *cm;		
	strct.volume=1.06 *cm3;
	breasts["age_1"] = strct;

	strct.a = 0.79 *cm;		
	strct.b = 0.79 *cm;		
	strct.c = 0.79 *cm;		
	strct.x0 = 5.73 *cm;//+/-
	strct.z0 = 30.31 *cm;		
	strct.volume=2.09 *cm3;
	breasts["age_5"] = strct;

	strct.a = 0.94 *cm;		
	strct.b = 0.94 *cm;		
	strct.c = 0.94 *cm;		
	strct.x0 = 6.95 *cm;//+/-
	strct.z0 = 37.73 *cm;		
	strct.volume=3.51 *cm3;
	breasts["age_10"] = strct;

	strct.a = 4.95 *cm;		
	strct.b = 4.35 *cm;		
	strct.c = 4.15 *cm;		
	strct.x0 = 8.63 *cm;//+/-
	strct.z0 = 46.87 *cm;		
	strct.volume=391 *cm3;
	breasts["age_15"] = strct;

	strct.a = 4.95 *cm;		
	strct.b = 4.35 *cm;		
	strct.c = 4.15 *cm;		
	strct.x0 = 10.00 *cm;//+/-
	strct.z0 = 52.00 *cm;		
	strct.volume=388 *cm3;
	breasts["adult"] = strct;
}
void MirdMapConstants::initKidneys(){
	kidneys_struct strct;
	strct.a = 1.79 *cm;		
	strct.b = 0.93 *cm;		
	strct.c = 1.70 *cm;		
	strct.x0 = 1.91 *cm;//+/-
	strct.y0 = 2.94 *cm;
	strct.z0 = 10.03 *cm;		
	strct.x1 = 0.71 *cm;
	strct.volume=22.0 *cm3;
	kidneys["newborn"] = strct;

	strct.a = 2.61 *cm;		
	strct.b = 1.25 *cm;		
	strct.c = 2.41 *cm;		
	strct.x0 = 2.64 *cm;//+/-
	strct.y0 = 3.90 *cm;
	strct.z0 = 14.25 *cm;		
	strct.x1 = 0.95 *cm;
	strct.volume=60.5 *cm3;
	kidneys["age_1"] = strct;

	strct.a = 3.20 *cm;		
	strct.b = 1.40 *cm;		
	strct.c = 3.20 *cm;		
	strct.x0 = 3.44 *cm;//+/-
	strct.y0 = 4.50 *cm;
	strct.z0 = 18.94 *cm;		
	strct.x1 = 1.31 *cm;
	strct.volume=111 *cm3;
	kidneys["age_5"] = strct;

	strct.a = 3.66 *cm;		
	strct.b = 1.47 *cm;		
	strct.c = 3.99 *cm;		
	strct.x0 = 4.17 *cm;//+/-
	strct.y0 = 5.04 *cm;
	strct.z0 = 23.59 *cm;		
	strct.x1 = 1.74 *cm;
	strct.volume=166 *cm3;
	kidneys["age_10"] = strct;

	strct.a = 4.05 *cm;		
	strct.b = 1.53 *cm;		
	strct.c = 4.96 *cm;		
	strct.x0 = 5.18 *cm;//+/-
	strct.y0 = 5.88 *cm;
	strct.z0 = 29.30 *cm;		
	strct.x1 = 2.48 *cm;
	strct.volume=238 *cm3;
	kidneys["age_15"] = strct;

	strct.a = 4.50 *cm;		
	strct.b = 1.50 *cm;		
	strct.c = 5.50 *cm;		
	strct.x0 = 6.00 *cm;//+/-
	strct.y0 = 6.00 *cm;
	strct.z0 = 32.50 *cm;		
	strct.x1 = 3.00 *cm;
	strct.volume=288 *cm3;
	kidneys["adult"] = strct;
}
void MirdMapConstants::initLungs(){
	lungs_struct strct;
	strct.a = 1.89 *cm;		
	strct.b = 3.68 *cm;		
	strct.c = 7.41 *cm;		
	strct.x0 = 2.70 *cm;//+/-
	strct.z0 = 13.42 *cm;		
	strct.x1r = -2.30 *cm;
	strct.y1r = 0.75 *cm;
	strct.z1r = 14.15 *cm;
	strct.z2r = 17.85 *cm;
	strct.x1l = 3.00 *cm;
	strct.y1l = 0.30 *cm;
	strct.z2l = 17.90 *cm;
	strct.volume=171 *cm3;
	lungs["newborn"] = strct;

	strct.a = 2.68 *cm;		
	strct.b = 4.88 *cm;		
	strct.c = 10.53 *cm;		
	strct.x0 = 3.74 *cm;//+/-
	strct.z0 = 19.08 *cm;		
	strct.x1r = -2.90 *cm;
	strct.y1r = 0.70 *cm;
	strct.z1r = 20.10 *cm;
	strct.z2r = 24.60 *cm;
	strct.x1l = 3.90 *cm;
	strct.y1l = 0.40 *cm;
	strct.z2l = 24.80 *cm;
	strct.volume=484 *cm3;
	lungs["age_1"] = strct;

	strct.a = 3.47 *cm;		
	strct.b = 5.63 *cm;		
	strct.c = 13.99 *cm;		
	strct.x0 = 4.87 *cm;//+/-
	strct.z0 = 25.35 *cm;		
	strct.x1r = -3.50 *cm;
	strct.y1r = 1.00 *cm;
	strct.z1r = 26.90 *cm;
	strct.z2r = 32.30 *cm;
	strct.x1l = 5.00 *cm;
	strct.y1l = 0.50 *cm;
	strct.z2l = 32.60 *cm;
	strct.volume=980 *cm3;
	lungs["age_5"] = strct;

	strct.a = 3.82 *cm;		
	strct.b = 6.30 *cm;		
	strct.c = 17.42 *cm;		
	strct.x0 = 5.91 *cm;//+/-
	strct.z0 = 31.57 *cm;		
	strct.x1r = -4.10 *cm;
	strct.y1r = 1.30 *cm;
	strct.z1r = 33.40 *cm;
	strct.z2r = 39.60 *cm;
	strct.x1l = 5.90 *cm;
	strct.y1l = 0.75 *cm;
	strct.z2l = 40.00 *cm;
	strct.volume=1530 *cm3;
	lungs["age_10"] = strct;

	strct.a = 4.09 *cm;		
	strct.b = 6.98 *cm;		
	strct.c = 20.55 *cm;		
	strct.x0 = 7.33 *cm;//+/-
	strct.z0 = 39.21 *cm;		
	strct.x1r = -5.00 *cm;
	strct.y1r = 1.20 *cm;
	strct.z1r = 41.60 *cm;
	strct.z2r = 48.50 *cm;
	strct.x1l = 7.00 *cm;
	strct.y1l = 0.70 *cm;
	strct.z2l = 49.00 *cm;
	strct.volume=2200 *cm3;
	lungs["age_15"] = strct;

	strct.a = 5.00 *cm;		
	strct.b = 7.50 *cm;		
	strct.c = 24.00 *cm;		
	strct.x0 = 8.50 *cm;//+/-
	strct.z0 = 43.50 *cm;		
	strct.x1r = -5.40 *cm;
	strct.y1r = 1.50 *cm;
	strct.z1r = 46.00 *cm;
	strct.z2r = 54.00 *cm;
	strct.x1l = 8.00 *cm;
	strct.y1l = 1.00 *cm;
	strct.z2l = 55.00 *cm;
	strct.volume=3380 *cm3;
	lungs["adult"] = strct;
}
void MirdMapConstants::initOvaries(){
	ovaries_struct strct;
	strct.a = 0.30 *cm;		
	strct.b = 0.22 *cm;		
	strct.c = 0.57 *cm;		
	strct.x0 = 1.91 *cm;//+/-
	strct.z0 = 4.63 *cm;		
	strct.volume=0.315 *cm3;
	ovaries["newborn"] = strct;

	strct.a = 0.38 *cm;		
	strct.b = 0.28 *cm;		
	strct.c = 0.77 *cm;		
	strct.x0 = 2.64 *cm;//+/-
	strct.z0 = 6.58 *cm;		
	strct.volume=0.686 *cm3;
	ovaries["age_1"] = strct;

	strct.a = 0.53 *cm;		
	strct.b = 0.35 *cm;		
	strct.c = 1.07 *cm;		
	strct.x0 = 3.44 *cm;//+/-
	strct.z0 = 8.74 *cm;		
	strct.volume=1.66 *cm3;
	ovaries["age_5"] = strct;

	strct.a = 0.66 *cm;		
	strct.b = 0.40 *cm;		
	strct.c = 1.36 *cm;		
	strct.x0 = 4.17 *cm;//+/-
	strct.z0 = 10.89 *cm;		
	strct.volume=3.01 *cm3;
	ovaries["age_10"] = strct;

	strct.a = 1.10*cm;//1.17 *cm;		
	strct.b = 0.50*cm;//0.58 *cm;		
	strct.c = 1.80 *cm;//1.80 *cm		
	strct.x0 = 5.18 *cm;//+/-
	strct.z0 = 13.52 *cm;		
	strct.volume=10.01 *cm3;//larger than adult! ok!
	ovaries["age_15"] = strct;

	strct.a = 1.00 *cm;		
	strct.b = 0.50 *cm;		
	strct.c = 2.00 *cm;		
	strct.x0 = 6.00 *cm;//+/-
	strct.z0 = 15.00 *cm;		
	strct.volume=8.38 *cm3;
	ovaries["adult"] = strct;
}
void MirdMapConstants::initTestes(){
	testes_struct strct;
	strct.a = 0.36 *cm;		
	strct.b = 0.42 *cm;		
	strct.c = 0.64 *cm;		
	strct.y0 = -2.58 *cm;	
	strct.volume=0.811 *cm3;
	testes["newborn"] = strct;

	strct.a = 0.41 *cm;		
	strct.b = 0.47 *cm;		
	strct.c = 0.72 *cm;		
	strct.y0 = -3.73 *cm;			
	strct.volume=1.16 *cm3;
	testes["age_1"] = strct;

	strct.a = 0.45 *cm;		
	strct.b = 0.52 *cm;		
	strct.c = 0.80 *cm;		
	strct.y0 = -4.98 *cm;			
	strct.volume=1.57 *cm3;
	testes["age_5"] = strct;

	strct.a = 0.47 *cm;		
	strct.b = 0.55 *cm;		
	strct.c = 0.84 *cm;		
	strct.y0 = -6.15 *cm;	
	strct.volume=1.82 *cm3;
	testes["age_10"] = strct;

	strct.a = 0.96 *cm;		
	strct.b = 1.10 *cm;		
	strct.c = 1.69 *cm;		
	strct.y0 = -7.10 *cm;			
	strct.volume=15 *cm3;
	testes["age_15"] = strct;

	strct.a = 1.30 *cm;		
	strct.b = 1.50 *cm;		
	strct.c = 2.30 *cm;		
	strct.y0 = -8.00 *cm;			
	strct.volume=37.6 *cm3;
	testes["adult"] = strct;
}
void MirdMapConstants::initUterus(){
	uterus_struct strct;
	strct.a = 0.83 *cm;		
	strct.b = 2.57 *cm;		
	strct.c = 0.49 *cm;		
	strct.y0 = -0.98 *cm;	
	strct.z0 = 4.32 *cm;	
	strct.y1 = -2.27 *cm;	
	strct.volume=3.70 *cm3;
	uterus["newborn"] = strct;

	strct.a = 0.61 *cm;		
	strct.b = 1.80 *cm;		
	strct.c = 0.36 *cm;		
	strct.y0 = -1.30 *cm;	
	strct.z0 = 6.14 *cm;	
	strct.y1 = -2.20 *cm;	
	strct.volume=1.40 *cm3;
	uterus["age_1"] = strct;

	strct.a = 0.78 *cm;		
	strct.b = 2.00 *cm;		
	strct.c = 0.47 *cm;		
	strct.y0 = -1.50 *cm;	
	strct.z0 = 8.16 *cm;	
	strct.y1 = -2.51 *cm;	
	strct.volume=2.60 *cm3;
	uterus["age_5"] = strct;

	strct.a = 0.91 *cm;		
	strct.b = 2.17 *cm;		
	strct.c = 0.57 *cm;		
	strct.y0 = -1.68 *cm;	
	strct.z0 = 10.16 *cm;	
	strct.y1 = -2.78 *cm;	
	strct.volume=4.00 *cm3;
	uterus["age_10"] = strct;

	strct.a = 2.47 *cm;		
	strct.b = 5.61 *cm;		
	strct.c = 1.55 *cm;		
	strct.y0 = -1.96 *cm;	
	strct.z0 = 12.62 *cm;	
	strct.y1 = -4.77 *cm;	
	strct.volume=76 *cm3;
	uterus["age_15"] = strct;

	strct.a = 2.62 *cm;		
	strct.b = 5.22 *cm;		
	strct.c = 1.57 *cm;		
	strct.y0 = -2.00 *cm;	
	strct.z0 = 14.00 *cm;	
	strct.y1 = -4.62 *cm;	
	strct.volume=76 *cm3;
	uterus["adult"] = strct;
}
void MirdMapConstants::initScapula(){
	scapula_struct strct;
	strct.a1 = 5.40 *cm;		
	strct.a2 = 6.04 *cm;		
	strct.b = 4.80 *cm;		
	strct.z1 = 15.71 *cm;
	strct.z2 = 20.77 *cm;
	strct.m1 = 0.39;		
	strct.m2 = 1.23;
	strct.volume=9.64 *cm3;
	scapula["newborn"] = strct;

	strct.a1 = 7.48 *cm;		
	strct.a2 = 8.36 *cm;		
	strct.b = 6.37 *cm;		
	strct.z1 = 22.32 *cm;
	strct.z2 = 29.52 *cm;
	strct.m1 = 0.37;		
	strct.m2 = 1.18;
	strct.volume=25.3 *cm3;
	scapula["age_1"] = strct;

	strct.a1 = 9.73 *cm;		
	strct.a2 = 10.88 *cm;		
	strct.b = 7.35 *cm;		
	strct.z1 = 29.67 *cm;
	strct.z2 = 39.23 *cm;
	strct.m1 = 0.33;		
	strct.m2 = 1.05;
	strct.volume=50.4 *cm3;
	scapula["age_5"] = strct;

	strct.a1 = 11.82 *cm;		
	strct.a2 = 13.20 *cm;		
	strct.b = 8.23 *cm;		
	strct.z1 = 36.94 *cm;
	strct.z2 = 48.84 *cm;
	strct.m1 = 0.30;		
	strct.m2 = 0.97;
	strct.volume=85.7 *cm3;
	scapula["age_10"] = strct;

	strct.a1 = 14.66 *cm;		
	strct.a2 = 16.36 *cm;		
	strct.b = 9.60 *cm;		
	strct.z1 = 45.88 *cm;
	strct.z2 = 60.67 *cm;
	strct.m1 = 0.28;		
	strct.m2 = 0.91;
	strct.volume=154 *cm3;
	scapula["age_15"] = strct;

	strct.a1 = 17.00 *cm;		
	strct.a2 = 19.00 *cm;		
	strct.b = 9.80 *cm;		
	strct.z1 = 50.90 *cm;
	strct.z2 = 67.30 *cm;
	strct.m1 = 0.25;		
	strct.m2 = 0.80;
	strct.volume=202 *cm3;
	scapula["adult"] = strct;
}
void MirdMapConstants::initUrinaryBladder(){
	urinaryBladder_struct strct;
	strct.a = 1.69 *cm;		
	strct.b = 1.82 *cm;		
	strct.c = 1.14 *cm;		
	strct.d = 0.10 *cm;
	strct.y0 = -2.21 *cm;
	strct.z0 = 2.47 *cm;
	strct.volume=14.67 *cm3;
	urinaryBladder["newborn"] = strct;

	strct.a = 2.35 *cm;		
	strct.b = 2.42 *cm;		
	strct.c = 1.64 *cm;		
	strct.d = 0.14 *cm;
	strct.y0 = -2.93 *cm;
	strct.z0 = 3.51 *cm;
	strct.volume=39.11 *cm3;
	urinaryBladder["age_1"] = strct;

	strct.a = 3.04 *cm;		
	strct.b = 2.77 *cm;		
	strct.c = 2.16 *cm;		
	strct.d = 0.17 *cm;
	strct.y0 = -3.38 *cm;
	strct.z0 = 4.66 *cm;
	strct.volume=76.2 *cm3;
	urinaryBladder["age_5"] = strct;

	strct.a = 3.61 *cm;		
	strct.b = 3.04 *cm;		
	strct.c = 2.63 *cm;		
	strct.d = 0.20 *cm;
	strct.y0 = -3.78 *cm;
	strct.z0 = 5.81 *cm;
	strct.volume=120.9 *cm3;
	urinaryBladder["age_10"] = strct;

	strct.a = 4.27 *cm;		
	strct.b = 3.38 *cm;		
	strct.c = 3.11 *cm;		
	strct.d = 0.23 *cm;
	strct.y0 = -4.41 *cm;
	strct.z0 = 7.21 *cm;
	strct.volume=188.5 *cm3;
	urinaryBladder["age_15"] = strct;

	strct.a = 4.958 *cm;		
	strct.b = 3.458 *cm;		
	strct.c = 3.458 *cm;		
	strct.d = 0.252 *cm;
	strct.y0 = -4.50 *cm;
	strct.z0 = 8.00 *cm;
	strct.volume=245.7 *cm3;
	urinaryBladder["adult"] = strct;
}
void MirdMapConstants::initSpleen(){
	spleen_struct strct;
	strct.a = 1.13 *cm;		
	strct.b = 1.00 *cm;		
	strct.c = 1.85 *cm;		
	strct.x0 = 3.54 *cm;
	strct.y0 = 1.42 *cm;
	strct.z0 = 11.42 *cm;
	strct.volume=8.76 *cm3;
	spleen["newborn"] = strct;

	strct.a = 1.65 *cm;		
	strct.b = 1.35 *cm;		
	strct.c = 2.63 *cm;		
	strct.x0 = 4.94 *cm;
	strct.y0 = 1.85 *cm;
	strct.z0 = 16.23 *cm;
	strct.volume=24.05 *cm3;
	spleen["age_1"] = strct;

	strct.a = 2.09 *cm;		
	strct.b = 1.52 *cm;		
	strct.c = 3.49 *cm;		
	strct.x0 = 6.40 *cm;
	strct.y0 = 2.25 *cm;
	strct.z0 = 21.57 *cm;
	strct.volume=46.4 *cm3;
	spleen["age_5"] = strct;

	strct.a = 2.43 *cm;		
	strct.b = 1.68 *cm;		
	strct.c = 4.35 *cm;		
	strct.x0 = 7.65 *cm;
	strct.y0 = 2.52 *cm;
	strct.z0 = 26.85 *cm;
	strct.volume=74.4 *cm3;
	spleen["age_10"] = strct;

	strct.a = 2.90 *cm;		
	strct.b = 1.88 *cm;		
	strct.c = 5.19 *cm;		
	strct.x0 = 9.49 *cm;
	strct.y0 = 2.94 *cm;
	strct.z0 = 33.35 *cm;
	strct.volume=119 *cm3;
	spleen["age_15"] = strct;

	strct.a = 3.50 *cm;		
	strct.b = 2.00 *cm;		
	strct.c = 6.00 *cm;		
	strct.x0 = 11.00 *cm;
	strct.y0 = 3.00 *cm;
	strct.z0 = 37.00 *cm;
	strct.volume=176 *cm3;
	spleen["adult"] = strct;
}
void MirdMapConstants::initStomach(){
	stomach_struct strct;
	strct.a = 1.20 *cm;		
	strct.b = 1.39 *cm;		
	strct.c = 2.34 *cm;		
	strct.d = 0.22 *cm;		
	strct.x0 = 2.54 *cm;
	strct.y0 = -1.96 *cm;
	strct.z0 = 10.80 *cm;
	strct.volume=16.37 *cm3;
	stomach["newborn"] = strct;

	strct.a = 1.85 *cm;		
	strct.b = 2.05 *cm;		
	strct.c = 3.51 *cm;	
	strct.d = 0.33 *cm;
	strct.x0 = 3.52 *cm;
	strct.y0 = -2.70 *cm;
	strct.z0 = 15.35 *cm;
	strct.volume=55.7 *cm3;
	stomach["age_1"] = strct;

	strct.a = 2.55 *cm;		
	strct.b = 2.40 *cm;		
	strct.c = 4.66 *cm;	
	strct.d = 0.45 *cm;
	strct.x0 = 4.58 *cm;
	strct.y0 = -3.15 *cm;
	strct.z0 = 20.40 *cm;
	strct.volume=119.4 *cm3;
	stomach["age_5"] = strct;

	strct.a = 3.14 *cm;		
	strct.b = 2.74 *cm;		
	strct.c = 5.81 *cm;
	strct.d = 0.53 *cm;
	strct.x0 = 5.56 *cm;
	strct.y0 = -3.51 *cm;
	strct.z0 = 25.40 *cm;
	strct.volume=209.8 *cm3;
	stomach["age_10"] = strct;

	strct.a = 3.43 *cm;		
	strct.b = 2.92 *cm;		
	strct.c = 7.16 *cm;	
	strct.d = 0.56 *cm;
	strct.x0 = 6.90 *cm;
	strct.y0 = -3.92 *cm;
	strct.z0 = 31.55 *cm;
	strct.volume=300 *cm3;
	stomach["age_15"] = strct;

	strct.a = 4.00 *cm;		
	strct.b = 3.00 *cm;		
	strct.c = 8.00 *cm;
	strct.d = 0.613 *cm;
	strct.x0 = 8.00 *cm;
	strct.y0 = -4.00 *cm;
	strct.z0 = 35.00 *cm;
	strct.volume=402 *cm3;
	stomach["adult"] = strct;
}
void MirdMapConstants::initPancreas(){
	pancreas_struct strct;
	strct.a = 4.32 *cm;		
	strct.b = 0.50 *cm;		
	strct.c = 0.87 *cm;			
	strct.x0 = -0.09 *cm;
	strct.z0 = 11.42 *cm;
	strct.x1 = 0.99 *cm;
	strct.volume=2.69 *cm3;
	pancreas["newborn"] = strct;

	strct.a = 6.85 *cm;		
	strct.b = 0.71 *cm;		
	strct.c = 1.41 *cm;			
	strct.x0 = -0.43 *cm;
	strct.z0 = 16.23 *cm;
	strct.x1 = 1.32 *cm;
	strct.volume=9.87 *cm3;
	pancreas["age_1"] = strct;

	strct.a = 9.16 *cm;		
	strct.b = 0.90 *cm;		
	strct.c = 1.92 *cm;			
	strct.x0 = -0.57 *cm;
	strct.z0 = 21.57 *cm;
	strct.x1 = 1.72 *cm;
	strct.volume=22.7 *cm3;
	pancreas["age_5"] = strct;

	strct.a = 10.09 *cm;		
	strct.b = 0.92 *cm;		
	strct.c = 2.17 *cm;			
	strct.x0 = -0.38 *cm;//!!
	strct.z0 = 26.85 *cm;
	strct.x1 = 2.15 *cm;
	strct.volume=28.9 *cm3;
	pancreas["age_10"] = strct;

	strct.a = 13.32 *cm;		
	strct.b = 1.14 *cm;		
	strct.c = 2.87 *cm;			
	strct.x0 = -0.72 *cm;//!!
	strct.z0 = 33.35 *cm;
	strct.x1 = 2.61 *cm;
	strct.volume=62.4 *cm3;
	pancreas["age_15"] = strct;

	strct.a = 16.00 *cm;		
	strct.b = 1.20 *cm;		
	strct.c = 3.30 *cm;			
	strct.x0 = -1.00 *cm;//!!
	strct.z0 = 37.00 *cm;
	strct.x1 = 3.00 *cm;
	strct.volume=90.7 *cm3;
	pancreas["adult"] = strct;
}
void MirdMapConstants::initPelvis(){
	pelvis_struct strct;
	strct.a1 = 3.59 *cm;		
	strct.b1 = 5.54 *cm;		
	strct.a2 = 3.81 *cm;			
	strct.b2 = 5.88 *cm;
	strct.y01 = -1.86 *cm;
	strct.y02 = -1.47 *cm;
	strct.y1 = 2.45 *cm;
	strct.z1 = 4.32 *cm;
	strct.z2 = 6.79 *cm;
	strct.volume=28.9 *cm3;
	pelvis["newborn"] = strct;

	strct.a1 = 4.97 *cm;		
	strct.b1 = 7.35 *cm;		
	strct.a2 = 5.28 *cm;			
	strct.b2 = 7.80 *cm;
	strct.y01 = -2.47 *cm;
	strct.y02 = -1.95 *cm;
	strct.y1 = 3.25 *cm;
	strct.z1 = 6.14 *cm;
	strct.z2 = 9.65 *cm;
	strct.volume=76.0 *cm3;
	pelvis["age_1"] = strct;

	strct.a1 = 6.47 *cm;		
	strct.b1 = 8.48 *cm;		
	strct.a2 = 6.87 *cm;			
	strct.b2 = 9.00 *cm;
	strct.y01 = -2.85 *cm;
	strct.y02 = -2.25 *cm;
	strct.y1 = 3.75 *cm;
	strct.z1 = 8.16 *cm;
	strct.z2 = 12.82 *cm;
	strct.volume=151 *cm3;
	pelvis["age_5"] = strct;

	strct.a1 = 7.85 *cm;		
	strct.b1 = 9.49 *cm;		
	strct.a2 = 8.34 *cm;			
	strct.b2 = 10.08 *cm;
	strct.y01 = -3.19 *cm;
	strct.y02 = -2.52 *cm;
	strct.y1 = 4.20 *cm;
	strct.z1 = 10.16 *cm;
	strct.z2 = 15.97 *cm;
	strct.volume=258 *cm3;
	pelvis["age_10"] = strct;

	strct.a1 = 9.75 *cm;		
	strct.b1 = 11.07 *cm;		
	strct.a2 = 10.35 *cm;			
	strct.b2 = 11.76 *cm;
	strct.y01 = -3.72 *cm;
	strct.y02 = -2.94 *cm;
	strct.y1 = 4.90 *cm;
	strct.z1 = 12.62 *cm;
	strct.z2 = 19.83 *cm;
	strct.volume=460 *cm3;
	pelvis["age_15"] = strct;

	strct.a1 = 11.30 *cm;		
	strct.b1 = 11.30 *cm;		
	strct.a2 = 12.00 *cm;			
	strct.b2 = 12.00 *cm;
	strct.y01 = -3.80 *cm;
	strct.y02 = -3.00 *cm;
	strct.y1 = 5.00 *cm;
	strct.z1 = 14.00 *cm;
	strct.z2 = 22.00 *cm;
	strct.volume=606 *cm3;
	pelvis["adult"] = strct;
}
void MirdMapConstants::initRibCage(){
	ribCage_struct strct;
	strct.a = 5.40 *cm;		
	strct.b = 4.80 *cm;		
	strct.d = 0.21 *cm;			
	strct.z1 = 10.86 *cm;
	strct.z2 = 20.75 *cm;
	strct.c = 0.43 *cm;
	strct.volume=34 *cm3;
	ribCage["newborn"] = strct;

	strct.a = 7.48 *cm;		
	strct.b = 6.37 *cm;		
	strct.d = 0.28 *cm;			
	strct.z1 = 15.44 *cm;
	strct.z2 = 29.47 *cm;
	strct.c = 0.61 *cm;
	strct.volume=87.4 *cm3;
	ribCage["age_1"] = strct;

	strct.a = 9.73 *cm;		
	strct.b = 7.35 *cm;		
	strct.d = 0.34 *cm;			
	strct.z1 = 20.53 *cm;
	strct.z2 = 39.16 *cm;
	strct.c = 0.81 *cm;
	strct.volume=174 *cm3;
	ribCage["age_5"] = strct;

	strct.a = 11.82 *cm;		
	strct.b = 8.23 *cm;		
	strct.d = 0.39 *cm;			
	strct.z1 = 25.43 *cm;
	strct.z2 = 48.89 *cm;
	strct.c = 1.02 *cm;
	strct.volume=295 *cm3;
	ribCage["age_10"] = strct;

	strct.a = 14.66 *cm;		
	strct.b = 9.60 *cm;		
	strct.d = 0.47 *cm;			
	strct.z1 = 31.67 *cm;
	strct.z2 = 60.65 *cm;
	strct.c = 1.26 *cm;
	strct.volume=531 *cm3;
	ribCage["age_15"] = strct;

	strct.a = 17.00 *cm;		
	strct.b = 9.80 *cm;		
	strct.d = 0.50 *cm;			
	strct.z1 = 35.10 *cm;
	strct.z2 = 67.30 *cm;
	strct.c = 1.40 *cm;
	strct.volume=694 *cm3;
	ribCage["adult"] = strct;
}
void MirdMapConstants::initUpperLargeIntestine(){
	upperLargeIntestine_struct strct;
	strct.a = 0.79 *cm;		
	strct.b = 1.23 *cm;		
	strct.d = 0.27 *cm;	
	strct.x0 = -2.70 *cm;	
	strct.y0 = -1.16 *cm;	
	strct.z1 = 4.46 *cm;
	strct.z2 = 7.41 *cm;	
	strct.volume=9.01 *cm3;
	upperLargeIntestine["newborn"] = strct;

	strct.a = 1.10 *cm;		
	strct.b = 1.63 *cm;		
	strct.d = 0.37 *cm;	
	strct.x0 = -3.74 *cm;	
	strct.y0 = -1.53 *cm;	
	strct.z1 = 6.34 *cm;
	strct.z2 = 10.53 *cm;	
	strct.volume=23.6 *cm3;
	upperLargeIntestine["age_1"] = strct;

	strct.a = 1.43 *cm;		
	strct.b = 1.88 *cm;		
	strct.d = 0.46 *cm;	
	strct.x0 = -4.87 *cm;	
	strct.y0 = -1.77 *cm;	
	strct.z1 = 8.42 *cm;
	strct.z2 = 13.99 *cm;	
	strct.volume=47 *cm3;
	upperLargeIntestine["age_5"] = strct;

	strct.a = 1.74 *cm;		
	strct.b = 2.10 *cm;		
	strct.d = 0.54 *cm;	
	strct.x0 = -5.91 *cm;	
	strct.y0 = -1.98 *cm;	
	strct.z1 = 10.49 *cm;
	strct.z2 = 17.42 *cm;	
	strct.volume=79.6 *cm3;
	upperLargeIntestine["age_10"] = strct;

	strct.a = 2.16 *cm;		
	strct.b = 2.45 *cm;		
	strct.d = 0.65 *cm;	
	strct.x0 = -7.33 *cm;	
	strct.y0 = -2.31 *cm;	
	strct.z1 = 13.03 *cm;
	strct.z2 = 21.63 *cm;	
	strct.volume=142.9 *cm3;
	upperLargeIntestine["age_15"] = strct;

	strct.a = 2.50 *cm;		
	strct.b = 2.50 *cm;		
	strct.d = 0.7085 *cm;	
	strct.x0 = -8.50 *cm;	
	strct.y0 = -2.36 *cm;	
	strct.z1 = 14.45 *cm;
	strct.z2 = 24.00 *cm;	
	strct.volume=187.5 *cm3;
	upperLargeIntestine["adult"] = strct;
}
void MirdMapConstants::initTransverseColon(){
	transverseColon_struct strct;
	strct.b = 1.23 *cm;		
	strct.c = 0.46 *cm;		
	strct.d = 0.18 *cm;	
	strct.y0 = -1.16 *cm;	
	strct.z0 = 7.87 *cm;	
	strct.x1 = 3.33 *cm;	
	strct.volume=11.84 *cm3;
	transverseColon["newborn"] = strct;

	strct.b = 1.63 *cm;		
	strct.c = 0.65 *cm;		
	strct.d = 0.26 *cm;	
	strct.y0 = -1.53 *cm;	
	strct.z0 = 11.18 *cm;	
	strct.x1 = 4.62 *cm;	
	strct.volume=30.7 *cm3;
	transverseColon["age_1"] = strct;

	strct.b = 1.88 *cm;		
	strct.c = 0.87 *cm;		
	strct.d = 0.33 *cm;	
	strct.y0 = -1.77 *cm;	
	strct.z0 = 14.86 *cm;	
	strct.x1 = 6.01 *cm;	
	strct.volume=61.8 *cm3;
	transverseColon["age_5"] = strct;

	strct.b = 2.10 *cm;		
	strct.c = 1.08 *cm;		
	strct.d = 0.40 *cm;	
	strct.y0 = -1.98 *cm;	
	strct.z0 = 18.51 *cm;	
	strct.x1 = 7.30 *cm;	
	strct.volume=104 *cm3;
	transverseColon["age_10"] = strct;

	strct.b = 2.45 *cm;		
	strct.c = 1.35 *cm;		
	strct.d = 0.49 *cm;	
	strct.y0 = -2.31 *cm;	
	strct.z0 = 22.99 *cm;	
	strct.x1 = 9.06 *cm;	
	strct.volume=188.3 *cm3;
	transverseColon["age_15"] = strct;

	strct.b = 2.50 *cm;		
	strct.c = 1.50 *cm;		
	strct.d = 0.527 *cm;	
	strct.y0 = -2.36 *cm;	
	strct.z0 = 25.50 *cm;	
	strct.x1 = 10.50 *cm;	
	strct.volume=248 *cm3;
	transverseColon["adult"] = strct;
}
void MirdMapConstants::initLowerLargeIntestine(){
	lowerLargeIntestine_struct strct;
	strct.a = 0.60 *cm;		
	strct.b = 1.04 *cm;		
	strct.d = 0.20 *cm;	
	strct.x1 = 2.94 *cm;	
	strct.mx = 0.2477 *cm;	
	strct.my = 1.225 *cm;
	strct.z1 = 2.69 *cm;	
	strct.z2 = 7.41 *cm;	
	strct.volume=9.25 *cm3;
	lowerLargeIntestine["newborn"] = strct;

	strct.a = 0.83 *cm;		
	strct.b = 1.38 *cm;		
	strct.d = 0.27 *cm;	
	strct.x1 = 4.07 *cm;	
	strct.mx = 0.3432 *cm;	
	strct.my = 1.625 *cm;
	strct.z1 = 3.82 *cm;	
	strct.z2 = 10.53 *cm;	
	strct.volume=24.1 *cm3;
	lowerLargeIntestine["age_1"] = strct;

	strct.a = 1.08 *cm;		
	strct.b = 1.60 *cm;		
	strct.d = 0.34 *cm;	
	strct.x1 = 5.30 *cm;	
	strct.mx = 0.4466 *cm;	
	strct.my = 1.875 *cm;
	strct.z1 = 5.08 *cm;	
	strct.z2 = 13.99 *cm;	
	strct.volume=48.4 *cm3;
	lowerLargeIntestine["age_5"] = strct;

	strct.a = 1.31 *cm;		
	strct.b = 1.79 *cm;		
	strct.d = 0.40 *cm;	
	strct.x1 = 6.43 *cm;	
	strct.mx = 0.5421 *cm;	
	strct.my = 2.100 *cm;
	strct.z1 = 6.33 *cm;	
	strct.z2 = 17.42 *cm;	
	strct.volume=81.7 *cm3;
	lowerLargeIntestine["age_10"] = strct;

	strct.a = 1.62 *cm;		
	strct.b = 2.09 *cm;		
	strct.d = 0.49 *cm;	
	strct.x1 = 7.98 *cm;	
	strct.mx = 0.6728 *cm;	
	strct.my = 2.450 *cm;
	strct.z1 = 7.86 *cm;	
	strct.z2 = 21.63 *cm;	
	strct.volume=146.5 *cm3;
	lowerLargeIntestine["age_15"] = strct;

	strct.a = 1.88 *cm;		
	strct.b = 2.13 *cm;		
	strct.d = 0.54 *cm;	
	strct.x1 = 9.25 *cm;	
	strct.mx = 0.7800 *cm;	
	strct.my = 2.500 *cm;
	strct.z1 = 8.72 *cm;	
	strct.z2 = 24.00 *cm;	
	strct.volume=191.9 *cm3;
	lowerLargeIntestine["adult"] = strct;
}
void MirdMapConstants::initSigmoidColon(){
	sigmoidColon_struct strct;
	strct.a = 0.50 *cm;		
	strct.b = 0.77 *cm;		
	strct.d = 0.25 *cm;	
	strct.x0 = 0.95 *cm;	
	strct.z0 = 2.69 *cm;	
	strct.R1 = 1.77 *cm;	
	strct.R2 = 0.92 *cm;	
	strct.volume=5.12 *cm3;
	sigmoidColon["newborn"] = strct;

	strct.a = 0.69 *cm;		
	strct.b = 1.02 *cm;		
	strct.d = 0.34 *cm;	
	strct.x0 = 1.32 *cm;	
	strct.z0 = 3.82 *cm;	
	strct.R1 = 2.51 *cm;	
	strct.R2 = 1.31 *cm;	
	strct.volume=13.27 *cm3;
	sigmoidColon["age_1"] = strct;

	strct.a = 0.88 *cm;		
	strct.b = 1.21 *cm;		
	strct.d = 0.42 *cm;	
	strct.x0 = 1.72 *cm;	
	strct.z0 = 5.08 *cm;	
	strct.R1 = 3.33 *cm;	
	strct.R2 = 1.75 *cm;	
	strct.volume=26.71 *cm3;
	sigmoidColon["age_5"] = strct;

	strct.a = 0.96 *cm;		
	strct.b = 1.50 *cm;		
	strct.d = 0.48 *cm;	
	strct.x0 = 2.09 *cm;	
	strct.z0 = 6.33 *cm;	
	strct.R1 = 4.15 *cm;	
	strct.R2 = 2.18 *cm;	
	strct.volume=45 *cm3;
	sigmoidColon["age_10"] = strct;

	strct.a = 1.18 *cm;		
	strct.b = 1.76 *cm;		
	strct.d = 0.59 *cm;	
	strct.x0 = 2.59 *cm;	
	strct.z0 = 7.86 *cm;	
	strct.R1 = 5.16 *cm;	
	strct.R2 = 2.70 *cm;	
	strct.volume=80.6 *cm3;
	sigmoidColon["age_15"] = strct;

	strct.a = 1.57 *cm;		
	strct.b = 1.57 *cm;		
	strct.d = 0.66 *cm;	
	strct.x0 = 3.00 *cm;	
	strct.z0 = 8.72 *cm;	
	strct.R1 = 5.72 *cm;	
	strct.R2 = 3.00 *cm;	
	strct.volume=106 *cm3;
	sigmoidColon["adult"] = strct;
}
void MirdMapConstants::initThyroid(){
	thyroid_struct strct;
	strct.R = 0.87 *cm;		
	strct.r = 0.40 *cm;		
	strct.c = 2.00 *cm;	
	strct.y0 = -2.41 *cm;	
	strct.volume=1.24 *cm3;
	thyroid["newborn"] = strct;

	strct.R = 0.97 *cm;		
	strct.r = 0.44 *cm;		
	strct.c = 2.21 *cm;	
	strct.y0 = -2.87 *cm;	
	strct.volume=1.71 *cm3;
	thyroid["age_1"] = strct;

	strct.R = 1.21 *cm;		
	strct.r = 0.55 *cm;		
	strct.c = 2.76 *cm;	
	strct.y0 = -3.31 *cm;	
	strct.volume=3.32 *cm3;
	thyroid["age_5"] = strct;

	strct.R = 1.60 *cm;		
	strct.r = 0.73 *cm;		
	strct.c = 3.63 *cm;	
	strct.y0 = -3.56 *cm;	
	strct.volume=7.62 *cm3;
	thyroid["age_10"] = strct;

	strct.R = 1.85 *cm;		
	strct.r = 0.83 *cm;		
	strct.c = 4.20 *cm;	
	strct.y0 = -3.91 *cm;	
	strct.volume=11.9 *cm3;
	thyroid["age_15"] = strct;

	strct.R = 2.20 *cm;		
	strct.r = 1.00 *cm;		
	strct.c = 5.00 *cm;	
	strct.y0 = -4.00 *cm;	
	strct.volume=19.9 *cm3;
	thyroid["adult"] = strct;
}
void MirdMapConstants::initThymus(){
	thymus_struct strct;
	strct.a = 1.76 *cm;		
	strct.b = 0.70 *cm;		
	strct.c = 2.10 *cm;		
	strct.y0 = -3.60 *cm;
	strct.z0 = 19.30 *cm;
	strct.volume=10.08 *cm3;
	thymus["newborn"] = strct;

	strct.a = 1.75 *cm;		
	strct.b = 1.00 *cm;		
	strct.c = 3.00 *cm;		
	strct.y0 = -4.75 *cm;
	strct.z0 = 27.00 *cm;
	strct.volume=22 *cm3;
	thymus["age_1"] = strct;

	strct.a = 1.85 *cm;		
	strct.b = 1.05 *cm;		
	strct.c = 3.50 *cm;		
	strct.y0 = -5.48 *cm;
	strct.z0 = 35.00 *cm;
	strct.volume=28.5 *cm3;
	thymus["age_5"] = strct;

	strct.a = 1.85 *cm;		
	strct.b = 1.00 *cm;		
	strct.c = 3.90 *cm;		
	strct.y0 = -6.13 *cm;	
	strct.z0 = 43.00 *cm;
	strct.volume=30.2 *cm3;
	thymus["age_10"] = strct;

	strct.a = 1.75 *cm;		
	strct.b = 0.93 *cm;		
	strct.c = 4.00 *cm;		
	strct.y0 = -7.15 *cm;
	strct.z0 = 52.00 *cm;
	strct.volume=27.3 *cm3;
	thymus["age_15"] = strct;

	strct.a = 1.50 *cm;		
	strct.b = 0.80 *cm;		
	strct.c = 4.00 *cm;		
	strct.y0 = -7.30 *cm;
	strct.z0 = 57.00 *cm;
	strct.volume=20.1 *cm3;
	thymus["adult"] = strct;
}
void MirdMapConstants::initLiver(){
	liver_struct strct;
	strct.a = 5.19 *cm;		
	strct.b = 4.25 *cm;		
	strct.xm = 8.45 *cm;		
	strct.ym = 10.90 *cm;
	strct.zm = 13.27 *cm;
	strct.z1 = 8.33 *cm;
	strct.z2 = 13.27 *cm;
	strct.volume=117 *cm3;
	liver["newborn"] = strct;

	strct.a = 7.20 *cm;		
	strct.b = 5.47 *cm;		
	strct.xm = 12.83 *cm;		
	strct.ym = 16.55 *cm;
	strct.zm = 18.86 *cm;
	strct.z1 = 11.84 *cm;
	strct.z2 = 18.86 *cm;
	strct.volume=281 *cm3;
	liver["age_1"] = strct;

	strct.a = 9.39 *cm;		
	strct.b = 6.30 *cm;		
	strct.xm = 16.27 *cm;		
	strct.ym = 20.34 *cm;
	strct.zm = 25.06 *cm;
	strct.z1 = 15.74 *cm;
	strct.z2 = 25.06 *cm;
	strct.volume=562 *cm3;
	liver["age_5"] = strct;

	strct.a = 11.43 *cm;		
	strct.b = 6.83 *cm;		
	strct.xm = 21.98 *cm;		
	strct.ym = 29.67 *cm;
	strct.zm = 31.21 *cm;
	strct.z1 = 19.59 *cm;
	strct.z2 = 31.21 *cm;
	strct.volume=853 *cm3;
	liver["age_10"] = strct;

	strct.a = 14.19 *cm;		
	strct.b = 7.84 *cm;		
	strct.xm = 31.51 *cm;		
	strct.ym = 44.75 *cm;
	strct.zm = 38.76 *cm;
	strct.z1 = 24.34 *cm;
	strct.z2 = 38.76 *cm;
	strct.volume=1350 *cm3;
	liver["age_15"] = strct;

	strct.a = 16.50 *cm;		
	strct.b = 8.00 *cm;		
	strct.xm = 35.00 *cm;		
	strct.ym = 45.00 *cm;
	strct.zm = 43.00 *cm;
	strct.z1 = 27.00 *cm;
	strct.z2 = 43.00 *cm;
	strct.volume=1830 *cm3;
	liver["adult"] = strct;
}
void MirdMapConstants::initHeart(){
	heart_struct strct;
	strct.alpha1 = 0.5942;		
	strct.beta1 = -0.6421;		
	strct.gamma1 = -0.4845;		
	strct.alpha2 = -0.3291;		
	strct.beta2 = 0.3556;		
	strct.gamma2 = -0.8748;		
	strct.alpha3 = 0.7340;		
	strct.beta3 = 0.6792;		
	strct.gamma3 = 0.0;		
	strct.vx = 3.71 *cm;		
	strct.avy = 2.16 *cm;		
	strct.lavz = 1.34 *cm;		
	strct.ravz = 3.02 *cm;		
	strct.ax = 2.33 *cm;		
	strct.tlvw = 0.56 *cm;		
	strct.trvw = 0.26 *cm;		
	strct.tax = 0.13 *cm;		
	strct.x0 = 0.42 *cm;
	strct.y0 = -1.08 *cm;
	strct.z0 = 16.05 *cm;
	strct.volume=59.5 *cm3;
	heart["newborn"] = strct;

	strct.alpha1 = 0.6009;		
	strct.beta1 = -0.6216;		
	strct.gamma1 = -0.5025;		
	strct.alpha2 = -0.3493;		
	strct.beta2 = 0.3613;		
	strct.gamma2 = -0.8646;		
	strct.alpha3 = 0.7190;		
	strct.beta3 = 0.6950;		
	strct.gamma3 = 0.0;		
	strct.vx = 4.67 *cm;		
	strct.avy = 2.72 *cm;		
	strct.lavz = 1.68 *cm;		
	strct.ravz = 3.80 *cm;		
	strct.ax = 2.93 *cm;		
	strct.tlvw = 1.71 *cm;		
	strct.trvw = 0.33 *cm;		
	strct.tax = 0.16 *cm;		
	strct.x0 =0.54 *cm;
	strct.y0 = -1.67 *cm;
	strct.z0 = 22.43 *cm;
	strct.volume=118.6 *cm3;
	heart["age_1"] = strct;

	strct.alpha1 = 0.6237;		
	strct.beta1 = -0.5721;		
	strct.gamma1 = -0.5327;		
	strct.alpha2 = -0.3926;		
	strct.beta2 = 0.3601;		
	strct.gamma2 = -0.8463;		
	strct.alpha3 = 0.6760;		
	strct.beta3 = 0.7369;		
	strct.gamma3 = 0.0;		
	strct.vx = 5.72 *cm;		
	strct.avy = 3.33 *cm;		
	strct.lavz = 2.06 *cm;		
	strct.ravz = 4.66 *cm;		
	strct.ax = 3.59 *cm;		
	strct.tlvw = 0.86 *cm;		
	strct.trvw = 0.40 *cm;		
	strct.tax = 0.20 *cm;		
	strct.x0 =0.77 *cm;
	strct.y0 = -1.70 *cm;
	strct.z0 = 29.60 *cm;
	strct.volume=218.3 *cm3;
	heart["age_5"] = strct;

	strct.alpha1 = 0.6345;		
	strct.beta1 = -0.5370;		
	strct.gamma1 = -0.5559;		
	strct.alpha2 = -0.4243;		
	strct.beta2 = 0.3591;		
	strct.gamma2 = -0.8312;		
	strct.alpha3 = 0.6460;		
	strct.beta3 = 0.7633;		
	strct.gamma3 = 0.0;		
	strct.vx = 6.73 *cm;		
	strct.avy = 3.92 *cm;		
	strct.lavz = 2.43 *cm;		
	strct.ravz = 5.48 *cm;		
	strct.ax = 4.23 *cm;		
	strct.tlvw = 1.02 *cm;		
	strct.trvw = 0.47 *cm;		
	strct.tax = 0.23 *cm;		
	strct.x0 =0.80 *cm;
	strct.y0 = -1.70 *cm;
	strct.z0 = 36.60 *cm;
	strct.volume=355 *cm3;
	heart["age_10"] = strct;

	strct.alpha1 = 0.6453;		
	strct.beta1 = -0.5134;		
	strct.gamma1 = -0.5658;		
	strct.alpha2 = -0.4428;		
	strct.beta2 = 0.3523;		
	strct.gamma2 = -0.8245;		
	strct.alpha3 = 0.6226;		
	strct.beta3 = 0.7825;		
	strct.gamma3 = 0.0;		
	strct.vx = 7.86 *cm;		
	strct.avy = 4.57 *cm;		
	strct.lavz = 2.83 *cm;		
	strct.ravz = 6.40 *cm;		
	strct.ax = 4.94 *cm;		
	strct.tlvw = 1.19 *cm;		
	strct.trvw = 0.55 *cm;		
	strct.tax = 0.27 *cm;		
	strct.x0 =0.86 *cm;
	strct.y0 = -2.10 *cm;
	strct.z0 = 45.10 *cm;
	strct.volume=565 *cm3;
	heart["age_15"] = strct;

	strct.alpha1 = 0.6751;		
	strct.beta1 = -0.4727;		
	strct.gamma1 = -0.5664;		
	strct.alpha2 = -0.4640;		
	strct.beta2 = 0.3249;		
	strct.gamma2 = -0.8241;		
	strct.alpha3 = 0.5736;		
	strct.beta3 = 0.8191;		
	strct.gamma3 = 0.0;		
	strct.vx = 8.60 *cm;		
	strct.avy = 5.00 *cm;		
	strct.lavz = 3.10 *cm;		
	strct.ravz = 7.00 *cm;		
	strct.ax = 5.40 *cm;		
	strct.tlvw = 1.30 *cm;		
	strct.trvw = 0.60 *cm;		
	strct.tax = 0.30 *cm;		
	strct.x0 =1.00 *cm;
	strct.y0 = -1.80 *cm;
	strct.z0 = 50.00 *cm;
	strct.volume=740 *cm3;
	heart["adult"] = strct;
}
void MirdMapConstants::initSmallIntestine(){
	smallIntestine_struct strct;
	strct.a = 3.59 *cm;		
	strct.b = 5.54 *cm;		
	strct.y0 = -1.86 *cm;		
	strct.y1 = -2.39 *cm;
	strct.y2 = 1.08 *cm;
	strct.z1 = 5.25 *cm;
	strct.z2 = 8.33 *cm;
	strct.volume=50.9 *cm3;
	smallIntestine["newborn"] = strct;

	strct.a = 4.97 *cm;		
	strct.b = 7.35 *cm;		
	strct.y0 = -2.47 *cm;		
	strct.y1 = -3.16 *cm;
	strct.y2 = 1.43 *cm;
	strct.z1 = 7.46 *cm;
	strct.z2 = 11.84 *cm;
	strct.volume=132 *cm3;
	smallIntestine["age_1"] = strct;

	strct.a = 6.47 *cm;		
	strct.b = 8.48 *cm;		
	strct.y0 = -2.85 *cm;		
	strct.y1 = -3.65 *cm;
	strct.y2 = 1.65 *cm;
	strct.z1 = 9.91 *cm;
	strct.z2 = 15.74 *cm;
	strct.volume=265 *cm3;
	smallIntestine["age_5"] = strct;

	strct.a = 7.85 *cm;		
	strct.b = 9.49 *cm;		
	strct.y0 = -3.19 *cm;		
	strct.y1 = -4.08 *cm;
	strct.y2 = 1.85 *cm;
	strct.z1 = 12.34 *cm;
	strct.z2 = 19.59 *cm;
	strct.volume=447 *cm3;
	smallIntestine["age_10"] = strct;

	strct.a = 9.75 *cm;		
	strct.b = 11.07 *cm;		
	strct.y0 = -3.72 *cm;		
	strct.y1 = -4.76 *cm;
	strct.y2 = 2.16 *cm;
	strct.z1 = 15.32 *cm;
	strct.z2 = 24.34 *cm;
	strct.volume=806 *cm3;
	smallIntestine["age_15"] = strct;

	strct.a = 11.30 *cm;		
	strct.b = 11.30 *cm;		
	strct.y0 = -3.80 *cm;		
	strct.y1 = -4.86 *cm;
	strct.y2 = 2.20 *cm;
	strct.z1 = 17.00 *cm;
	strct.z2 = 27.00 *cm;
	strct.volume=1060 *cm3;
	smallIntestine["adult"] = strct;
}
void MirdMapConstants::initEsophagus(){
	esophagus_struct strct;
	strct.a = 0.35 *cm;		
	strct.b = 0.21 *cm;		
	strct.d = 0.14 *cm;		
	strct.y0 = 1.15 *cm;		
	strct.z2 = 13.27 *cm;		
	strct.z3 = 21.6 *cm;		
	strct.r = 0.23 *cm;		
	strct.x1p = 0.00 *cm;		
	strct.x2p = 3.03 *cm;
	strct.z1 = 12.93 *cm;
	strct.alpha1 = 0.571348;		
	strct.beta1 = -0.791789;		
	strct.gamma1 = -0.215943;		
	strct.alpha2 = 0.810922;		
	strct.beta2 = 0.585154;		
	strct.gamma2 = 0.0;		
	strct.alpha3 = 0.126360;		
	strct.beta3 = -0.175113;		
	strct.gamma3 = 0.976406;	
	esophagus["newborn"] = strct;

	strct.a = 0.51 *cm;		
	strct.b = 0.27 *cm;		
	strct.d = 0.19 *cm;		
	strct.y0 = 1.33 *cm;		
	strct.z2 = 18.86 *cm;		
	strct.z3 = 30.7 *cm;		
	strct.r = 0.34 *cm;		
	strct.x1p = 0.00 *cm;		
	strct.x2p = 3.86 *cm;
	strct.z1 = 18.40 *cm;
	strct.alpha1 = 0.639520;		
	strct.beta1 = -0.732178;		
	strct.gamma1 = -0.234370;		
	strct.alpha2 = 0.753155;		
	strct.beta2 = 0.657843;		
	strct.gamma2 = 0.0;		
	strct.alpha3 = 0.154178;		
	strct.beta3 = -0.176517;		
	strct.gamma3 = 0.972148;
	esophagus["age_1"] = strct;

	strct.a = 0.65 *cm;		
	strct.b = 0.32 *cm;		
	strct.d = 0.22 *cm;		
	strct.y0 = 1.69 *cm;		
	strct.z2 = 25.06 *cm;		
	strct.z3 = 40.8 *cm;		
	strct.r = 0.44 *cm;		
	strct.x1p = 0.00 *cm;		
	strct.x2p = 4.84 *cm;
	strct.z1 = 24.53 *cm;
	strct.alpha1 = 0.663545;		
	strct.beta1 = -0.701214;		
	strct.gamma1 = -0.260782;		
	strct.alpha2 = 0.726347;		
	strct.beta2 = 0.687328;		
	strct.gamma2 = 0.0;		
	strct.alpha3 = 0.179243;		
	strct.beta3 = -0.189418;		
	strct.gamma3 = 0.965398;
	esophagus["age_5"] = strct;

	strct.a = 0.79 *cm;		
	strct.b = 0.36 *cm;		
	strct.d = 0.25 *cm;		
	strct.y0 = 2.04 *cm;		
	strct.z2 = 31.21 *cm;		
	strct.z3 = 50.8 *cm;		
	strct.r = 0.52 *cm;		
	strct.x1p = 0.00 *cm;		
	strct.x2p = 5.75 *cm;
	strct.z1 = 30.62 *cm;
	strct.alpha1 = 0.678998;		
	strct.beta1 = -0.677776;		
	strct.gamma1 = -0.282102;		
	strct.alpha2 = 0.706470;		
	strct.beta2 = 0.707743;		
	strct.gamma2 = 0.0;		
	strct.alpha3 = 0.199655;		
	strct.beta3 = -0.199296;		
	strct.gamma3 = 0.959385;
	esophagus["age_10"] = strct;

	strct.a = 1.05 *cm;		
	strct.b = 0.40 *cm;		
	strct.d = 0.28 *cm;		
	strct.y0 = 2.29 *cm;		
	strct.z2 = 38.76 *cm;		
	strct.z3 = 63.1 *cm;		
	strct.r = 0.64 *cm;		
	strct.x1p = 0.00 *cm;		
	strct.x2p = 7.07 *cm;
	strct.z1 = 38.08 *cm;
	strct.alpha1 = 0.708385;		
	strct.beta1 = -0.637547;		
	strct.gamma1 = -0.302860;		
	strct.alpha2 = 0.668965;		
	strct.beta2 = 0.743294;		
	strct.gamma2 = 0.0;		
	strct.alpha3 = 0.225114;		
	strct.beta3 = -0.202603;		
	strct.gamma3 = 0.953035;
	esophagus["age_15"] = strct;

	strct.a = 1.17 *cm;		
	strct.b = 0.42 *cm;		
	strct.d = 0.30 *cm;		
	strct.y0 = 2.575 *cm;		
	strct.z2 = 43.00 *cm;		
	strct.z3 = 70.0 *cm;		
	strct.r = 0.70 *cm;		
	strct.x1p = 0.10 *cm;		
	strct.x2p = 7.80 *cm;
	strct.z1 = 42.30 *cm;
	strct.alpha1 = 0.736084;		
	strct.beta1 = -0.604969;		
	strct.gamma1 = -0.303634;		
	strct.alpha2 = 0.634945;		
	strct.beta2 = 0.772557;		
	strct.gamma2 = 0.0;		
	strct.alpha3 = 0.234575;		
	strct.beta3 = -0.192791;		
	strct.gamma3 = 0.952789;
	esophagus["adult"] = strct;
}
void MirdMapConstants::initGallBladder(){
	gallBladder_struct strct;
	strct.alpha1 = 0.9292;		
	strct.beta1 = 0.0;		
	strct.gamma1 = -0.3695;		
	strct.alpha2 = -0.1018;		
	strct.beta2 = 0.9613;		
	strct.gamma2 = -0.2559;		
	strct.alpha3 = 0.3553;		
	strct.beta3 = 0.2754;		
	strct.gamma3 = 0.8933;
	strct.r1 = 0.458 *cm;		
	strct.r2 = 0.500 *cm;		
	strct.s = 0.0 ;		
	strct.h = 3.10 *cm;		
	strct.x0 = -0.67 *cm;		
	strct.y0 = -1.75 *cm;		
	strct.z0 = 8.68 *cm;		
	strct.volume = 2.43 *cm3;					
	gallBladder["newborn"] = strct;
		
	strct.alpha1 = 0.9770;		
	strct.beta1 = 0.0;		
	strct.gamma1 = -0.2123;		
	strct.alpha2 = -0.0348;		
	strct.beta2 = 0.9866;		
	strct.gamma2 = -0.1594;		
	strct.alpha3 = 0.2105;		
	strct.beta3 = 0.1632;		
	strct.gamma3 = 0.9639;
	strct.r1 = 0.884 *cm;		
	strct.r2 = 0.937 *cm;		
	strct.s = 0.2275 ;		
	strct.h = 3.54 *cm;		
	strct.x0 = -0.71 *cm;		
	strct.y0 = -2.08 *cm;		
	strct.z0 = 13.16 *cm;		
	strct.volume = 5.50 *cm3;	
	gallBladder["age_1"] = strct;

	strct.alpha1 = 0.9814;		
	strct.beta1 = 0.0;		
	strct.gamma1 = -0.1921;		
	strct.alpha2 = -0.0291;		
	strct.beta2 = 0.9884;		
	strct.gamma2 = -0.1490;		
	strct.alpha3 = 0.1898;		
	strct.beta3 = 0.1518;		
	strct.gamma3 = 0.9700;
	strct.r1 = 1.414 *cm;		
	strct.r2 = 1.499 *cm;		
	strct.s = 0.2275 ;		
	strct.h = 5.66 *cm;		
	strct.x0 = -0.59 *cm;//!!!
	strct.y0 = -2.40 *cm;		
	strct.z0 = 17.49 *cm;		
	strct.volume = 22.50 *cm3;	
	gallBladder["age_5"] = strct;

	strct.alpha1 = 0.9722;		
	strct.beta1 = 0.0;		
	strct.gamma1 = -0.2342;		
	strct.alpha2 = -0.0400;		
	strct.beta2 = 0.9853;		
	strct.gamma2 = -0.1661;		
	strct.alpha3 = 0.2307;		
	strct.beta3 = 0.1709;		
	strct.gamma3 = 0.9579;
	strct.r1 = 1.768 *cm;		
	strct.r2 = 1.874 *cm;		
	strct.s = 0.2275 ;		
	strct.h = 7.07 *cm;		
	strct.x0 = -1.69 *cm;
	strct.y0 = -2.69 *cm;		
	strct.z0 = 21.77 *cm;		
	strct.volume = 44 *cm3;	
	gallBladder["age_10"] = strct;

	strct.alpha1 = 0.9550;		
	strct.beta1 = 0.0;		
	strct.gamma1 = -0.2964;		
	strct.alpha2 = -0.0606;		
	strct.beta2 = 0.9789;		
	strct.gamma2 = -0.1952;		
	strct.alpha3 = 0.2903;		
	strct.beta3 = 0.2044;		
	strct.gamma3 = 0.9349;
	strct.r1 = 1.916 *cm;		
	strct.r2 = 2.031 *cm;		
	strct.s = 0.2275 ;		
	strct.h = 7.66 *cm;		
	strct.x0 = -3.98 *cm;
	strct.y0 = -3.14 *cm;		
	strct.z0 = 27.04 *cm;		
	strct.volume = 56 *cm3;
	gallBladder["age_15"] = strct;

	strct.alpha1 = 0.9615;		
	strct.beta1 = 0.0;		
	strct.gamma1 = -0.2748;		
	strct.alpha2 = -0.0574;		
	strct.beta2 = 0.9779;		
	strct.gamma2 = -0.2008;		
	strct.alpha3 = 0.2687;		
	strct.beta3 = 0.2090;		
	strct.gamma3 = 0.9403;
	strct.r1 = 2.000 *cm;		
	strct.r2 = 2.120 *cm;		
	strct.s = 0.2275 ;		
	strct.h = 8.00 *cm;		
	strct.x0 = -4.50 *cm;
	strct.y0 = -3.20 *cm;		
	strct.z0 = 30.00 *cm;		
	strct.volume = 63.7 *cm3;
	gallBladder["adult"] = strct;
}
void MirdMapConstants::initClavicles(){
	clavicles_struct strct;
	strct.y0 = 0.73 *cm;		
	strct.z1 = 21.06 *cm;		
	strct.R = 5.07 *cm;		
	strct.r = 0.2833 *cm;		
	strct.cotan1 = 5.5868;		
	strct.cotan2 = 0.38510;		
	strct.volume = 2.62 *cm3;					
	clavicles["newborn"] = strct;
	
	strct.y0 = 1.38 *cm;		
	strct.z1 = 29.93 *cm;		
	strct.R = 7.14 *cm;		
	strct.r = 0.3930 *cm;		
	strct.cotan1 = 5.6814;		
	strct.cotan2 = 0.43161;		
	strct.volume = 6.85 *cm3;
	clavicles["age_1"] = strct;

	strct.y0 = 3.14 *cm;		
	strct.z1 = 39.78 *cm;		
	strct.R = 9.80 *cm;		
	strct.r = 0.4491 *cm;		
	strct.cotan1 = 5.9977;		
	strct.cotan2 = 0.56391;		
	strct.volume = 13.7 *cm3;
	clavicles["age_5"] = strct;

	strct.y0 = 4.93 *cm;		
	strct.z1 = 49.53 *cm;		
	strct.R = 12.40 *cm;		
	strct.r = 0.5981 *cm;		
	strct.cotan1 = 6.2581;		
	strct.cotan2 = 0.65708;		
	strct.volume = 23.2 *cm3;
	clavicles["age_10"] = strct;

	strct.y0 = 7.22 *cm;		
	strct.z1 = 61.52 *cm;		
	strct.R = 15.93 *cm;		
	strct.r = 0.7274 *cm;		
	strct.cotan1 = 6.4852;		
	strct.cotan2 = 0.73137;		
	strct.volume = 41.6 *cm3;
	clavicles["age_15"] = strct;

	strct.y0 = 11.10 *cm;		
	strct.z1 = 68.25 *cm;		
	strct.R = 20.00 *cm;		
	strct.r = 0.7883 *cm;		
	strct.cotan1 = 7.0342;		
	strct.cotan2 = 0.89415;		
	strct.volume = 54.7 *cm3;
	clavicles["adult"] = strct;
}
void MirdMapConstants::initFacial(){
	facial_struct strct;
	strct.a1 = 4.17 *cm;		
	strct.b1 = 5.43 *cm;		
	strct.d = 0.07 *cm;		
	strct.z1 = 2.16 *cm;		
	strct.z5 = 8.18 *cm;		
	strct.volume = 6.13 *cm3;					
	facial["newborn"] = strct;

	strct.a1 = 5.73 *cm;		
	strct.b1 = 7.44 *cm;		
	strct.d = 0.14 *cm;		
	strct.z1 = 2.93 *cm;		
	strct.z5 = 11.18 *cm;		
	strct.volume = 22.8 *cm3;					
	facial["age_1"] = strct;

	strct.a1 = 6.68 *cm;		
	strct.b1 = 8.60 *cm;		
	strct.d = 0.58 *cm;		
	strct.z1 = 3.30 *cm;		
	strct.z5 = 12.57 *cm;		
	strct.volume = 114 *cm3;	
	facial["age_5"] = strct;

	strct.a1 = 6.93 *cm;		
	strct.b1 = 8.90 *cm;		
	strct.d = 0.74 *cm;		
	strct.z1 = 3.61 *cm;		
	strct.z5 = 13.73 *cm;		
	strct.volume = 161 *cm3;	
	facial["age_10"] = strct;

	strct.a1 = 6.92 *cm;		
	strct.b1 = 8.91 *cm;		
	strct.d = 1.10 *cm;		
	strct.z1 = 3.79 *cm;		
	strct.z5 = 14.05 *cm;		
	strct.volume = 234 *cm3;	
	facial["age_15"] = strct;

	strct.a1 = 7.00 *cm;		
	strct.b1 = 9.00 *cm;		
	strct.d = 1.40 *cm;		
	strct.z1 = 4.00 *cm;		
	strct.z5 = 14.73 *cm;		
	strct.volume = 305 *cm3;	
	facial["adult"] = strct;
}