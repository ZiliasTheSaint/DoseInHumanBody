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
//D. Fulea, National Institute of Public Health Bucharest, Cluj-Napoca Regional Center, Romania
//
#ifndef HumanPhantomSD_h
#define HumanPhantomSD_h 1

#include "G4VSensitiveDetector.hh"
#include "HumanPhantomHit.hh"
#include "globals.hh"
#include <map>

class G4Step;

class HumanPhantomConstruction;

class HumanPhantomSD : public G4VSensitiveDetector
{
  public:
      HumanPhantomSD(G4String);
     ~HumanPhantomSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*, G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
 
 private:
      HumanPhantomHitsCollection* collection;
	  
	  HumanPhantomConstruction*    phantom;     //pointer to the geometry

	   bool inThyroid(G4double, G4double, G4double);
   bool inLiver(G4double, G4double, G4double);
   bool inHeart(G4double, G4double, G4double);
   bool inLeftClavicle(G4double, G4double, G4double);
   bool inRightClavicle(G4double, G4double, G4double);
   bool inGallBladder(G4double, G4double, G4double);
   bool inSmallIntestine(G4double, G4double, G4double);
   bool inThymus(G4double, G4double, G4double);
   bool inLeftTesticle(G4double, G4double, G4double);
   bool inRightTesticle(G4double, G4double, G4double);
   bool inEsophagus(G4double, G4double, G4double);
   bool inFacial(G4double, G4double, G4double);
};
#endif

