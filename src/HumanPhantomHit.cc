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
//
// $Id: G4HumanPhantomHit.cc,v 1.11 2007-05-15 14:45:35 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// D. Fulea, National Institute of Public Health Bucharest, Cluj-Napoca Regional Center, Romania

#include "HumanPhantomHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<HumanPhantomHit> HumanPhantomHitAllocator;


HumanPhantomHit::HumanPhantomHit() {}

HumanPhantomHit::~HumanPhantomHit() {}

//const HumanPhantomHit& right - trimis prin referinta nu prin valoare
//const asigura ca obiectul HumanPhantomHit nu poate fi modificat aici

/*
cpp theory:
...daca ai fi trimis prin valoare,
//trimiterea prin valoare inseamna ca programul face o copie a obiectului right iar pentru copie ar fi
//trebuit sa apeleze constructorul care e fix acesta->crash infinit loop!
//'right' s-a nascut in HumanDetectorSD aici: HumanPhantomHit* newHit = new HumanPhantomHit();
//prin apelarea empty constructor!!! Acest constructor intra apoi in scena prin collection->insert(newHit);
Collection e definita in hh: typedef G4THitsCollection<HumanPhantomHit> HumanPhantomHitsCollection;
ca alias la clase Geant4! si e pointata in HumanDetectorSD.hh: HumanPhantomHitsCollection* collection;

Referinta nu e un pointer si nu e obiectul in sine ci un ALIAS!?!
Compilatorul ascunde adresa din spate, permitandu-ti sa lucrezi cu ea folosind punctul
(right.edep), ca si cum ai avea obiectul real.
Pe scurt: In hardware cand se executa codul, o referinta e o adresa de memorie, un pointer deghizat.
PC-ul are nevoie de el sa ajunga la obiectul original.
Cand se apeleaza right.edep compilatorul executa o dereferentiere implicita...citeste adresa din spatele aliasului
right -> gaseste originalul.

daca ar fi pointer s-ar declara asa HumanPhantomHit* right
Pointerul este o variabila separata care contine o adresa de memorie. Pentru a ajunge la date, esti
obligat sa folosesti sageata: right->edep.
*/
HumanPhantomHit::HumanPhantomHit(const HumanPhantomHit& right)
  : G4VHit()
{
  bodyPartID = right.bodyPartID;
  edep = right.edep;
  mass = right.mass;//---------------
}

const HumanPhantomHit& HumanPhantomHit::operator=(const HumanPhantomHit& right)
{
  bodyPartID = right.bodyPartID;
  edep  = right.edep;
  mass = right.mass;//---------------
  return *this;
}

G4int HumanPhantomHit::operator==(const HumanPhantomHit& right) const
{
  return (this==&right) ? 1 : 0;
}

void HumanPhantomHit::Draw()
{
}

void HumanPhantomHit::Print()
{
    G4cout << "Energy deposit: " << G4BestUnit(edep,"Energy")
	 << "BodyPartID: " << bodyPartID 
	 << "BodyPartMass: " << G4BestUnit(mass,"Mass") 
	 << G4endl;
}

