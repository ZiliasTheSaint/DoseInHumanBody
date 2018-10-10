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

