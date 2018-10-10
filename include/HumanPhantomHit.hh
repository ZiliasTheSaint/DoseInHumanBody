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
// $Id: G4HumanPhantomHit.hh,v 1.10 2007-05-15 14:46:54 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//D. Fulea, National Institute of Public Health Bucharest, Cluj-Napoca Regional Center, Romania

#ifndef HumanPhantomHit_h
#define HumanPhantomHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class HumanPhantomHit : public G4VHit
{
public:

  HumanPhantomHit();
  ~HumanPhantomHit();
  HumanPhantomHit(const HumanPhantomHit&);
  const HumanPhantomHit& operator=(const HumanPhantomHit&);
  G4int operator==(const HumanPhantomHit&) const;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  void Draw();
  void Print();

public:  
  void SetBodyPartID  (G4String bodyPartName) { bodyPartID = bodyPartName;};
  void SetEdep (G4double de) { edep = de; };
  void SetMass (G4double ma) { mass = ma; };//===============
      
  G4String GetBodyPartID() { return bodyPartID; };
  G4double GetEdep()    { return edep; }; 
  G4double GetMass()    { return mass; }; //================
      
private:
  G4String bodyPartID;
  G4double      edep;
  G4double      mass;//=========
};

typedef G4THitsCollection<HumanPhantomHit> HumanPhantomHitsCollection;

extern G4Allocator<HumanPhantomHit> HumanPhantomHitAllocator;

inline void* HumanPhantomHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) HumanPhantomHitAllocator.MallocSingle();
  return aHit;
}

inline void HumanPhantomHit::operator delete(void *aHit)
{
  HumanPhantomHitAllocator.FreeSingle((HumanPhantomHit*) aHit);
}

#endif
