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
#ifndef MaleBuilder_h
#define MaleBuilder_h 1

#include "PhantomBuilder.hh"

class PhantomBuilder;
class MaleBuilder: public PhantomBuilder
{
public:
  MaleBuilder();
  ~MaleBuilder();

  G4double getMassOfLeftBreast(){return 1.0;};
  G4double getMassOfRightBreast(){return 1.0;};
  G4double getMassOfLeftOvary(){return 1.0;};
  G4double getMassOfRightOvary(){return 1.0;};
  G4double getMassOfUterus(){return 1.0;};

  
  void BuildLeftTesticle(const G4String&,G4bool,G4bool);
  void BuildRightTesticle(const G4String&,G4bool,G4bool);  
  G4double getMassOfLeftTesticle(){return volOfLeftTesticle*rhoOfLeftTesticle;};
  G4double getMassOfRightTesticle(){return volOfRightTesticle*rhoOfRightTesticle;};

  //testes are NOT in trunk,legs or head (maybe here for someones!!)! HA HA HA!
protected:
	G4double volOfLeftTesticle;G4double rhoOfLeftTesticle;
	G4double volOfRightTesticle;G4double rhoOfRightTesticle;
};
#endif
