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
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
//

#include "G4VisAttributes.hh"
#include "HumanPhantomColour.hh"
#include "G4Colour.hh"
HumanPhantomColour::HumanPhantomColour()
{
  white = G4Colour(1.0, 1.0, 1.0);
  pink = G4Colour(0.94, 0.5, 0.5);
  grey = G4Colour(0.46, 0.53, 0.6);
  yellow = G4Colour(1.0, 1.0, 0.);
  blue = G4Colour(0.25,0.41, 0.88 );
  lightBlue = G4Colour(0.28, 0.82, 0.8);
  green = G4Colour(0., 1., 0.);
  brown = G4Colour(0.5, 0.5, 0.);
  purple = G4Colour(0.85,0.44,0.84);
  red = G4Colour(1.0, 0.0, 0.0);
  orange = G4Colour(1.,0.5,0.); 
  black =  G4Colour(0.,0.,0.); 
}

HumanPhantomColour::~HumanPhantomColour()
{ 

}

G4Colour HumanPhantomColour::GetColour(const G4String& colourName)
{
  
 // Returns the colour
  if (colourName == "pink") return pink;
  else if(colourName == "white") return white;
  else if (colourName == "grey") return grey;
  else if (colourName == "yellow") return yellow;
  else if (colourName == "blue") return blue;
  else if (colourName == "lightBlue") return lightBlue;
  else if (colourName == "green") return green;
  else if (colourName == "brown") return brown;
  else if (colourName == "purple") return purple;
  else if (colourName == "red") return red;
  else if (colourName == "orange") return orange;
  else if  (colourName == "black") return black; 
else {G4cout<< colourName << "does not exist !!!"<< G4endl; return white;}
}
