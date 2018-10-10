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
//D. Fulea, National Institute of Public Health Bucharest, Cluj-Napoca Regional Center, Romania
//

#ifndef HumanPhantomRunAction_h
#define HumanPhantomRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <map>
#include <vector>
//============
class HumanPhantomConstruction;
class HumanPhantomPrimaryGeneratorAction;//!!!!!!!!!!!!!!
class XRayBuilder;
//=========
class HumanPhantomRunAction : public G4UserRunAction
{
  public:
    HumanPhantomRunAction();
   ~HumanPhantomRunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
    void Fill(G4String bodypartName, G4double energyDeposit);

	///////////
	void fillPerEvent(G4double);
	//===============

	//void FillMass(G4String bodypartName, G4double mass);//==================

	//void SetPhantomTotalMass(G4double mass){
	//	phantomTotalMass=mass;
	//}
	
	G4double basicAnalyze(G4double x, G4double x2, G4int n){
	  G4double temp=x;G4double temp2=x2;G4int nevents=n;
	  temp=temp/nevents;temp2=temp2/nevents;
	  temp2 = temp2 - temp*temp;
	  if (nevents>1) temp2=temp2/(nevents-1.0);
	  if (temp2 >0.)
	   temp2 = sqrt(temp2); 
	  //else temp2 = 99.99;//never!
	  //=========percent
	  if (temp!=0.0){
	   temp2 = std::min(100.0*temp2/temp,99.9);
      } //else temp2 = 99.9;//no score means no score not necessarly an error!
	  return temp2;
  }

private:
    void totalRunEnergyDeposit(G4int n);
    std::map<std::string,G4double> energyTotal;
	
	//std::map<std::string,G4double> bodypartMass;
	std::map<std::string,double> bodypartWeightMap;
	std::map<std::string,G4double> bodypartMassMap;
	std::map<std::string,G4double> energy;//=======
    std::map<std::string,G4double> energy2;//=======
	std::map<std::string,G4double> dose;//=======
	    
	HumanPhantomConstruction*    phantom;     //pointer to the geometry
	HumanPhantomPrimaryGeneratorAction* gun; //pointer to source

	G4double totalPhantomEnergy;
	G4double totalPhantomEnergy2;
	
	bool isGoodEnoughToPrint(G4String);
	
	void prepareAliasSampling(int nsbin, std::vector<double> fs_array);
	void prepareAliasSampling(int nsbin, double* fs_array);
	int ndata;
	double* enV_array;
	double* ws_array;
	int* ibin_array;
	double photonsPerDAP;
	double kermaPer_mas_at75cm;
	bool dummy;
	XRayBuilder* xraybuilder;

	G4double sumEDet, sum2EDet;

	double linInt(double x1, double y1, double x2, double y2,double x);
	int findNearestDataIndexFromArray(double* a, int arraySize, double value, bool lowerThanValue);
	
	double getScaleFactor(int numberOfEvents);
};
#endif





