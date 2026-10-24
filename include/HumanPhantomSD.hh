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

	  //->Grok Ai help
	  // Call after age / scale / sex are set. Also called lazily on first hit.
    void CacheOrganParameters();

 
 private:
      HumanPhantomHitsCollection* collection;
	  
	  HumanPhantomConstruction*    phantom;     //pointer to the geometry

   /*bool inThyroid(G4double, G4double, G4double);
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
   bool inFacial(G4double, G4double, G4double);*/

	  //->Grok ai help
	  bool paramsCached;

    // Thyroid
    G4double thy_ct, thy_c, thy_rr, thy_r, thy_y0;

    // Liver
    G4double liv_a, liv_b, liv_xm, liv_ym, liv_zm, liv_z1, liv_z2;

    // Heart
    G4double hrt_alpha1, hrt_beta1, hrt_gamma1;
    G4double hrt_alpha2, hrt_beta2, hrt_gamma2;
    G4double hrt_alpha3, hrt_beta3, hrt_gamma3;
    G4double hrt_x0, hrt_y0, hrt_z0;
    G4double hrt_vx, hrt_avy, hrt_lavz, hrt_ravz;
    G4double hrt_ax, hrt_tlvw, hrt_trvw, hrt_taw;

    // Clavicles
    G4double clv_cot1, clv_cot2, clv_y0, clv_z1, clv_rr, clv_r;

    // Gall bladder
    G4double glb_alpha1, glb_beta1, glb_gamma1;
    G4double glb_alpha2, glb_beta2, glb_gamma2;
    G4double glb_alpha3, glb_beta3, glb_gamma3;
    G4double glb_r1, glb_r2, glb_s, glb_h;
    G4double glb_x0, glb_y0, glb_z0;

    // Small intestine
    G4double smi_a, smi_b, smi_y0, smi_y1, smi_y2, smi_z1, smi_z2;

    // Thymus
    G4double thm_a, thm_b, thm_c, thm_y0, thm_z0;

    // Testes
    G4double tst_a, tst_b, tst_c, tst_y0;

    // Esophagus
    G4double esp_alpha1, esp_beta1, esp_gamma1;
    G4double esp_alpha2, esp_beta2, esp_gamma2;
    G4double esp_alpha3, esp_beta3, esp_gamma3;
    G4double esp_a, esp_y0, esp_b, esp_d;
    G4double esp_z1, esp_z2, esp_z3, esp_r, esp_xp1, esp_xp2;

    // Facial
    G4double fac_ct, fac_ch1;
    G4double fac_aS, fac_bS, fac_cS, fac_dS;
    G4double fac_a1, fac_b1, fac_d, fac_z1, fac_z5;
	  //max efficiency const - read only function! these functions cannot alter class member like collection or phantom
	  bool inThyroid       (G4double X, G4double Y, G4double Z) const;
    bool inLiver         (G4double X, G4double Y, G4double Z) const;
    bool inHeart         (G4double X, G4double Y, G4double Z) const;
    bool inLeftClavicle  (G4double X, G4double Y, G4double Z) const;
    bool inRightClavicle (G4double X, G4double Y, G4double Z) const;
    bool inGallBladder   (G4double X, G4double Y, G4double Z) const;
    bool inSmallIntestine(G4double X, G4double Y, G4double Z) const;
    bool inThymus        (G4double X, G4double Y, G4double Z) const;
    bool inLeftTesticle  (G4double X, G4double Y, G4double Z) const;
    bool inRightTesticle (G4double X, G4double Y, G4double Z) const;
    bool inEsophagus     (G4double X, G4double Y, G4double Z) const;
    bool inFacial        (G4double X, G4double Y, G4double Z) const;
};
#endif

