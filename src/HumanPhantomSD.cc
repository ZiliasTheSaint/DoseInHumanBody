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
#include "HumanPhantomSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "HumanPhantomConstruction.hh"
#include "G4RunManager.hh"

#include "MirdMapConstants.hh"

HumanPhantomSD::HumanPhantomSD(G4String name)
  :G4VSensitiveDetector(name),
  paramsCached(false)//Grok ai!!!!-> initialization with false here!!!!!
{
 G4String HCname;
 collectionName.insert(HCname="HumanPhantomCollection");

 phantom = (HumanPhantomConstruction*)
             G4RunManager::GetRunManager()->GetUserDetectorConstruction();
}

HumanPhantomSD::~HumanPhantomSD()
{ 	
}

void HumanPhantomSD::Initialize(G4HCofThisEvent* HCE)
{
  collection = new HumanPhantomHitsCollection
                          (SensitiveDetectorName,collectionName[0]); 
  static G4int HCID = -1;
  if(HCID<0)
  { 
  HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( HCID, collection ); 
}

//->Grok ai added

// ---------------------------------------------------------------------------
// Cache all organ parameters once (no allocations inside ProcessHits)
// ---------------------------------------------------------------------------
void HumanPhantomSD::CacheOrganParameters()
{
  if (!phantom) return;

  const G4double scaleXY  = phantom->GetScaleXY();
  const G4double scaleZ   = phantom->GetScaleZ();
  const G4String ageGroup = phantom->GetAgeGroup();

  MirdMapConstants mmc;   // on stack – destroyed automatically
  /*
  CPP theory: 
  MirdMapConstants* mmc = new MirdMapConstants();
	mmc->initLiver(); //mmc este un pointer la un obiect, are steluta (*) si e localizat in memoria globala heap.
	pointerii au nevoie sa fie DELETE manual

	MirdMapConstants mmc; // Creat direct pe stiva (stack), fara steluta (*)
mmc.initLiver();      // Se apeleaza cu PUNCT...NU e pointer la obiect ci e OBIECTUL pur si simplu.
Stiva (stack) e memoria alocata ACESTEI functii nu e heap-ul global! 
Pointerii ajung in heap si de aia tre DELETE.

static MirdMapConstants mmc; // Creat o singura data in memoria globala
mmc.initLiver();             // Se apeleaza tot cu PUNCT. e OBIECT
OBIECTELE nu au nevoie de delete...c++ le delete auto la finalul functiei!

in c we control where objects resides...in stack or in heap!

Sintax
auto trk = mmc.GetTrunkMap();
auto tells compilator to figure out data type of trk based on what mmc.GetTrunkMap() returns
thy_ct = scaleZ * trk[ageGroup].Ct;
type of thy_ct is in hh!!!
  */

  // --- Thyroid ---
  mmc.initTrunk();
  auto trk = mmc.GetTrunkMap();
  thy_ct = scaleZ * trk[ageGroup].Ct;

  mmc.initThyroid();
  auto thy = mmc.GetThyroidMap();
  thy_c  = scaleZ  * thy[ageGroup].c;
  thy_rr = scaleXY * thy[ageGroup].R;
  thy_r  = scaleXY * thy[ageGroup].r;
  thy_y0 = scaleXY * thy[ageGroup].y0;
  thy_y0 = -thy_y0;   // MOVED IN THE BACK TO MATCH..NECK ADDED

  // --- Liver ---
  mmc.initLiver();
  auto lvr = mmc.GetLiverMap();
  liv_a  = scaleXY * lvr[ageGroup].a;
  liv_b  = scaleXY * lvr[ageGroup].b;
  liv_xm = scaleXY * lvr[ageGroup].xm;
  liv_ym = scaleXY * lvr[ageGroup].ym;
  liv_zm = scaleZ  * lvr[ageGroup].zm;
  liv_z1 = scaleZ  * lvr[ageGroup].z1;
  liv_z2 = scaleZ  * lvr[ageGroup].z2;

  // --- Heart ---
  mmc.initHeart();
  auto hrt = mmc.GetHeartMap();
  hrt_alpha1 = hrt[ageGroup].alpha1;
  hrt_beta1  = hrt[ageGroup].beta1;
  hrt_gamma1 = hrt[ageGroup].gamma1;
  hrt_alpha2 = hrt[ageGroup].alpha2;
  hrt_beta2  = hrt[ageGroup].beta2;
  hrt_gamma2 = hrt[ageGroup].gamma2;
  hrt_alpha3 = hrt[ageGroup].alpha3;
  hrt_beta3  = hrt[ageGroup].beta3;
  hrt_gamma3 = hrt[ageGroup].gamma3;
  hrt_x0   = scaleXY * hrt[ageGroup].x0;
  hrt_y0   = scaleXY * hrt[ageGroup].y0;
  hrt_z0   = scaleZ  * hrt[ageGroup].z0;
  hrt_vx   = scaleXY * hrt[ageGroup].vx;
  hrt_avy  = scaleXY * hrt[ageGroup].avy;
  hrt_lavz = scaleZ  * hrt[ageGroup].lavz;
  hrt_ravz = scaleZ  * hrt[ageGroup].ravz;
  hrt_ax   = scaleXY * hrt[ageGroup].ax;
  hrt_tlvw = scaleXY * hrt[ageGroup].tlvw;
  hrt_trvw = scaleXY * hrt[ageGroup].trvw;
  hrt_taw  = scaleXY * hrt[ageGroup].tax;

  // --- Clavicles ---
  mmc.initClavicles();
  auto clv = mmc.GetClaviclesMap();
  clv_cot1 = clv[ageGroup].cotan1;
  clv_cot2 = clv[ageGroup].cotan2;
  clv_y0   = scaleXY * clv[ageGroup].y0;
  clv_z1   = scaleZ  * clv[ageGroup].z1;
  clv_rr   = scaleXY * clv[ageGroup].R;
  clv_r    = scaleZ  * clv[ageGroup].r;

  // --- Gall bladder ---
  mmc.initGallBladder();
  auto glb = mmc.GetGallBladderMap();
  glb_alpha1 = glb[ageGroup].alpha1;
  glb_beta1  = glb[ageGroup].beta1;
  glb_gamma1 = glb[ageGroup].gamma1;
  glb_alpha2 = glb[ageGroup].alpha2;
  glb_beta2  = glb[ageGroup].beta2;
  glb_gamma2 = glb[ageGroup].gamma2;
  glb_alpha3 = glb[ageGroup].alpha3;
  glb_beta3  = glb[ageGroup].beta3;
  glb_gamma3 = glb[ageGroup].gamma3;
  glb_r1 = scaleZ  * glb[ageGroup].r1;
  glb_r2 = scaleXY * glb[ageGroup].r2;
  glb_s  = glb[ageGroup].s;
  glb_h  = scaleZ  * glb[ageGroup].h;
  glb_x0 = scaleXY * glb[ageGroup].x0;
  glb_y0 = scaleXY * glb[ageGroup].y0;
  glb_z0 = scaleZ  * glb[ageGroup].z0;

  // --- Small intestine ---
  mmc.initSmallIntestine();
  auto smi = mmc.GetSmallIntestineMap();
  smi_a  = scaleXY * smi[ageGroup].a;
  smi_b  = scaleXY * smi[ageGroup].b;
  smi_y0 = scaleXY * smi[ageGroup].y0;
  smi_y1 = scaleXY * smi[ageGroup].y1;
  smi_y2 = scaleXY * smi[ageGroup].y2;
  smi_z1 = scaleZ  * smi[ageGroup].z1;
  smi_z2 = scaleZ  * smi[ageGroup].z2;

  // --- Thymus ---
  mmc.initThymus();
  auto thm = mmc.GetThymusMap();
  thm_a  = scaleXY * thm[ageGroup].a;
  thm_b  = scaleXY * thm[ageGroup].b;
  thm_c  = scaleZ  * thm[ageGroup].c;
  thm_y0 = scaleXY * thm[ageGroup].y0;
  thm_z0 = scaleZ  * thm[ageGroup].z0;

  // --- Testes ---
  mmc.initTestes();
  auto tst = mmc.GetTestesMap();
  tst_a  = scaleXY * tst[ageGroup].a;
  tst_b  = scaleXY * tst[ageGroup].b;
  tst_c  = scaleZ  * tst[ageGroup].c;
  tst_y0 = scaleXY * tst[ageGroup].y0;

  // --- Esophagus ---
  mmc.initEsophagus();
  auto esp = mmc.GetEsophagusMap();
  esp_alpha1 = esp[ageGroup].alpha1;
  esp_beta1  = esp[ageGroup].beta1;
  esp_gamma1 = esp[ageGroup].gamma1;
  esp_alpha2 = esp[ageGroup].alpha2;
  esp_beta2  = esp[ageGroup].beta2;
  esp_gamma2 = esp[ageGroup].gamma2;
  esp_alpha3 = esp[ageGroup].alpha3;
  esp_beta3  = esp[ageGroup].beta3;
  esp_gamma3 = esp[ageGroup].gamma3;
  esp_a   = scaleXY * esp[ageGroup].a;
  esp_y0  = scaleXY * esp[ageGroup].y0;
  esp_b   = scaleXY * esp[ageGroup].b;
  esp_d   = scaleXY * esp[ageGroup].d;
  esp_z1  = scaleZ  * esp[ageGroup].z1;
  esp_z2  = scaleZ  * esp[ageGroup].z2;
  esp_z3  = scaleZ  * esp[ageGroup].z3;
  esp_r   = scaleXY * esp[ageGroup].r;
  esp_xp1 = scaleXY * esp[ageGroup].x1p;
  esp_xp2 = scaleXY * esp[ageGroup].x2p;

  // --- Facial ---
  mmc.initTrunk();
  trk = mmc.GetTrunkMap();
  fac_ct = scaleZ * trk[ageGroup].Ct;

  mmc.initHead();
  auto hed = mmc.GetHeadMap();
  fac_ch1 = scaleZ * hed[ageGroup].Ch1;

  mmc.initBrain();
  auto brn = mmc.GetBrainMap();
  fac_aS = scaleXY * brn[ageGroup].a;
  fac_bS = scaleXY * brn[ageGroup].b;
  fac_cS = scaleZ  * brn[ageGroup].c;

  mmc.initSkull();
  auto skl = mmc.GetSkullMap();
  fac_dS = scaleXY * skl[ageGroup].d;

  mmc.initFacial();
  auto fcl = mmc.GetFacialMap();
  fac_a1 = scaleXY * fcl[ageGroup].a1;
  fac_b1 = scaleXY * fcl[ageGroup].b1;
  fac_d  = scaleXY * fcl[ageGroup].d;
  fac_z1 = scaleZ  * fcl[ageGroup].z1;
  fac_z5 = scaleZ  * fcl[ageGroup].z5;

  paramsCached = true;
}

G4bool HumanPhantomSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{  
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep==0.) return false;

   //-> Grok - Lazy cache – safe if age/scale were set after SD construction
  if (!paramsCached) CacheOrganParameters();
  
  G4double X=aStep->GetTrack()->GetPosition().x();
  G4double Y=aStep->GetTrack()->GetPosition().y();
  G4double Z=aStep->GetTrack()->GetPosition().z();

  
  
  G4String bodypartName = aStep->GetPreStepPoint()->GetTouchable()
   			->GetVolume()->GetLogicalVolume()->GetName();

  //==============here we must check inThyroid....
  //problem: virtual male genitalia region!!!!!!!!not score as hits!!!!must add it to SD!!
		if (inThyroid(X,Y,Z))
		   bodypartName="logicalThyroid";
		else if (inLiver(X,Y,Z))
		   bodypartName="logicalLiver";//////
		else if (inHeart(X,Y,Z))
		   bodypartName="logicalHeart";////////////
		else if (inLeftClavicle(X,Y,Z))
		   bodypartName="logicalLeftClavicle";/////////
		else if (inRightClavicle(X,Y,Z))
		   bodypartName="logicalRightClavicle";///////////
		else if (inGallBladder(X,Y,Z))
		   bodypartName="logicalGallBladder";/////
        else if (inSmallIntestine(X,Y,Z))
		   bodypartName="logicalSmallIntestine";////
		else if (inThymus(X,Y,Z))
		   bodypartName="logicalThymus";
		else if (inEsophagus(X,Y,Z))
		   bodypartName="logicalEsophagus";////////////////
		else if (inLeftTesticle(X,Y,Z)){
			//G4String sex = phantom->GetPhantomSex();
			if(phantom->GetPhantomSex()=="Male")//if(sex=="Male")
		      bodypartName="logicalLeftTesticle";////
		}
		else if (inRightTesticle(X,Y,Z)){
			//G4String sex = phantom->GetPhantomSex();
			if(phantom->GetPhantomSex()=="Male")//if(sex=="Male")
		      bodypartName="logicalRightTesticle";/////
		}
		else if (inFacial(X,Y,Z)){
	      bodypartName="logicalFacialBones";/////
		}

  //=================
  G4double bodypartMassd = 0.0;
  if (bodypartName=="logicalmaleGenitalia"){
	  //virtual
	  bodypartMassd =1.0 * kg;//
  }else
	bodypartMassd=phantom->GetBodypartMass(bodypartName);  
  //==================

  HumanPhantomHit* newHit = new HumanPhantomHit();
  newHit->SetEdep(edep);
  newHit->SetBodyPartID(bodypartName);
  newHit->SetMass(bodypartMassd);
  collection->insert(newHit);
  
 return true;
}

void HumanPhantomSD::EndOfEvent(G4HCofThisEvent*)
{

 /*G4int NbHits = collection->entries();
      G4cout << "\n-------->Hits Collection: in this event they are " << NbHits 
             << " hits in the tracker chambers: " << G4endl;
      for (G4int i=0;i<NbHits;i++) (*collection)[i]->Print();*/
}

//==================================================================
//Grok -> all below functions are const!!!!!!
bool HumanPhantomSD::inThyroid(G4double X, G4double Y, G4double Z) const{	
	/*
	G4double scaleXY=phantom->GetScaleXY();
	G4double scaleZ=phantom->GetScaleZ();
	*/

	//CONCLUSION: xMIRD=x;yMIRD=-z;zMIRD=y;
	const G4double x = X, y = -Z, z = Y;//G4double x=X;G4double y=-Z;G4double z=Y;
	//------------------------------------------
	/*G4String ageGroup=phantom->GetAgeGroup();*/
	//MirdMapConstants* mmc = new MirdMapConstants();
	//->Gemini static to be reused!!!!! NEW keyword for static is not necessary in cpp!
	/*static MirdMapConstants mmc;
	mmc.initTrunk();//mmc->initTrunk();
  std::map<std::string,trunk_struct> trk = mmc.GetTrunkMap();//mmc->GetTrunkMap();  
  G4double bt = scaleXY*trk[ageGroup].Bt;
  G4double ct = scaleZ*trk[ageGroup].Ct;  
  mmc.initThyroid();//mmc->initThyroid();
  std::map<std::string,thyroid_struct> thy = mmc.GetThyroidMap();  //mmc->GetThyroidMap();  
  G4double c = scaleZ*thy[ageGroup].c;
  G4double rr = scaleXY*thy[ageGroup].R;
  G4double r = scaleXY*thy[ageGroup].r;
  G4double y0 = scaleXY*thy[ageGroup].y0;  
  //delete mmc;
    y0=-y0;//MOVED IN THE BACK TO MATCH..NECK ADDED!!!!!!!!!!!
	//-------------------------------
	//G4double ct=scaleZ*70.0*cm;//top trunk coordonate
	//G4double c=scaleZ*5.00*cm;
	//G4double rr = scaleXY*2.2*cm;// R
	//G4double r = scaleXY*1.0*cm;// r
	//G4double y0 = scaleXY*4.0*cm;//-scaleXY*4.0*cm;//MOVED IN THE BACK TO MATCH..NECK ADDED!!!!!!!!!!!
	//---------------------------
	double tau=0.;//adimensional
	*/
	if (z - thy_ct >= 0. && z - thy_ct <= thy_c) {//if (z - ct >= 0. && z - ct <= c) {
		double tau = 0.;//added grok also z - thy_ct >= 0 is guaranteed if here!
			if (z - thy_ct <= 0.25 * thy_c) {//if (z - ct >= 0. && z - ct <= 0.25 * c) {
				tau = 1. + ((std::sqrt(2.) - 2.) / 2.) * ((z - thy_ct) / (0.25 * thy_c));
				//tau = 1. + ((std::sqrt(2.) - 2.) / 2.) * ((z - ct) / (0.25 * c));// System.out.println("tau "+tau);
			} else {
				//unnecesary if below
				//if (z - ct >= 0.25 * c && z - ct <= c) {
				//tau = ((-1. + 2. * std::sqrt(2.)) / 3.)	+ ((2. - std::sqrt(2.)) / 2.) * ((z - ct) / (0.75 * c));
				tau = ((-1. + 2. * std::sqrt(2.)) / 3.)  + ((2. - std::sqrt(2.)) / 2.) * ((z - thy_ct) / (0.75 * thy_c));
				//}
			}
			// -----------test------------------------------
			//G4double d = std::pow((x), 2) + std::pow((y - y0), 2);
			const G4double d = x*x + (y - thy_y0)*(y - thy_y0);
			/*if (d <= rr * rr) {
				if (d >= r * r) {
					if (y <= y0) {
						G4double d1 = std::pow(
								(y - y0 - std::abs(x)), 2);
						G4double d2 = 2
								* tau
								* tau
								* (std::pow((x), 2) + std::pow(
										(y - y0), 2));
						if (d1 >= d2)
							return true;
					}
				}
			}*/
			if (d <= thy_rr*thy_rr && d >= thy_r*thy_r && y <= thy_y0) {
				const G4double d1 = (y - thy_y0 - std::abs(x)) * (y - thy_y0 - std::abs(x));
				const G4double d2 = 2. * tau * tau * d;
				if (d1 >= d2) return true;
			}
		}
	return false;
}

bool HumanPhantomSD::inLiver(G4double X, G4double Y, G4double Z) const{
	/*G4double scaleXY=phantom->GetScaleXY();
	G4double scaleZ=phantom->GetScaleZ();*/
	
	//CONCLUSION: xMIRD=x;yMIRD=-z;zMIRD=y;
	const G4double x = X, y = -Z, z = Y;//G4double x=X;G4double y=-Z;G4double z=Y;
	//-----------------------------------
	/*G4String ageGroup=phantom->GetAgeGroup();
	//MirdMapConstants* mmc = new MirdMapConstants();
	//->Gemini static to be reused!!!!! NEW keyword for static is not necessary in cpp!
	static MirdMapConstants mmc;

  //mmc->initTrunk();
  //std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  //G4double bt = scaleXY*trk[ageGroup].Bt;
  //G4double ct = scaleZ*trk[ageGroup].Ct;  
  mmc.initLiver();//mmc->initLiver();
  std::map<std::string,liver_struct> lvr = mmc.GetLiverMap();//mmc->GetLiverMap();  
  G4double a = scaleXY*lvr[ageGroup].a;
  G4double b = scaleXY*lvr[ageGroup].b;
  G4double xm = scaleXY*lvr[ageGroup].xm;
  G4double ym = scaleXY*lvr[ageGroup].ym;  
  G4double zm = scaleZ*lvr[ageGroup].zm;  
  G4double z1 = scaleZ*lvr[ageGroup].z1;  
  G4double z2 = scaleZ*lvr[ageGroup].z2;  
  //delete mmc;
	//-------------------------------
	//G4double a = scaleXY*16.50*cm;
	//G4double b = scaleXY*8.00*cm;
	//G4double xm = scaleXY*35.00*cm;
	//G4double ym = scaleXY*45.00*cm;
	//G4double zm = scaleZ*43.00*cm;
	//G4double z1 = scaleZ*27.00*cm;
	//G4double z2 = scaleZ*43.00*cm;
	*/
	/*if (z >= z1 && z <= z2) {
		double dd = std::pow((x) / a, 2)
				+ std::pow((y) / b, 2);
		if (dd <= 1.) {
			double ddd = (x / xm) + (y / ym)
					- (z / zm);
			if (ddd <= -1.)
				return true;
		}
	}

	return false;*/
	if (z >= liv_z1 && z <= liv_z2) {
    const double dd = (x/liv_a)*(x/liv_a) + (y/liv_b)*(y/liv_b);
    if (dd <= 1.) {
      const double ddd = (x/liv_xm) + (y/liv_ym) - (z/liv_zm);
      if (ddd <= -1.) return true;
    }
  }
  return false;
}

bool HumanPhantomSD::inHeart(G4double X, G4double Y, G4double Z) const{
	/*G4double scaleXY=phantom->GetScaleXY();
	G4double scaleZ=phantom->GetScaleZ();*/
	
	//CONCLUSION: xMIRD=x;yMIRD=-z;zMIRD=y;
	const G4double x = X, y = -Z, z = Y;//G4double x=X;G4double y=-Z;G4double z=Y;
	//--------------------------------
	/*G4String ageGroup=phantom->GetAgeGroup();
	//MirdMapConstants* mmc = new MirdMapConstants();
	//->Gemini static to be reused!!!!! NEW keyword for static is not necessary in cpp!
	static MirdMapConstants mmc;
  //mmc->initTrunk();
  //std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  //G4double bt = scaleXY*trk[ageGroup].Bt;
  //G4double ct = scaleZ*trk[ageGroup].Ct;  
  mmc.initHeart();//mmc->initHeart();
  std::map<std::string,heart_struct> hrt = mmc.GetHeartMap(); //mmc->GetHeartMap();  
  double alfa1 = hrt[ageGroup].alpha1;
  double beta1 = hrt[ageGroup].beta1;
  double gama1 = hrt[ageGroup].gamma1;
  double alfa2 = hrt[ageGroup].alpha2;
  double beta2 = hrt[ageGroup].beta2;
  double gama2 = hrt[ageGroup].gamma2;
  double alfa3 = hrt[ageGroup].alpha3;
  double beta3 = hrt[ageGroup].beta3;
  double gama3 = hrt[ageGroup].gamma3;
  G4double x0 = scaleXY*hrt[ageGroup].x0;
  G4double y0 = scaleXY*hrt[ageGroup].y0;
  G4double z0 = scaleZ*hrt[ageGroup].z0;
  G4double vx = scaleXY*hrt[ageGroup].vx;
  G4double avy = scaleXY*hrt[ageGroup].avy;
  G4double lavz = scaleZ*hrt[ageGroup].lavz;
  G4double ravz = scaleZ*hrt[ageGroup].ravz;
  G4double ax = scaleXY*hrt[ageGroup].ax;
  G4double tlvw =scaleXY*hrt[ageGroup].tlvw;
  G4double trvw =scaleXY*hrt[ageGroup].trvw;
  G4double taw =scaleXY*hrt[ageGroup].tax;  
  //delete mmc;
	//-------------------------------
	*/

	/*G4double x1 = 1.0*cm+alfa1 * (x - x0) + beta1 * (y - y0) + gama1
				* (z - z0);
	G4double y1 = 1.0*cm+alfa2 * (x - x0) + beta2 * (y - y0) + gama2
				* (z - z0);
	G4double z1 = alfa3 * (x - x0) + beta3 * (y - y0) + gama3
				* (z - z0);*/
  const G4double x1 = 1.0*cm + hrt_alpha1*(x - hrt_x0) + hrt_beta1*(y - hrt_y0) + hrt_gamma1*(z - hrt_z0);
  const G4double y1 = 1.0*cm + hrt_alpha2*(x - hrt_x0) + hrt_beta2*(y - hrt_y0) + hrt_gamma2*(z - hrt_z0);
  const G4double z1 =          hrt_alpha3*(x - hrt_x0) + hrt_beta3*(y - hrt_y0) + hrt_gamma3*(z - hrt_z0);
		// left ventricle-WALL+CONTENTS
	if (x1 >= 0.) {
		/*double dd = std::pow((x1) / vx, 2) + std::pow((y1) / avy, 2)
				+ std::pow((z1) / lavz, 2);
		if (dd <= 1.)
			return true;
		*/
		const double dd = (x1/hrt_vx)*(x1/hrt_vx) + (y1/hrt_avy)*(y1/hrt_avy) + (z1/hrt_lavz)*(z1/hrt_lavz);
		if (dd <= 1.) return true;
	}
		// right ventricle -WALL+CONTENTS and dd previous>1---is taken into account
	if (x1 >= 0. && z1 < 0.) {
		/*double dd = std::pow((x1) / vx, 2) + std::pow((y1) / avy, 2)
				+ std::pow((z1) / ravz, 2);
		//portion of left ventricle is excluded so :
		double ddd = std::pow((x1) / vx, 2) + std::pow((y1) / avy, 2)
				+ std::pow((z1) / lavz, 2);//should be >1!!------>redundant coding: it is already returned true!
		if ((dd <= 1.)&&(ddd>1.))
			return true;*/
		const double dd  = (x1/hrt_vx)*(x1/hrt_vx) + (y1/hrt_avy)*(y1/hrt_avy) + (z1/hrt_ravz)*(z1/hrt_ravz);
		const double ddd = (x1/hrt_vx)*(x1/hrt_vx) + (y1/hrt_avy)*(y1/hrt_avy) + (z1/hrt_lavz)*(z1/hrt_lavz);
		if (dd <= 1. && ddd > 1.) return true;
	}
		// left atrium-WALL+CONTENTS
		// part 1
	if (x1 < 0. && z1 >= 0.) {
		/*double dd = std::pow((x1) / ax, 2) + std::pow((y1) / avy, 2)
				+ std::pow((z1) / lavz, 2);
		if (dd <= 1)
			return true;*/
		const double dd = (x1/hrt_ax)*(x1/hrt_ax) + (y1/hrt_avy)*(y1/hrt_avy) + (z1/hrt_lavz)*(z1/hrt_lavz);
		if (dd <= 1.) return true;
	}
		// part2
	if (x1 < 0. && z1 < 0.) {
		/*double dd = std::pow((x1) / ax, 2) + std::pow((y1) / avy, 2)
				+ std::pow((z1) / (lavz - tlvw + taw), 2);
		if (dd <= 1)
			return true;*/
		const double dd = (x1/hrt_ax)*(x1/hrt_ax) + (y1/hrt_avy)*(y1/hrt_avy)
                    + (z1/(hrt_lavz - hrt_tlvw + hrt_taw))*(z1/(hrt_lavz - hrt_tlvw + hrt_taw));
		if (dd <= 1.) return true;
	}
	 // right atrium-WALL+CONTENTS dd previous>1
	if (x1 < 0. && z1 < 0.) {
		/*double dd = std::pow((x1) / ax, 2) + std::pow((y1) / avy, 2)
				+ std::pow((z1) / ravz, 2);
		//exclusion of common volume (with left atrium)
		double ddd = std::pow((x1) / ax, 2) + std::pow((y1) / avy, 2)
				+ std::pow((z1) / (lavz-tlvw+taw), 2);//should be >1!!------>redundant coding: it is already returned true!

		if ((dd <= 1.)&&(ddd>1.))
			return true;*/
		const double dd  = (x1/hrt_ax)*(x1/hrt_ax) + (y1/hrt_avy)*(y1/hrt_avy) + (z1/hrt_ravz)*(z1/hrt_ravz);
		const double ddd = (x1/hrt_ax)*(x1/hrt_ax) + (y1/hrt_avy)*(y1/hrt_avy)
                     + (z1/(hrt_lavz - hrt_tlvw + hrt_taw))*(z1/(hrt_lavz - hrt_tlvw + hrt_taw));
		if (dd <= 1. && ddd > 1.) return true;
	}
	
	return false;
}

//bool HumanPhantomSD::inLeftClavicle(G4double X, G4double Y, G4double Z) const{
bool HumanPhantomSD::inLeftClavicle(G4double X, G4double Y, G4double Z) const{
	/*G4double scaleXY=phantom->GetScaleXY();
	G4double scaleZ=phantom->GetScaleZ();*/
	
	//CONCLUSION: xMIRD=x;yMIRD=-z;zMIRD=y;
	const G4double x = X, y = -Z, z = Y;//G4double x=X;G4double y=-Z;G4double z=Y;
	//-----------------------------
	/*G4String ageGroup=phantom->GetAgeGroup();
	//MirdMapConstants* mmc = new MirdMapConstants();
	//->Gemini static to be reused!!!!! NEW keyword for static is not necessary in cpp!
	static MirdMapConstants mmc;

  //mmc->initTrunk();
  //std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  //G4double bt = scaleXY*trk[ageGroup].Bt;
  //G4double ct = scaleZ*trk[ageGroup].Ct;  
  mmc.initClavicles();//mmc->initClavicles();
  std::map<std::string,clavicles_struct> clv = mmc.GetClaviclesMap(); //mmc->GetClaviclesMap();  
  double cot1 = clv[ageGroup].cotan1;
  double cot2 = clv[ageGroup].cotan2;
  G4double y0 = scaleXY*clv[ageGroup].y0;
  G4double z1 = scaleZ*clv[ageGroup].z1;
  G4double rr = scaleXY*clv[ageGroup].R;  
  G4double r = scaleZ*clv[ageGroup].r;//Z or XY??????
  //delete mmc;
	//-------------------------------
	*//*G4double y0 = scaleXY*11.10*cm;
	G4double z1 = scaleZ*68.25*cm;
	G4double rr = scaleXY*20.00*cm;////////////////////
	G4double r = scaleXY*0.7883*cm;//scaleZ*0.7883*cm;/////////////////////torus=>see lowerLargeIntestine
	//G4double r = 0.7883*cm;/////////////////////
	double cot1 = 7.0342;
	double cot2 = 0.89415;*/
		// -------------------

	if (y < 0.) {
		/*double cot = (y0 - y) / x;//LEFT//Math.abs(coord[0]);
		if (cot >= cot2 && cot <= cot1) {
			G4double dd = std::pow((z - z1), 2)
					+ std::pow(
							(rr - std::sqrt(x * x
									+ (y - y0) * (y - y0))),
							2);
			if (dd <= r * r)
				return true;
		}*/
		const double cot = (clv_y0 - y) / x;
		if (cot >= clv_cot2 && cot <= clv_cot1) {
			const G4double rad = std::sqrt(x*x + (y - clv_y0)*(y - clv_y0));
			const G4double dd  = (z - clv_z1)*(z - clv_z1) + (clv_rr - rad)*(clv_rr - rad);
			if (dd <= clv_r * clv_r) return true;
		}
	}

	return false;
}

bool HumanPhantomSD::inRightClavicle(G4double X, G4double Y, G4double Z)const{
	/*G4double scaleXY=phantom->GetScaleXY();
	G4double scaleZ=phantom->GetScaleZ();*/
	
	//CONCLUSION: xMIRD=x;yMIRD=-z;zMIRD=y;
	const G4double x = X, y = -Z, z = Y;//G4double x=X;G4double y=-Z;G4double z=Y;
	//-----------------------
	/*G4String ageGroup=phantom->GetAgeGroup();
	//MirdMapConstants* mmc = new MirdMapConstants();
	//->Gemini static to be reused!!!!! NEW keyword for static is not necessary in cpp!
	static MirdMapConstants mmc;

  //mmc->initTrunk();
  //std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  //G4double bt = scaleXY*trk[ageGroup].Bt;
  //G4double ct = scaleZ*trk[ageGroup].Ct;  
  mmc.initClavicles();//mmc->initClavicles();
  std::map<std::string,clavicles_struct> clv = mmc.GetClaviclesMap();//mmc->GetClaviclesMap();  
  double cot1 = clv[ageGroup].cotan1;
  double cot2 = clv[ageGroup].cotan2;
  G4double y0 = scaleXY*clv[ageGroup].y0;
  G4double z1 = scaleZ*clv[ageGroup].z1;
  G4double rr = scaleXY*clv[ageGroup].R;  
  G4double r = scaleZ*clv[ageGroup].r;//Z or XY??????
  //delete mmc;
	//-------------------------------
	*//*G4double y0 = scaleXY*11.10*cm;
	G4double z1 = scaleZ*68.25*cm;
	G4double rr = scaleXY*20.00*cm;///////////////////////////////////////
	G4double r = scaleXY*0.7883*cm;//scaleZ*0.7883*cm;///////////////////////////////////////
	double cot1 = 7.0342;
	double cot2 = 0.89415;*/
		// -------------------

	if (y < 0.) {
		/*double cot = (y0 - y) /(-x);//RIGHT//Math.abs(coord[0]);
		if (cot >= cot2 && cot <= cot1) {
			G4double dd = std::pow((z - z1), 2)
					+ std::pow(
							(rr - std::sqrt(x * x
									+ (y - y0) * (y - y0))),
							2);
			if (dd <= r * r)
				return true;
		}*/
		const double cot = (clv_y0 - y) / (-x);
		if (cot >= clv_cot2 && cot <= clv_cot1) {
			const G4double rad = std::sqrt(x*x + (y - clv_y0)*(y - clv_y0));
			const G4double dd  = (z - clv_z1)*(z - clv_z1) + (clv_rr - rad)*(clv_rr - rad);
			if (dd <= clv_r * clv_r) return true;
		}
	}

	return false;
}

bool HumanPhantomSD::inGallBladder(G4double X, G4double Y, G4double Z) const{
	/*G4double scaleXY=phantom->GetScaleXY();
	G4double scaleZ=phantom->GetScaleZ();*/
	
	//CONCLUSION: xMIRD=x;yMIRD=-z;zMIRD=y;
	const G4double x = X, y = -Z, z = Y;//G4double x=X;G4double y=-Z;G4double z=Y;
	//----------------------
	/*G4String ageGroup=phantom->GetAgeGroup();
	//MirdMapConstants* mmc = new MirdMapConstants();
	//->Gemini static to be reused!!!!! NEW keyword for static is not necessary in cpp!
	static MirdMapConstants mmc;
  //mmc->initTrunk();
  //std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  //G4double bt = scaleXY*trk[ageGroup].Bt;
  //G4double ct = scaleZ*trk[ageGroup].Ct;  
  
  mmc.initGallBladder();//mmc->initGallBladder();
  std::map<std::string,gallBladder_struct> glb = mmc.GetGallBladderMap();//mmc->GetGallBladderMap();  
  double alfa1 = glb[ageGroup].alpha1;
  double beta1 = glb[ageGroup].beta1;
  double gama1 = glb[ageGroup].gamma1;
  double alfa2 = glb[ageGroup].alpha2;
  double beta2 = glb[ageGroup].beta2;
  double gama2 = glb[ageGroup].gamma2;
  double alfa3 = glb[ageGroup].alpha3;
  double beta3 = glb[ageGroup].beta3;
  double gama3 = glb[ageGroup].gamma3;
  G4double r1 = scaleZ*glb[ageGroup].r1;
  G4double r2 = scaleXY*glb[ageGroup].r2;
  double s = glb[ageGroup].s;
  G4double h = scaleZ*glb[ageGroup].h;  
  G4double x0 = scaleXY*glb[ageGroup].x0;
  G4double y0 = scaleXY*glb[ageGroup].y0;
  G4double z0 = scaleZ*glb[ageGroup].z0;  
  //delete mmc;
	//-------------------------------
	*//*double alfa1 = 0.9615;
	double beta1 = 0.0;
	double gama1 = -0.2748;

	double alfa2 = -0.0574;
	double beta2 = 0.9779;
	double gama2 = -0.2008;

	double alfa3 = 0.2687;
	double beta3 = 0.2090;
	double gama3 = 0.9403;
	
	//G4double r1 = 2.000*cm;
	G4double r2 = scaleXY*2.120*cm;
	double s = 0.2275;
	G4double h = scaleZ*8.00*cm;
	G4double x0 = -scaleXY*4.50*cm;
	G4double y0 = -scaleXY*3.20*cm;
	G4double z0 = scaleZ*30.00*cm;*/
		// ------------------------

	/*G4double x1 = alfa1 * (x - x0) + beta1 * (y - y0) + gama1
				* (z - z0);
	G4double y1 = alfa2 * (x - x0) + beta2 * (y - y0) + gama2
				* (z - z0);
	G4double z1 = alfa3 * (x - x0) + beta3 * (y - y0) + gama3
				* (z - z0);*/

	const G4double x1 = glb_alpha1*(x - glb_x0) + glb_beta1*(y - glb_y0) + glb_gamma1*(z - glb_z0);
	const G4double y1 = glb_alpha2*(x - glb_x0) + glb_beta2*(y - glb_y0) + glb_gamma2*(z - glb_z0);
	const G4double z1 = glb_alpha3*(x - glb_x0) + glb_beta3*(y - glb_y0) + glb_gamma3*(z - glb_z0);
	
	// hemispherical part
	if (z1 < 0.) {
		/*G4double dd = std::pow((x1), 2) + std::pow((y1), 2)
					+ std::pow((z1), 2);
		if (dd <= r2 * r2)
				return true;*/
		if (x1*x1 + y1*y1 + z1*z1 <= glb_r2*glb_r2) return true;
	}
	// conical part
	if (z1 >= 0. && z1 <= glb_h) {//if (z1 >= 0. && z1 <= h) {//a height
		/*G4double dd = std::pow((x1), 2) + std::pow((y1), 2);//a surface
		G4double ddd = std::pow((r2 - s * z1), 2);//+wall!
		//G4double ddd = std::pow((r1 - s * z1), 2);//contents
		if (dd <= ddd)
			return true;*/
		const G4double dd  = x1*x1 + y1*y1;
		const G4double ddd = (glb_r2 - glb_s * z1) * (glb_r2 - glb_s * z1);
		if (dd <= ddd) return true;
	}
	
	return false;	
}

bool HumanPhantomSD::inSmallIntestine(G4double X, G4double Y, G4double Z) const{
	/*G4double scaleXY=phantom->GetScaleXY();
	G4double scaleZ=phantom->GetScaleZ();*/
	
    //CONCLUSION: xMIRD=x;yMIRD=-z;zMIRD=y;
	const G4double x = X, y = -Z, z = Y;//G4double x=X;G4double y=-Z;G4double z=Y;
	//-----------------------------------
	/*G4String ageGroup=phantom->GetAgeGroup();
	//MirdMapConstants* mmc = new MirdMapConstants();
	//->Gemini static to be reused!!!!! NEW keyword for static is not necessary in cpp!
	static MirdMapConstants mmc;
  //mmc->initTrunk();
  //std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  //G4double bt = scaleXY*trk[ageGroup].Bt;
  //G4double ct = scaleZ*trk[ageGroup].Ct;  
  mmc.initSmallIntestine();//mmc->initSmallIntestine();
  std::map<std::string,smallIntestine_struct> smi = mmc.GetSmallIntestineMap();//mmc->GetSmallIntestineMap();  
  G4double a = scaleXY*smi[ageGroup].a;
  G4double b = scaleXY*smi[ageGroup].b;
  G4double y0 = scaleXY*smi[ageGroup].y0;  
  G4double y1 = scaleXY*smi[ageGroup].y1;  
  G4double y2 = scaleXY*smi[ageGroup].y2;  
  G4double z1 = scaleZ*smi[ageGroup].z1;  
  G4double z2 = scaleZ*smi[ageGroup].z2;  
  //delete mmc;
	//-------------------------------
	*//*G4double a = scaleXY*11.30*cm;
	G4double b = scaleXY*11.30*cm;
	G4double y0 = -scaleXY*3.80*cm;
	G4double y1 = -scaleXY*4.86*cm;
	G4double y2 = scaleXY*2.20*cm;
	G4double z1 = scaleZ*17.00*cm;
	G4double z2 = scaleZ*27.00*cm;*/

	/*if (y >= y1 && y <= y2) {
		if (z >= z1 && z <= z2) {
			double dd = std::pow((x) / a, 2)
					+ std::pow((y - y0) / b, 2);
			if (dd <= 1.)
				return true;
		}
	}*/
	if (y >= smi_y1 && y <= smi_y2 && z >= smi_z1 && z <= smi_z2) {
		const double dd = (x/smi_a)*(x/smi_a) + ((y - smi_y0)/smi_b)*((y - smi_y0)/smi_b);
		if (dd <= 1.) return true;
	}

	return false;
}

bool HumanPhantomSD::inThymus(G4double X, G4double Y, G4double Z) const{
	/*G4double scaleXY=phantom->GetScaleXY();
	G4double scaleZ=phantom->GetScaleZ();*/
	
	//CONCLUSION: xMIRD=x;yMIRD=-z;zMIRD=y;
	const G4double x = X, y = -Z, z = Y;//G4double x=X;G4double y=-Z;G4double z=Y;
	//---------------------------
	/*G4String ageGroup=phantom->GetAgeGroup();
	//MirdMapConstants* mmc = new MirdMapConstants();
	//->Gemini static to be reused!!!!! NEW keyword for static is not necessary in cpp!
	static MirdMapConstants mmc;
  //mmc->initTrunk();
  //std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  //G4double bt = scaleXY*trk[ageGroup].Bt;
  //G4double ct = scaleZ*trk[ageGroup].Ct;  
  mmc.initThymus();//mmc->initThymus();
  std::map<std::string,thymus_struct> thm = mmc.GetThymusMap();//mmc->GetThymusMap();  
  G4double a = scaleXY*thm[ageGroup].a;
  G4double b = scaleXY*thm[ageGroup].b;
  G4double c = scaleZ*thm[ageGroup].c;
  G4double y0 = scaleXY*thm[ageGroup].y0;
  G4double z0 = scaleZ*thm[ageGroup].z0;
  //delete mmc;
	//-------------------------------
	//G4double a = scaleXY*1.50*cm;
	//G4double b = scaleXY*0.80*cm;
	//G4double c = scaleZ*4.00*cm;
	//G4double y0 = -scaleXY*7.30*cm;
	//G4double z0 = scaleZ*57.00*cm;
*/
	/*double dd = std::pow((x) / a, 2)
			+ std::pow((y - y0) / b, 2)
			+ std::pow((z - z0) / c, 2);
	if (dd <= 1.)
		return true;

	return false;*/
	const double dd = (x/thm_a)*(x/thm_a)
                  + ((y - thm_y0)/thm_b)*((y - thm_y0)/thm_b)
                  + ((z - thm_z0)/thm_c)*((z - thm_z0)/thm_c);
	return (dd <= 1.);
}

bool HumanPhantomSD::inLeftTesticle(G4double X, G4double Y, G4double Z) const{
	/*G4double scaleXY=phantom->GetScaleXY();
	G4double scaleZ=phantom->GetScaleZ();*/
	//G4cout << " LEFT TESTICLE " << G4endl;
	//CONCLUSION: xMIRD=x;yMIRD=-z;zMIRD=y;
	const G4double x = X, y = -Z, z = Y;//G4double x=X;G4double y=-Z;G4double z=Y;
	//-----------------------
	/*G4String ageGroup=phantom->GetAgeGroup();
	//MirdMapConstants* mmc = new MirdMapConstants();
	//->Gemini static to be reused!!!!! NEW keyword for static is not necessary in cpp!
	static MirdMapConstants mmc;

  //mmc->initTrunk();
  //std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  //G4double bt = scaleXY*trk[ageGroup].Bt;
  //G4double ct = scaleZ*trk[ageGroup].Ct;  
  mmc.initTestes();//mmc->initTestes();
  std::map<std::string,testes_struct> tst = mmc.GetTestesMap(); //mmc->GetTestesMap();  
  G4double a = scaleXY*tst[ageGroup].a;
  G4double b = scaleXY*tst[ageGroup].b;
  G4double c = scaleZ*tst[ageGroup].c;
  G4double y0 = scaleXY*tst[ageGroup].y0;  
  //delete mmc;
	//-------------------------------
	//G4double a = scaleXY*1.30*cm;
	//G4double b = scaleXY*1.50*cm;
	//G4double c = scaleZ*2.30*cm;
	//G4double y0 = -scaleXY*8.00*cm;
	*/

	/*
	double d = std::pow((x - a) / a, 2)
			+ std::pow((y - y0) / b, 2)
			+ std::pow((z + c) / c, 2);// left testes	
	
	if (d <= 1.)
		return true;

	return false;*/
	const double d = ((x - tst_a)/tst_a)*((x - tst_a)/tst_a)
                 + ((y - tst_y0)/tst_b)*((y - tst_y0)/tst_b)
                 + ((z + tst_c)/tst_c)*((z + tst_c)/tst_c);
	return (d <= 1.);
}

bool HumanPhantomSD::inRightTesticle(G4double X, G4double Y, G4double Z) const{
	/*G4double scaleXY=phantom->GetScaleXY();
	G4double scaleZ=phantom->GetScaleZ();*/
	
	//CONCLUSION: xMIRD=x;yMIRD=-z;zMIRD=y;
	const G4double x = X, y = -Z, z = Y;//G4double x=X;G4double y=-Z;G4double z=Y;
	//-------------------------------
	/*G4String ageGroup=phantom->GetAgeGroup();
	//MirdMapConstants* mmc = new MirdMapConstants();
	//->Gemini static to be reused!!!!! NEW keyword for static is not necessary in cpp!
	static MirdMapConstants mmc;
  //mmc->initTrunk();
  //std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  //G4double bt = scaleXY*trk[ageGroup].Bt;
  //G4double ct = scaleZ*trk[ageGroup].Ct;  
  mmc.initTestes();//mmc->initTestes();
  std::map<std::string,testes_struct> tst = mmc.GetTestesMap();//mmc->GetTestesMap();  
  G4double a = scaleXY*tst[ageGroup].a;
  G4double b = scaleXY*tst[ageGroup].b;
  G4double c = scaleZ*tst[ageGroup].c;
  G4double y0 = scaleXY*tst[ageGroup].y0;  
  //delete mmc;
	//-------------------------------
	//G4double a = scaleXY*1.30*cm;
	//G4double b = scaleXY*1.50*cm;
	//G4double c = scaleZ*2.30*cm;
	//G4double y0 = -scaleXY*8.00*cm;
	*/
	/*
	double d1 = std::pow((x + a) / a, 2)
			+ std::pow((y - y0) / b, 2)
			+ std::pow((z + c) / c, 2);// right testes
	//if (y<-70*cm) return true;///yes in air!!!
	if (d1 <= 1.)
		return true;

	return false;*/
	 const double d1 = ((x + tst_a)/tst_a)*((x + tst_a)/tst_a)
                  + ((y - tst_y0)/tst_b)*((y - tst_y0)/tst_b)
                  + ((z + tst_c)/tst_c)*((z + tst_c)/tst_c);
	return (d1 <= 1.);
}

bool HumanPhantomSD::inEsophagus(G4double X, G4double Y, G4double Z) const{
	/*G4double scaleXY=phantom->GetScaleXY();
	G4double scaleZ=phantom->GetScaleZ();*/
	//CONCLUSION: xMIRD=x;yMIRD=-z;zMIRD=y;
	 const G4double x = X, y = -Z, z = Y;//G4double x=X;G4double y=-Z;G4double z=Y;
	//--------------------------------
	/*G4String ageGroup=phantom->GetAgeGroup();
	//MirdMapConstants* mmc = new MirdMapConstants();
	//->Gemini static to be reused!!!!! NEW keyword for static is not necessary in cpp!
	static MirdMapConstants mmc;

  //mmc->initTrunk();
  //std::map<std::string,trunk_struct> trk = mmc->GetTrunkMap();  
  //G4double bt = scaleXY*trk[ageGroup].Bt;
  //G4double ct = scaleZ*trk[ageGroup].Ct;  
  mmc.initEsophagus();//mmc->initEsophagus();
  std::map<std::string,esophagus_struct> esp = mmc.GetEsophagusMap();//mmc->GetEsophagusMap();  
  double alfa1 = esp[ageGroup].alpha1;
  double beta1 = esp[ageGroup].beta1;
  double gama1 = esp[ageGroup].gamma1;
  double alfa2 = esp[ageGroup].alpha2;
  double beta2 = esp[ageGroup].beta2;
  double gama2 = esp[ageGroup].gamma2;
  double alfa3 = esp[ageGroup].alpha3;
  double beta3 = esp[ageGroup].beta3;
  double gama3 = esp[ageGroup].gamma3;
  G4double a = scaleXY*esp[ageGroup].a;
  G4double y0 = scaleXY*esp[ageGroup].y0;
  G4double b = scaleXY*esp[ageGroup].b;
  G4double d = scaleXY*esp[ageGroup].d;  
  G4double z1 = scaleZ*esp[ageGroup].z1;
  G4double z2 = scaleZ*esp[ageGroup].z2;
  G4double z3 = scaleZ*esp[ageGroup].z3;
  G4double r = scaleXY*esp[ageGroup].r;
  G4double xp1 =scaleXY*esp[ageGroup].x1p;
  G4double xp2 =scaleXY*esp[ageGroup].x2p;
  //delete mmc;
	//-------------------------------
	*//*double alfa1 = 0.736084;
	double beta1 = -0.604969;
	double gama1 = -0.303634;

	double alfa2 = 0.634945;
	double beta2 = 0.772557;
	double gama2 = 0.0;

	double alfa3 = 0.234575;
	double beta3 = -0.192791;
	double gama3 = 0.952789;
	
	G4double a = scaleXY*1.17*cm;
	G4double b = scaleXY*0.42*cm;
	G4double d = scaleXY*0.30*cm;
	G4double y0 = scaleXY*2.575*cm;
	G4double z2 = scaleZ*43.00*cm;
	G4double z3 = scaleZ*70.00*cm;
	G4double r = scaleXY*0.70*cm;//scaleZ*0.70*cm;///////////////////////////
	G4double xp1 = scaleZ*0.10*cm;//scaleXY*0.10*cm;//x'1
	G4double xp2 = scaleZ*7.80*cm;//scaleXY*7.80*cm;//////////////////////
	G4double z1 = scaleZ*42.30*cm;*/
		// ------------------------

	/*G4double xp = alfa1 * (x) + beta1 * (y - y0) + gama1
				* (z - z1);
	G4double yp = alfa2 * (x) + beta2 * (y - y0) + gama2
				* (z - z1);
	G4double zp = alfa3 * (x) + beta3 * (y - y0) + gama3
				* (z - z1);*/
	const G4double xp = esp_alpha1*x + esp_beta1*(y - esp_y0) + esp_gamma1*(z - esp_z1);
	const G4double yp = esp_alpha2*x + esp_beta2*(y - esp_y0) + esp_gamma2*(z - esp_z1);
	const G4double zp = esp_alpha3*x + esp_beta3*(y - esp_y0) + esp_gamma3*(z - esp_z1);

	//==============================thoracic part:
	if (z >= esp_z2 && z <= esp_z3) {//if (z >= z2 && z <= z3) {
		/*double dd = std::pow(x/a, 2) + std::pow((y-y0)/b, 2);
		double dd1 = std::pow(x/(a-d), 2) + std::pow((y-y0)/(b-d), 2);
		if(dd<=1. && dd1>=1.)
			return true;*/
		const double dd  = (x/esp_a)*(x/esp_a) + ((y - esp_y0)/esp_b)*((y - esp_y0)/esp_b);
		const double dd1 = (x/(esp_a - esp_d))*(x/(esp_a - esp_d))
                     + ((y - esp_y0)/(esp_b - esp_d))*((y - esp_y0)/(esp_b - esp_d));
		if (dd <= 1. && dd1 >= 1.) return true;
	}
	//===============================Abdominal part
	if (xp >= esp_xp1 && xp <= esp_xp2) {//if (xp >= xp1 && xp <= xp2) {//height=>scaleZ
		/*G4double dd = std::pow((yp), 2) + std::pow((zp), 2);//surface=>xy scale
		G4double ddd = std::pow((r), 2);
		if (dd <= ddd)
			return true;*/
		if (yp*yp + zp*zp <= esp_r*esp_r) return true;
	}
	
	return false;	
}

bool HumanPhantomSD::inFacial(G4double X, G4double Y, G4double Z) const{
	/*G4double scaleXY=phantom->GetScaleXY();
	G4double scaleZ=phantom->GetScaleZ();*/
	//CONCLUSION: xMIRD=x;yMIRD=-z;zMIRD=y;
	const G4double x = X, y = -Z, z = Y;//G4double x=X;G4double y=-Z;G4double z=Y;
	//--------------------------------
	/*G4String ageGroup=phantom->GetAgeGroup();
	//MirdMapConstants* mmc = new MirdMapConstants();
	//->Gemini static to be reused!!!!! NEW keyword for static is not necessary in cpp!
	static MirdMapConstants mmc;

  mmc.initTrunk();//mmc->initTrunk();
  std::map<std::string,trunk_struct> trk = mmc.GetTrunkMap();//mmc->GetTrunkMap();    
  G4double ct = scaleZ*trk[ageGroup].Ct;  
  mmc.initHead();//mmc->initHead();
  std::map<std::string,head_struct> hed = mmc.GetHeadMap();//mmc->GetHeadMap();    
  G4double ch1 = scaleZ*hed[ageGroup].Ch1; 
  mmc.initBrain();//mmc->initBrain();
  std::map<std::string,brain_struct> brn = mmc.GetBrainMap();//mmc->GetBrainMap();    
  G4double aS = scaleXY*brn[ageGroup].a;
  G4double bS = scaleXY*brn[ageGroup].b;
  G4double cS = scaleZ*brn[ageGroup].c;
  mmc.initSkull();//mmc->initSkull();
  std::map<std::string,skull_struct> skl = mmc.GetSkullMap();//mmc->GetSkullMap();    
  G4double dS = scaleXY*skl[ageGroup].d; 
  mmc.initFacial();//mmc->initFacial();
  std::map<std::string,facial_struct> fcl = mmc.GetFacialMap(); //mmc->GetFacialMap();    
  G4double a1 = scaleXY*fcl[ageGroup].a1;
  G4double b1 = scaleXY*fcl[ageGroup].b1;  
  G4double d = scaleXY*fcl[ageGroup].d;  
  G4double z1 = scaleZ*fcl[ageGroup].z1;
  G4double z5 = scaleZ*fcl[ageGroup].z5;  
  //delete mmc;*/
	/*
	G4double a2 = aS+dS;
	G4double b2 = bS+dS;
	G4double c2 = cS+dS;
	*/
	const G4double a2 = fac_aS + fac_dS;
	const G4double b2 = fac_bS + fac_dS;
	const G4double c2 = fac_cS + fac_dS;
	
	/*if (y<=0.)
		if (ct+z1<=z && z<=ct+z5){
			double zzz = ct+ch1;
			double dd = std::pow(x/a2, 2) + std::pow(y/b2, 2)+ std::pow((z-zzz)/c2, 2);
			if (dd>1.){//not in corresponding skull/brain!!
				double dd1 = std::pow(x/a1, 2) + std::pow(y/b1, 2);
				double dd2 = std::pow(x/(a1-d), 2) + std::pow(y/(b1-d), 2);
				if (dd1<=1. && dd2>=1.)
					return true;
			}
		}		
	//}
		*/
	if (y <= 0. && fac_ct + fac_z1 <= z && z <= fac_ct + fac_z5) {
		const double zzz = fac_ct + fac_ch1;
		const double dd  = (x/a2)*(x/a2) + (y/b2)*(y/b2) + ((z - zzz)/c2)*((z - zzz)/c2);
		if (dd > 1.) {   // outside corresponding skull/brain
			const double dd1 = (x/fac_a1)*(x/fac_a1) + (y/fac_b1)*(y/fac_b1);
			const double dd2 = (x/(fac_a1 - fac_d))*(x/(fac_a1 - fac_d))
                       + (y/(fac_b1 - fac_d))*(y/(fac_b1 - fac_d));
		if (dd1 <= 1. && dd2 >= 1.) return true;
		}
	}
	return false;	
}