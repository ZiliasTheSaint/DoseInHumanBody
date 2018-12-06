//D. Fulea, National Institute of Public Health Bucharest, Cluj-Napoca Regional Center, Romania
//#include"XRayBuilder.hh"
#include"XRAYBuilder.hh"

string XRayBuilder::filename = " ";
const string XRayBuilder::dataS="Data";
const string XRayBuilder::data_filterS = "FILTERS";
const string XRayBuilder::data_kermaS = "KERMA";
const string XRayBuilder::data_wS = "SPECW";
const string XRayBuilder::data_mS = "SPECMO";
const string XRayBuilder::data_rS = "SPECRH";
const string XRayBuilder::data_constantS = "CONSTANT";
const string XRayBuilder::data_rippleS = "RIPPLE";
const string XRayBuilder::defaultSpectrumExt = ".SPC";
const string XRayBuilder::filterKermaExt = ".CSV";
const string XRayBuilder::rippleExt_R0 = ".R0";
const string XRayBuilder::rippleExt_R = ".R";
#ifdef _WIN32
#define FILE_SEPARATOR   "\\"
#else
#define FILE_SEPARATOR   "/"
#endif

XRayBuilder::XRayBuilder(double kv, double filtration, double anodAngle, double ripple, string anodMaterial)
	:ndata(0)
{
	this->kv=kv;
	this->filtration=filtration;
	this->anodAngle=anodAngle;
	this->ripple=ripple;
	this->anodMaterial=anodMaterial;

	if (anodMaterial=="W") ianod=0;
	else if (anodMaterial=="MO") ianod=1;
	else if (anodMaterial=="RH") ianod=2;

	//-------------
	MAXITER = 10000;//maximum iteration used in HVL evaluation
	THI_HVL1 = 120.0;
	TMM=new double[20];
	NOATT=0;
	//----------
	generateSpectrum();

	G4cout<<"\n" <<"Spectrum generated for: " << G4endl
	<<"kVp: "<<kv<<"; filtration: "<<filtration<< "; anode angle: "
		<<anodAngle<<"; ripple: "<<ripple << "; anode material: "<<anodMaterial << G4endl;
	G4cout<<"\n-------------------------------------"<< G4endl;
}

XRayBuilder::~XRayBuilder()
{;}

void XRayBuilder::generateSpectrum(){

	int kvint=kv;//make sure is an int
	string ikv = intToString(kvint);
	string kvS = "";
	if (kv < 100.0)
		kvS = "0" + ikv;
	else
		kvS = ikv;
	
	int ianod_file = anodAngle;
	string ianod_fileS=intToString(ianod_file);
	string uanodS = "";
	if (anodAngle >= 10) {
		if (ianod == 0)
			uanodS = ianod_fileS + "0";
		else if (ianod == 1)
			uanodS = ianod_fileS + "1";
		else if (ianod == 2)
			uanodS = ianod_fileS + "2";
	} else {
		if (ianod == 0)
		   uanodS = "0" + ianod_fileS +"0";
		else if (ianod == 1)
		   uanodS = "0" + ianod_fileS +"1";
		else if (ianod == 2)
		   uanodS = "0" + ianod_fileS + "2";
	}
	
	string filename = kvS + uanodS;	
	readSpectrum(filename);
	
	filename = "KERMAIR";// allways
	readKerma(filename);	

	filename="AL";//allways
	readAttCoef(filename,1);//allways
	TMM[0]=filtration;//allways
		
	//G4cout<<"\n NDATAAAAAAAAAAAAAA= "<<ndata<<" EMAX "<<en_array[ndata-1]
	//<<" ; YS "<<YS[ndata-1]<<" ; KERMAFactor "<<KVAL[ndata-1]<<" ; ATT[n-1][0] "<<ATT[ndata-1][0]<<" \n";

	buildSpectra();
	// in order to compute kermain air=>hvlcal
	computeHVL1("AL", false);  
}

void XRayBuilder::readSpectrum(string filename){
	string filenam = "";
	string ext = "";
	string subdir = "";
	string dir = "";

	string rippleS=intToString(ripple);

	if (ripple == 0) {
		ext = defaultSpectrumExt;
		subdir = data_constantS;
	} else if (ripple < 10) {
		ext = rippleExt_R0 + rippleS;
		subdir = data_rippleS;
	} else {
		ext = rippleExt_R + rippleS;
		subdir = data_rippleS;
	}
	
	if (ianod == 0) {
		dir = data_wS;
	} else if (ianod == 1) {
		dir = data_mS;
	} else if (ianod == 2) {
		dir = data_rS;
	}
	
	string file_sep=FILE_SEPARATOR;
	filenam = dataS + file_sep + dir + file_sep + subdir + file_sep
				+ filename + ext;

	//G4cout<<"\n filepath= " <<filenam;
	int number_of_lines = 0;
    string line;
    ifstream myfile(filenam.c_str());//filenam);//THIS IS REQUIRED FOR LINUX C COMPILER!!!!

	if (!myfile.is_open()){
		  G4cout<<"Error opening file SPECTRUM!!!!!!!!!!!!!!!! \n";
		  return;
	}

	ndata=0;//initialize ndata
	int dataCounter=0;//in line data counter
	double temp;
	vector<double> enV;
	vector<double> pV;

	while (myfile.good()){
		number_of_lines++;
		getline(myfile, line);//read line and put it in line variable.
        
		char* str=new char[line.length()+1];
		strcpy(str,line.c_str());//copy line into str.
		//---------------------------
		if (number_of_lines==1){
			//G4cout<<"\n"<<line<<"\n";//print some info to know what file we just opened!
		}
		//-----------------------------
		if (number_of_lines>1){
			char* pch;
			pch=strtok(str," ,");//space and comma

			while (pch != NULL)
			{
				dataCounter++;
				if(dataCounter==1){
					istringstream buffer(pch);
					buffer>>temp;
					enV.push_back(temp);
				}
				if(dataCounter==2){
					istringstream buffer(pch);
					buffer>>temp;
					pV.push_back(temp);
					ndata=ndata+1;
					dataCounter=0;//reset data counter
				}
				pch = strtok (NULL, " ,");//next data and also break the loop if next data = null
			}
			//G4cout <<enV.back()<<"  "<<pV.back()<<"\n";
		}
		
	}//while (myfile.good())
	
	//G4cout << "\nNdata= " << ndata<<"\n";
	//G4cout << "\nNumber of lines in text file [SPECTRUM]: " << number_of_lines<<"\n";

	en_array=new double[ndata];
	YS=new double[ndata];
	XVAL=new double[ndata];
	ATT=new double*[ndata];
	for (int i=0;i<ndata;i++){
		en_array[i]=enV[i];//Energy in keV
		YS[i]=pV[i];//number of photons/(mAs*mm2) at 750 mm distance from X-Ray tube. Unattenuated spectra!
		XVAL[i]=pV[i];

			//prepare attenuation matrix====
	
		//for(int i=0;i<ndata;i++)
		ATT[i]=new double[20];
	}

	DE=en_array[1]-en_array[0];//delta energy=width of a bin

	myfile.close();
	//G4cout << "\nLast 2 energies: " << en_array[ndata-2]<<"  "<<en_array[ndata-1]<<"\n";
	//--------------------	 
}

void XRayBuilder::readKerma(string filename){
	string filenam = "";	
	string file_sep=FILE_SEPARATOR;	
	filenam = dataS + file_sep + data_kermaS + file_sep + filename
				+ filterKermaExt;
	//G4cout<<"\n filepath= " <<filenam;
	int number_of_lines = 0;
    string line;
    ifstream myfile(filenam.c_str());//filenam);

	if (!myfile.is_open()){
		  G4cout<<"Error opening file KERMA!!!!!!!!!!!!!!!! \n";
		  return;
	}

	int ndatak=0;
	int dataCounter=0;//in line data counter
	double temp;
	vector<double> enV;//same range as readSpectra
	vector<double> pV;

	while (myfile.good()){
		number_of_lines++;
		getline(myfile, line);//read line and put it in line variable.
        
		char* str=new char[line.length()+1];
		strcpy(str,line.c_str());//copy line into str.

		if (number_of_lines>0){
			char* pch;
			pch=strtok(str," ,");//space and comma

			while (pch != NULL)
			{
				dataCounter++;
				if(dataCounter==1){
					istringstream buffer(pch);
					buffer>>temp;
					enV.push_back(temp);
				}
				if(dataCounter==2){
					istringstream buffer(pch);
					buffer>>temp;
					pV.push_back(temp);
					ndatak++;
					dataCounter=0;//reset data counter
				}
				pch = strtok (NULL, " ,");//next data and also break the loop if next data = null
			}
			//G4cout <<enV.back()<<"  "<<pV.back()<<"\n";
		}
		
	}//while (myfile.good())
	
	//G4cout << "\nNumber of lines in text file [KERMA]: " << number_of_lines<<"\n";

	KVAL=new double[ndatak];
	for (int i=0;i<ndatak;i++)
		KVAL[i]=pV[i];//photonsToKermainAir conversion factors; KVAL[uGy*mm2/photons]=>YS*kVAL=>KERMA[uGy]/(mAs) at 750 mm

	myfile.close();
	//G4cout << "\nLast 2 pkerma: " << kval_array[ndatak-2]<<"  "<<kval_array[ndatak-1]<<"\n";
	//--------------------	 
}

void XRayBuilder::readAttCoef(string filename, int j)
{
	string filenam = "";// datas+file_sep+dataspec+file_sep+enerfile+defaultext;
	string file_sep=FILE_SEPARATOR;	
	filenam = dataS + file_sep + data_filterS + file_sep + filename
				+ filterKermaExt;

	int number_of_lines = 0;
    string line;
    ifstream myfile(filenam.c_str());//filenam);

	if (!myfile.is_open()){
		  G4cout<<"Error opening file FILTERS!!!!!!!!!!!!!!!! \n";
		  return;
	}

	NOATT++;//Number of attenuators is increased each time attcoef are read!

	int ndatak=0;
	int dataCounter=0;//in line data counter
	double temp;
	vector<double> enV;//same range as readSpectra and readKerma
	vector<double> pV;

	while (myfile.good()){
		number_of_lines++;
		getline(myfile, line);//read line and put it in line variable.
        
		char* str=new char[line.length()+1];
		strcpy(str,line.c_str());//copy line into str.

		if (number_of_lines>0){
			char* pch;
			pch=strtok(str," ,");//space and comma

			while (pch != NULL)
			{
				dataCounter++;
				if(dataCounter==1){
					istringstream buffer(pch);
					buffer>>temp;
					enV.push_back(temp);
				}
				if(dataCounter==2){
					istringstream buffer(pch);
					buffer>>temp;
					pV.push_back(temp);
					ndatak++;
					dataCounter=0;//reset data counter
				}
				pch = strtok (NULL, " ,");//next data and also break the loop if next data = null
			}
			//G4cout <<enV.back()<<"  "<<pV.back()<<"\n";
		}
		
	}//while (myfile.good())
	
	//G4cout << "\nNumber of lines in text file [FILTER]: " << number_of_lines<<" ndatak "<<ndatak<<"\n";
	//----------------------------------
	for(int i=0;i<ndata;i++)//limited to ndata!!!!!!!!!!!
		ATT[i][j-1]=pV[i];

	myfile.close();	
}

void XRayBuilder::buildSpectra(){
    //initial YS==at spectrum Read
	double d = kv / DE;// number of energies!!
	int N = d;//.intValue();
	NEFF = N;//store it
	YF=new double[N];/////////////////////INITIALIZATION
	for (int i = 1; i <= N; i++) {
		YF[i - 1] = 1.0;//initialization
	}
	for (int j = 1; j <= NOATT; j++) {//multiple attenuators
		for (int i = 1; i <= N; i++) {
			if (ATT[i - 1][j - 1] * TMM[j - 1] > 1440) {//set YS
				YS[i - 1] = 0.0;
			}
			else {
				YS[i - 1] = XVAL[i - 1]	* exp(-ATT[i - 1][j - 1] * TMM[j - 1]);//attenuation
			}

			if (j == 1) {
				YF[i - 1] = YS[i - 1];//new photon flux
			} else {
				if (XVAL[i - 1] != 0.0)
					YF[i - 1] = YF[i - 1] * YS[i - 1] / XVAL[i - 1];//old YF * exp(...); Note: exp(...)=YS/XVAL!! OK!
				else
					YF[i - 1] = 0.0;
			}
		}
	}

	TOT_PHOTONS = 0.0;
	for (int i = 1; i <= N; i++) {
		TOT_PHOTONS = TOT_PHOTONS + YF[i - 1];
	}
}

void XRayBuilder::resetYS(){
	for (int i = 1; i <= NEFF; i++) {
			YS[i - 1] = YF[i - 1];
	}
}

void XRayBuilder::computeHVL1(string filename, bool tissB) {

		double t1 = 0.0;
		double t2 = 0.0;
		double t3 = 0.0;
		MEANSPECTRUMENERGY = 0.0;
		KERMAPERMASAT750MM = 0.0;
		int kx = 0;
		double tt = 0.0;
		double tlo = 0.0;
		double thi = 20.0;// MAXIMUM SEAK!!!
		// String filename="AL";
		readAttCoef(filename);// =>ATTVAL

		resetYS();

		while (true) {
			t1 = 0.0;
			t2 = 0.0;
			t3 = 0.0;
			kx = kx + 1;
			for (int l = 1; l <= NEFF; l++) {
				t1 = t1 + YS[l - 1] * KVAL[l - 1];// sum yiKermai=>uGy/mAs at 750 mm
				t2 = t2 + YS[l - 1] * DE * l;// sum yiEi=>keV/(mAs*mm2)
				t3 = t3 + YS[l - 1];// sum yi=>photons/(mAs*mm2) at 750mm
			}
			if (kx == 1) {
				KERMAPERMASAT750MM = t1;//=>uGy/mAs at 750 mm
				MEANSPECTRUMENERGY = t2 / t3;//=>keV ok!
				thi = THI_HVL1;
			} else {
				if (t1 > KERMAPERMASAT750MM / 2.0) {
					tlo = tt;
				} else {
					thi = tt;
				}
				if (t1 < (KERMAPERMASAT750MM / 2.0) * 1.00005
						&& t1 > (KERMAPERMASAT750MM / 2.0) * 0.99995)
					break;
			}
			if (tissB) {
				HVL1_TISS = (tlo + thi) / 2.0;// 1st guess
				tt = HVL1_TISS;
			} else {
				HVL1 = (tlo + thi) / 2.0;// 1st guess
				tt = HVL1;
			}
			
			for (int i = 1; i <= NEFF; i++) {
				if (ATTVAL[i - 1] * tt > 1440) {
					YS[i - 1] = 0.0;
				} else {
					YS[i - 1] = YF[i - 1] * exp(-ATTVAL[i - 1] * tt);
				}
			}
			if (kx > MAXITER) {
				break;
			}
		}

		G4cout <<"\n--------------Spectrum summary------------"<< G4endl;
		G4cout <<"uGy/mAs at 750 mm in air= "<<KERMAPERMASAT750MM<< G4endl;
		G4cout <<"Mean energy in keV = "<<MEANSPECTRUMENERGY<< G4endl;
		G4cout <<"Photons/(mm2*uGy) = "<<Get_photons_per_mm2_per_uGy()<< G4endl;
		if (!tissB)
			G4cout <<"HVL1 in mm"<<filename<<" = "<<HVL1<< G4endl;		
		else
			G4cout <<"HVL1 in mm"<<filename<<" = "<<HVL1_TISS<< G4endl;
	}

void XRayBuilder::readAttCoef(string filename)
{
	string filenam = "";
	string file_sep=FILE_SEPARATOR;	
	filenam = dataS + file_sep + data_filterS + file_sep + filename
				+ filterKermaExt;

	int number_of_lines = 0;
    string line;
    ifstream myfile(filenam.c_str());//filenam);

	if (!myfile.is_open()){
		  G4cout<<"Error opening file FILTERS SINGLE CALL!!!!!!!!!!!!!!!! \n";
		  return;
	}

	int ndatak=0;
	int dataCounter=0;//in line data counter
	double temp;
	vector<double> enV;//same range as readSpectra and readKerma
	vector<double> pV;

	while (myfile.good()){
		number_of_lines++;
		getline(myfile, line);//read line and put it in line variable.
        
		char* str=new char[line.length()+1];
		strcpy(str,line.c_str());//copy line into str.

		if (number_of_lines>0){
			char* pch;
			pch=strtok(str," ,");//space and comma

			while (pch != NULL)
			{
				dataCounter++;
				if(dataCounter==1){
					istringstream buffer(pch);
					buffer>>temp;
					enV.push_back(temp);
				}
				if(dataCounter==2){
					istringstream buffer(pch);
					buffer>>temp;
					pV.push_back(temp);
					ndatak++;
					dataCounter=0;//reset data counter
				}
				pch = strtok (NULL, " ,");//next data and also break the loop if next data = null
			}
			//G4cout <<enV.back()<<"  "<<pV.back()<<"\n";
		}
		
	}//while (myfile.good())
	
	//G4cout << "\nNumber of lines in text file [FILTER]: " << number_of_lines<<" ndatak "<<ndatak<<"\n";
	//----------------------------------
	ATTVAL=new double[ndatak];
	for(int i=0;i<ndatak;i++)
		ATTVAL[i]=pV[i];

	myfile.close();	
}
	