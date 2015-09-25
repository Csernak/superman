//============================================================================
// Name        : PrimeCycles.cpp
// Author      : Gergely Gyebroszki + Gabor Csernak
// Version     :
// Copyright   : 2014 (C) gyebro.com
// Description : Symbolic Prime Cycle and periodic Orbit Generator for the PD controlled oscillator
//====================, 2*========================================================
// Lenyegeben ugyanazt adja, mint a regi Maple alapu szamitas, ami brute force modon vegigment minden index permutacion. 
// Ott az osszes per. pontot kapjuk meg, azaz pl. a 4 periodusu pontoknal megtalal n1 db. 1 periodusut, 2*n2 db ket periodusut es 4*n4 db. 
// negy periodusut, mig ez a program csak a prim ciklusokhoz tartozo per. palyakat (n1,n2,n3,n4,stb) adja ki. Ezt figyelembe veve a ket 
// program a megfelelo szamu pontot (palyat) talalja meg.
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
using namespace std;

#include "primecycle.h"
#include "cmd_line_arg_utils.h"
using namespace MicroChaosTools;
#define PLOT 0
#define RUNS_ON_CONDOR

int main(int argc, char *argv[]) {
/**************** ALAPVETO INFOK, PARAMETEREK BEOLVASASA ***************************/
	cout << "Prime Cycle and Periodic Orbit Generator" << endl;
	for (int i=0; i<argc; i++) {
		cout << argv[i] << " ";
	};
	cout << endl;
	CommandLineArgumentParser clap(argc, argv);
	string symbolString; //abc
	string lengthString; //periodus hossz
	string orderString; //Milyen sorrend
	string PparamString;
	string DparamString;
	string alphaparamString;
	string betaparamString;
	string analString; //Legyen-e analizis
	bool ok = true;
	//Nem kell feltetlenul megadni az abc-t, a tobbi alapjan javasol.
	clap.GetString("-n",symbolString); 
	//Meg kell adni, hogy milyen periodushosszig menjen el:
	ok &= clap.GetString("-l",lengthString);
	//Be kell adni a lekepezes parametereit is:
	ok &= clap.GetString("-P",PparamString);
	ok &= clap.GetString("-D",DparamString);
	ok &= clap.GetString("-a",alphaparamString);
	ok &= clap.GetString("-b",betaparamString);
	clap.GetString("-o",orderString);
	clap.GetString("-c",analString);

	if (!ok) {
		cout << "usage:\n"
				"PrimeCycles -n (number of symbols) -l [cycle length] -o (order) -c (analysis)\n"
				"-P [parameter P]\n"
				"-D [parameter D]\n"
				"-a [parameter alpha]\n"
				"-b [parameter beta]\n"
				" (number of symbols): optional, positive integer specifying the length of alphabet\n"
				" [cycle length]: positive integer specifying the length of cycles\n"
				" (order): optional, 0 for ascending, 1 for descending order, default = 0\n"
				" (analysis): optional, 0 for analysis, 1 for no analysis, default = 0\n";				
		return 1;
	}
//************** Bemeneti parameterek valtozoba toltese:	
	int order = atoi(orderString.c_str());
	enum ordering o;
	if (order == 1) {
		o = descending;
	}
	else {
		o = ascending;
	}
    
	int an = atoi(analString.c_str());
	bool anal;
	if (an == 1) {
		anal = false;
	}
	else {
		anal = true;
	}
	
   	
//***************Lekepezes parameterei:		
	size_t length = atoi(lengthString.c_str());
	double P = (double)atof(PparamString.c_str());
#ifdef RUNS_ON_CONDOR
	// P 0 es 48 kozott valtozik egyesevel
	P = 0.7 + P*0.1;
#endif
	double D = (double)atof(DparamString.c_str());
	double alpha = (double)atof(alphaparamString.c_str());
	double beta = (double)atof(betaparamString.c_str());
	int elmmax = floor(P/(P-alpha*alpha)); //Elso virtualis fixpont indexe (utolso lehetseges attraktor indexe), elméleti érték
	int mmax; 
//***********Uj segedvaltozok:
	double gamma = sqrt(alpha*alpha+beta*beta);
	double e = exp(beta);
	double s = sinh(gamma);
	double c = cosh(gamma);
	double param[4] = {P,D,alpha,beta};
//************Stabilitas szamitashoz:    
	bool stable = true; 
	double Pmin = alpha*alpha;
	double Pmax;
	double Dvalt = ((-3*e + 2*c*c*e + e*e*e + c*(-1 + e*e))*gamma +  beta*(1 - 2*c*e - 3*e*e)*s)/((-1 + 2*c*e - e*e)*s);
	double Pmax1 = (alpha*alpha*(gamma-e*e*gamma-D*e*s))/(gamma-c*e*gamma+beta*e*s);
	double Pmax2 = (alpha*alpha*((1+2*c*e+e*e)*gamma-2*D*e*s))/(gamma-e*e*gamma+2*beta*e*s);
	if(D < Dvalt) {
		Pmax = Pmax1;
	} else {
		Pmax = Pmax2;
	}

	double Dmax = (alpha*alpha*(1 + 2*c*e + e*e)*gamma + (-1 + e*e)*gamma*P - 2*beta*e*P*s)/(2*alpha*alpha*e*s);
	double Dmin = (alpha*alpha*(gamma-e*e*gamma)+P*((-1+c*e)*gamma-beta*e*s))/(alpha*alpha*e*s);

	if(P < Pmin) {
		cout << "Unstable, too small P\n";
		cout << "Pmin: " << Pmin << ", Pmax: " << Pmax << ", Dmin: "<< Dmin << ", Dmax: " << Dmax << '\n';
		stable = false;
	}
	if(P > Pmax) {
		cout << "Unstable, too large P\n";
		cout << "Pmin: " << Pmin << ", Pmax: " << Pmax << ", Dmin: "<< Dmin << ", Dmax: " << Dmax << '\n';
		stable = false;
	}    
	if(D > Dmax) {
		cout << "Unstable, too large D\n";      
		cout << "Pmin: " << Pmin << ", Pmax: " << Pmax << ", Dmin: "<< Dmin << ", Dmax: " << Dmax << '\n';
		stable = false;
	} 
	if(D < Dmin) {
		cout << "Unstable, too small D\n";
		cout << "Pmin: " << Pmin << ", Pmax: " << Pmax << ", Dmin: "<< Dmin << ", Dmax: " << Dmax << '\n';
		stable = false;
	} 
	if(!stable) {exit(-1);}	
	
/******************* Stabilitas szamitas vege ********************************/
//1 < P < 3: alphabet = (unsigned char)(2*elmmax+1); length = 7
//3 < P < 4: alphabet = 15; length = 7
//4 < P < 5: alphabet = 19; length = 7
//5 < P < 5.2: alphabet = 25; length = 8
//5.2 < P < 5.5: alphabet = 33; length = 8

//Ha megadjuk az abc-t a -n kapcsolóval, akkor azt veszi figyelembe, különben a 2 mmax+1-et. 
//Ez csak a parameterek beolvasasa utan johet, mert kell az mmax.
	int abc = atoi(symbolString.c_str());
	unsigned char alphabet; //char volt
	if (abc > 1)  {
		alphabet = (unsigned char)abc;
		mmax = (int)((abc-1)/2);
	} else {
		alphabet = (unsigned char)(2*elmmax+1);
		mmax = elmmax;
	}
	
//*** Mappanev keszitese
	stringstream folder;
	folder << "run_alpha_" << alpha << "_beta_" << beta << "_P_" << P << "_D_" << D;
	string folderstr = folder.str().c_str();
	stringstream foldercmd;
	foldercmd << "mkdir " << folderstr;
	system(foldercmd.str().c_str());
	
//******************** INICIALIZALAS:
	SystemInfo SI(length,param,alphabet); //A matrixok eloallitasa
	
	
//********************** Kimeneti uzenetek es a prim ciklusok szamitasa:
	cout << "Generating Periodic Orbits for\n";
	cout << " N = " << (int)alphabet << " symbols"; 
	if(abc > 1) { 
		cout << " (given)\n";
		cout << "Suggested alphabet for these parameters: N = " << 2*elmmax +1 << '\n';
	} else{cout << " (size of absorbing sphere)\n";}
	cout << " up to l = " << length << " length.\n";
	cout << "Number of prime cycles = " << numberOfAllPrimeCycles(alphabet,length) << endl << endl;
	
/**************** FAJLOK **************************/
	// Nem a prim ciklusokat, hanem az igaziakat mentjuk el. 
	ofstream outfile[length]; // Minden ciklushosszhoz mas fajl azonositot keszit.
	stringstream ss[length]; // Minden ciklushosszhoz elmenti a palya pontok koordinatait; egy sor kimarad ket ciklus kozott.
	// A fajl neve cyclesn.out, ahol n a ciklushossz.
	for(int i = 0; i < length;i++) {
		ss[i] << folderstr << "/cycles" << i+1 << ".out";
		outfile[i].open(ss[i].str().c_str());
		outfile[i] << "#P = " << P << ", D = " << D << ", alpha = " << alpha << ", beta = " << beta << '\n';
	}
	// Ebbe az adott hosszhoz tartozo megtalalt palyak szamat irjuk ki:
	stringstream ss2; 
	ofstream outfile2;
	ss2 << folderstr << "/cyclenum" << (int)alphabet << ".out";
	outfile2.open(ss2.str().c_str());
	
	// Ebbe kerul az analizis eredmenye
	stringstream ananev; 
	ofstream anafile;
	ananev << folderstr << "/analysis.dat";
	anafile.open(ananev.str().c_str());
	if(anal) {
		anafile << "#P = " << P << ", D = " << D << ", alpha = " << alpha << ", beta = " << beta << ", mmax = " << mmax << '\n';
		anafile << "#Cycle center positions: x,y, cycle length: l\n";	 
		anafile << "#x\ty\tl\terror\tmaxix\tminix\tmaxiy\tminiy\tmaxjump_index\tmaxjump\tnumberofbandsvisited\tmin_index\tmax_index\n";
	} else {
		anafile << "No analysis was done.\n";
	}
	
//********************* ITT JON A SZAMITAS *********************	
	unsigned long long int primecount;
	unsigned long long int pccount;		
	
	for (size_t i=1; i<=length; i++) {
		primecount = numberOfPrimeCycles(alphabet,i);
		cout << "Generating cycles of length = " << i << ".\nNumber of primes: " << primecount;
		clock_t t = clock();
		pccount = generatePrimeCycles(alphabet,i,o,outfile[i-1],anal,anafile,SI);
		t = clock() - t;
		cout << "\nNumber of cycles: " << pccount<< '\n';
		cout << (double)(pccount)/(double)(primecount)*100 << " % of prime cycles exists.\n";
		cout << "Done in " << ((float)t)/CLOCKS_PER_SEC << " sec.\n\n";		
		outfile2 << i << '\t' << pccount << '\n';		
		outfile[i-1].close();		
	} // for ciklus a palyahosszokra	
	outfile2.close();
	// if(anal)
	{anafile.close();}
	fflush(stdout); //enelkul nem irta ki az eredmenyeket

	
	//******************************************* ABRAZOLAS
	bool abra = true; //legyen-e abra
	bool perpont = true; //prim periodikus palyak pontjai
	bool szimul = true; //szimulacio
	bool kapcsvonal = true; //kapcsolovonalak
	bool percenter = true; //periodikus palyak kozeppontjai
	bool firstitem = true; //elsokent rajzolt elem-e

	if(abra && (perpont || szimul || kapcsvonal || percenter)) {
		//parameter fajl letrehozasa
		ofstream gnufile; //a gnuplothoz
		stringstream gnu; //Minden ciklushosszhoz elmenti a palya pontok koordinatait; egy sor kimarad ket ciklus kozott.
		gnu << folderstr << "/gnupar";
		gnufile.open(gnu.str().c_str());
		gnufile << "#P = " << P << ", D = " << D << ", alpha = " << alpha << ", beta = " << beta << '\n';
		//parameter fajl letrehozasa az analizishez
		ofstream gnufile2; //a gnuplothoz
		stringstream gnu2; //Minden ciklushosszhoz elmenti a palya pontok koordinatait; egy sor kimarad ket ciklus kozott.
		gnu2 << folderstr << "/gnupar2";
		gnufile2.open(gnu2.str().c_str());
		gnufile2 << "#P = " << P << ", D = " << D << ", alpha = " << alpha << ", beta = " << beta << '\n';
		fflush(stdout);
		//A szimulacio legyen az elso
		if(szimul) {
			ofstream szimfile; //a gnuplothoz
			stringstream szim; //Minden ciklushosszhoz elmenti a palya pontok koordinatait; egy sor kimarad ket ciklus kozott.
			szim << folderstr << "/szimul.dat";
			szimfile.open(szim.str().c_str());
			szimfile << "#P = " << P << ", D = " << D << ", alpha = " << alpha << ", beta = " << beta << '\n';
			int SIMLENGTH = 10000/elmmax;
			vec y(2);
			double maxx = 0;
			double maxy = 0;
			
			for(int i = -elmmax; i <= elmmax; i++)
			{
				y(0) = double(i/P*1.001); y(1) = 0;
				for(int j = 0; j < SIMLENGTH; j++) {
					y = SI.U[1]*y+int(trunc(P*y(0)+D*y(1)))*SI.UB[0];
					if(j > SIMLENGTH/10) {
						szimfile << y(0) << '\t' << y(1) << '\n';
						if(y(0) > maxx){maxx = y(0);}
						if(y(1) > maxy){maxy = y(1);}
					}
				}//szimulacio
				szimfile << '\n';
			}//for a kapcs. vonalakra
			szimfile.close();
			if(firstitem) {
				gnufile << "plot ["<< -maxx << ':' << maxx <<"][" << -maxy << ':'<< maxy << "] 'szimul.dat'";
				firstitem = false;
			} else {
			  gnufile << ", 'szimul.dat'";	  
			}
		}//szimul
		if(perpont){
			if(firstitem){gnufile << "plot "; firstitem = false;}
			else{gnufile << ", ";} 
			for(size_t i = 1; i <= length; i++) {
				gnufile << "'cycles" << i << ".out'";
				if(i < length) {
					gnufile << ", "; 
				}
			}//for
		}//if perpont
		if(kapcsvonal){
			int megmennyit = 3;
			if(firstitem) {
				gnufile << "plot ";
				firstitem = false;
			} else {
				gnufile << ", ";
			}
			for(int i = -mmax-megmennyit; i <= mmax+megmennyit; i++) {
				if(i != 0) {
					gnufile << i/D << '-' << P/D << "*x title ''";	  
				}
				if(i < mmax+megmennyit && i !=0) {
					gnufile << ", "; 
				}
			}//for
		}//kapcsvonal

		//******* GNUPLOT HIVASA ************
		gnufile << ";\n"; //Abrazolo sor vege
		gnufile << "set term post\n set output 'perpontok.ps'\n replot \n";
		gnufile << "pause -1 'Press ENTER to return...';\n";
		gnufile.close();
		char run[36] = ("gnuplot gnupar \0");
		if(PLOT){ cout << system(run); }
		if(anal){
			gnufile2 << "plot 'analysis.dat'";
			gnufile2 << ";\n"; //Abrazolo sor vege
			gnufile2 << "set term post\n set output 'perpontok.ps'\n replot \n";
			gnufile2 << "pause -1 'Press ENTER to return...';\n";
			char run2[36] = ("gnuplot gnupar2\0");
			if(PLOT){cout << system(run2); }
		}
		gnufile2.close();
	}// abrazolas
	
	return 0;
} // main VEGE
