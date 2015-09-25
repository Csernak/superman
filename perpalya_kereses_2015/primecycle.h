//============================================================================
// Name        : primecycle.cpp
// Author      : Gergely Gyebroszki
// Version     :
// Copyright   : 2014 (C) gyebro.com
// Description : Prime Cycle Generator
//============================================================================

#ifndef PRIMECYCLE_H
#define PRIMECYCLE_H

#include "cycle.h"
#include "math_extensions.h"
#include "armadillo"

using namespace std;
using namespace arma;

namespace MicroChaosTools {
  
  /********************- SYSTEMINFO ************************/
class SystemInfo //A rendszer parameterei, matrixok kulonbozo hatvanyai
  {
  public:
    vector<mat> U; //Az U^i, az U matrix és hatvanyai
    vector<vec> UB; //U^i*b vektorok
    vector<mat> IUinv; //(Imat-U^n) inverzei
    double param[4];
    double P;
    double D;
    double alpha;
    double beta;
    int mmax;
    int elmmax;
    
    //Konstruktor; a maximalis periodus hossz es a parametereket tartalmazo tomb a bemenete.
    //Csinalhatna a stabilitasszamitast is!!
    SystemInfo(size_t length,double* params,unsigned char alphabet){
    mat Um; //Az inicializalashoz kellenek ezek.
    vec B;
    //double 
    P = params[0];
    //double 
    D = params[1];
    //double 
    alpha = params[2];
    //double 
    beta = params[3];
   /*
    param[0] = P;
    param[1] = D;
    param[2] = alpha;
    param[3] = beta;*/
    double gamma = sqrt(alpha*alpha+beta*beta);
    double e = exp(beta);
    double s = sinh(gamma);
    double c = cosh(gamma);
    elmmax = floor(P/(P-alpha*alpha)); 
    mmax = (int)((alphabet-1)/2);
    
    if(mmax < elmmax)
    {cout << "Alphabet is shorter than theoretical!\n";}
    
    //Lenullazzuk a matrixokat:
    U.resize(0); //Ennek 0., 1., 2., i., stb eleme: I, U, U^2, U^i, stb
    UB.resize(0); //Ennek elemei: b, U b, U^(i) b
    IUinv.resize(0); //Ennek elemei: (I-U^(i+1))^(-1)
    
/*MATRIXOK GENERALASA*/
//U inicializalasa:
    Um << 1 << 0 <<endr << 0 << 1 << endr; //Az U elso (nulla indexu) eleme az egysegmatrix lesz
    U.push_back(Um); //Betesszuk U-ba
    
      //U masodik (1 indexu) eleme maga az U matrix:	
    Um << gamma* c+beta* s  << s << endr 
    << alpha*alpha* s << gamma* c-beta* s << endr;
    Um = Um/(e*gamma);
    U.push_back(Um);	
    
//UB inicializalasa:
	B << (e-c)/(e*alpha*alpha)-(beta*s)/(e*gamma*alpha*alpha)<< endr << -s/(e*gamma) << endr;
	UB.push_back(B); //Elso (nulla indexu) elem: b
	UB.push_back(U[1]*B); //Masodik (1 indexu) elem: U b.

//IUinv inicializalasa:	
	IUinv.push_back(inv(U[0]-U[1]));//Elso (0 indexu) elem: (I-U)^(-1) 
	
     for(size_t i = 2; i <= length; i++){
	U.push_back(U[1]*U[i-1]); //U-nak a 2 indexu (3. elemetol) a length indexu (length+1-edik elemeig)  
	UB.push_back(U[1]*UB[i-1]); //UB indexei hasonloan
	IUinv.push_back(inv(U[0]-U[i])); //IUinv 1 indexu (2.) elemetol (I-U^2) a length-1 indexu (length-edik) elemig, (I-U^length) 
      }//for
      //cout << U << '\n' << UB << '\n' << IUinv << '\n';
    }//konstruktor vege
    
};//SystemInfo vege
 
 /************************* PERIODIC ORBIT ***********************************/
  
  //Ez tartalmazza a periodikus palyakat
template <class T> class PeriodicOrbit
{
public:
  bool isvalid;
  vector<vec> palya; //koordinata vektorok vektora
   vec ce; 
 /******************* KONSTRUKTOR *************************/
 PeriodicOrbit(Cycle<T> c, SystemInfo SI){ 
    //KONSTRUKTOR; kap egy c prim ciklust es a rendszer jellemzoit.
    //Ez alapjan legyartja a lehetseges palya pontokat es kilep, ha nem jo valamelyik.
   ce.resize(c.length());
    for (size_t i=0; i<c.length(); i++) {
     ce(i) = int(c(i))-SI.mmax; 
    }//for
    isvalid = true;
    palya.resize(0);   
    isvalid = test(c, SI); //A test() fv. gyartja le a pontokat, az hivja meg az ellenort().
    };//Konstruktor vege

 /********************** TEST ******************************/
 //template <class T> 
 bool test(Cycle<T> c, SystemInfo SI) {
    vec Vektor(2);
    bool isvalid = true;
    Vektor << 0 << endr << 0 << endr; //0,0 oszlopvektor
    vector<mat> y; //Ebben lesznek az adott per. palya pontjai
    y.resize(0);
    vec cm(c.length()); //A per. palya indexeket egesz tombbe kell tenni, negativ is lehet!
     //cm.resize(0);   
    //c elemei 0-tol (alphabet-1)-ig lehetnek indexelve.
    //Nekunk az kell, hogy -mmax-tol mmax-ig legyenek, ezert 
    // c minden elemebol ki kell vonni mmax-ot:   
    for (size_t i=0; i<c.length(); i++) {
     cm(i) = int(c(i))-SI.mmax; 
    }//for
    vec y0(2); 
   
     // MOST gyartjuk le a periodikus palyak pontjait:
        for(int k = 1; k <= c.length(); k++)
	{
            Vektor += (cm(c.length()-k))*SI.UB[k-1];       
	} //for k, ciklus a j-edik prim ciklus elemeire
  
  //Ezek a lehetseges periodikus palyak kezdopontjai, ez eppen a j-edik:
    y0 = SI.IUinv[c.length()-1]*Vektor;

   y.push_back(y0);

  if(int(trunc(SI.P*y0(0)+SI.D*y0(1))) != cm(0)){
  
    isvalid = false; return isvalid;}
  //  %fixitt = Perpalyak(j,1)/alpha^2
   // y 1 = (I ? U p ) ?1 m 1 U p?1 + m 2 U p?2 + . . . + m p U 0 b
    //%Meg ellenorizni kell, hogy ott vannak-e, ahol kell:    
  
  if(isvalid){
        for(int p=1; p <= c.length()-1; p++)
	{
	  y0 = SI.U[1]*y[p-1]+cm(p-1)*SI.UB[0];
	  y.push_back(y0);
	
	    if(int(trunc(SI.P*y0(0)+SI.D*y0(1))) != cm(p)){
	      isvalid = false; return isvalid;}        
    }//for
    
    }//Ha az elso pont oke.
    if(isvalid) //%Ha van rendben levo ciklus, akkor fajlba irjuk a koordinatakat is (es majd azt is, hogy melyik savban vannak, lasd print):
     {
      //PeriodicOrbit<T>::
      for(int i = 0; i < c.length(); i++)
      {
	palya.push_back(y[i]);
      }
    
      return isvalid;
    }
    
}//test vege   
    
//Palyak kozeppontjainak, egyeb adatainak kiiratasa
void analyse(ofstream& out, Cycle<T> c,SystemInfo SI)
{  //Az eredmenyek analizise
  double szummx;
  double szummy;
  double maxix = 0;
  double minix = 0;
  double maxiy = 0;
  double miniy = 0;
  vec y0(2);
//Palya kozeppontok 
  for (size_t i=0; i < palya.size(); i++) {
    szummx+=PeriodicOrbit<T>::palya[i](0);
    szummy+=PeriodicOrbit<T>::palya[i](1);
    } 
    szummx = szummx/double(palya.size());
    szummy = szummy/double(palya.size());
    out << szummx << '\t' << szummy << "\t" << palya.size();// << '\n';
  
//Kumulalt hiba ellenorzes
    y0 = PeriodicOrbit<T>::palya[0];
    double error = 0;
      //A periodikus palya atmeroje
    maxix = y0(0);
    minix = y0(0);
    maxiy = y0(1);
    miniy = y0(1);
    for(int i = 1; i < c.length(); i++)
    {
    y0 = SI.U[1]*y0+trunc(SI.P*y0(0)+SI.D*y0(1))*SI.UB[0];
    error += sqrt((y0(0)-PeriodicOrbit<T>::palya[i](0))*(y0(0)-PeriodicOrbit<T>::palya[i](0))+(y0(1)-PeriodicOrbit<T>::palya[i](1))*(y0(1)-PeriodicOrbit<T>::palya[i](1)));
      maxix = max(maxix,PeriodicOrbit<T>::palya[i](0));
    minix = min(minix,PeriodicOrbit<T>::palya[i](0));
    maxiy = max(maxiy,PeriodicOrbit<T>::palya[i](1));
    miniy = min(miniy,PeriodicOrbit<T>::palya[i](1));
      
    }
    out << '\t' << error;
    out << '\t' << maxix;
    out << '\t' << minix;
    out << '\t' << maxiy;
    out << '\t' << miniy;
    if(error > 1e-10) {cout << "Error: " << error << '\n';}

//Mekkora a legnagyobb ugras egy palya menten
    int jump_index = 0;
    int maxjump_index = abs(c(0)-c(c.length()-1));  
    double jump = 0;
    double maxjump  = norm(PeriodicOrbit<T>::palya[0]-PeriodicOrbit<T>::palya[c.length()-1]);
    for(size_t i = 0; i < c.length()-1; i++)
    {
    jump_index = abs(c(i)-c(i+1));
    if(jump_index > maxjump_index){maxjump_index = jump_index;}
    jump = norm(PeriodicOrbit<T>::palya[i]-PeriodicOrbit<T>::palya[i+1]);;
    if(jump > maxjump){maxjump = jump;}          
    }
    out << '\t' << maxjump_index << '\t' << maxjump;  
    
//Hany savot erint a palya
/*vec cm(c.length());
cm.fill(0);
for(size_t i = 0; i < c.length(); i++)
{
  cm(i) = c(i);
}
*/
vec cm2 =  unique(ce);
//Ennyi savot erint:
out << '\t' << cm2.n_rows;

//Minimalis es maximalis indexu sav:
    int max_index = -2*SI.mmax;
    int min_index = 2*SI.mmax;  
    
    for(size_t i = 0; i < c.length(); i++)
    {
    if(ce(i) > max_index) {max_index = ce(i);};
    if(ce(i) < min_index) {min_index = ce(i);};     
    }
    out << '\t' << min_index << '\t' << max_index;   
  
    out  << '\n';
    fflush(stdout);
    return;
  }//analyse vege
  //*****************************************************************
  
  //Palyapontok és a sáv kiirasa  
  void print(ofstream& out) 
  {
    for (size_t i=0; i < PeriodicOrbit<T>::palya.size(); i++) {
      out << palya[i](0) << '\t' << /*PeriodicOrbit<T>::*/palya[i](1) << '\t' << ce(i) << '\n';
    }
    out << '\n';
  }//print vege


  
};//PeriodicOrbit class vege
  

/**
 * Number of Prime Cycles with length n in case of full N-symbol dynamics
 * See ChaosBook Ch 15.7
 */
template < class T >
unsigned long long int numberOfPrimeCycles(T alphabet, size_t length) {
	vector<size_t> divs = divisors(length);
	unsigned long long int N = (unsigned long long int)alphabet;
	unsigned long long int numPC = uintpow(N,length);
	for(size_t i=0; i<divs.size(); i++) {
		numPC -= divs[i]*numberOfPrimeCycles(alphabet,divs[i]);
	}
	numPC /= length;
	return numPC;
};//numberOfPrimeCycles vege

/**
 * Number of Prime Cycles with length up to n, in case of full N-symbol dyn.
 */
template < class T >
unsigned long long int numberOfAllPrimeCycles(T alphabet, size_t length) {
	unsigned long long int allPC = 0;
	for (size_t i=1; i<=length; i++) {
		allPC += numberOfPrimeCycles(alphabet,i);
	}
	return allPC;
};//numberOfAllPrimeCycles vege


enum ordering {
	ascending,
	descending
}; //ordering vege


//EZ SZÁMÍT IGAZÁN
template < class T >
unsigned long long int generatePrimeCycles(T alphabet, size_t length, ordering o, ofstream& output,bool anal, ofstream& anafile,SystemInfo SI) {
	Cycle<T> c(alphabet,length);		
	c.fill(0);
	unsigned long long int generated = 0;
	if(o == ascending) {
		if(c.isLowestCyclicPermutation()) 
		{
		  PeriodicOrbit<T> po(c,SI);
		  if(po.isvalid)
		  {
			po.print(output);			
			if(anal){po.analyse(anafile,c,SI);};
			generated++;
		  }
		
		}
		while(c.increase()){
			if(c.isLowestCyclicPermutation() ) 
			{
			      PeriodicOrbit<T> po2(c,SI);
			      if(po2.isvalid)
			      {
				      po2.print(output);
				      if(anal){po2.analyse(anafile,c,SI);};
				      generated++;
			      }
			}
		}
	}
	else if (o == descending) {
		if(c.isHighestCyclicPermutation()) {
		  PeriodicOrbit<T> po(c,SI);
		  if(po.isvalid)
		  {
			po.print(output);
			if(anal){po.analyse(anafile,c,SI);};
			generated++;
		  }
		}
		while(c.increase()){
			if(c.isHighestCyclicPermutation() ) {
				 PeriodicOrbit<T> po2(c,SI);
			      if(po2.isvalid)
			      {
				po2.print(output);
				if(anal){po2.analyse(anafile,c,SI);};
				generated++;
			      }
			
			}
		}
	}
	return generated;
};//generatePrimeCycles vege


} /*end of namespace MicroChaosTools */

#endif
