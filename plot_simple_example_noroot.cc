#include <iostream>
#include <fstream>
#include <cmath>
#include "water_air_simpifield_formula.h"
#include "other_materials_ai.h"

using namespace std;

const int Nf = 100;
const double frmin = 20.;      //min frequency in Hz
const double frmax = 20000.;  //max frequency in Hz
double df = log(frmax/frmin)/Nf;        //frequency step inverse exponent

int main()
{
  //Parameter array for double*,double* functions
  double *airpars = new double[3];
  airpars[0] = 20; //TCelsiusAir;
  airpars[1] = 30; //Humidity
  airpars[2] = 1;  //Pressure
  //Printing results
  cout << "Printing air absorption, original and simplified formula" << endl;
  for(int ifr=0;ifr<Nf;ifr++)
  {
    double fr = frmin*exp(ifr*df); // running logarithmically from 20 Hz to 20000 kHz
    cout << fr << "\t" << AFuncAir(&fr,airpars) << "\t" << air_absorption(&fr,airpars) << endl;
  }
  cout << endl;
  delete airpars;
  
  //Parameter array for double*,double* functions
  double *waterpars = new double[4];
  waterpars[0] = 8;   //TCelsiusWater;
  waterpars[1] = 35;  //Salinity
  waterpars[2] = 500; //Depth
  waterpars[3] = 8;   //pH
  //Printing results
  cout << "Printing water absorption, original and simplified formula" << endl;
  for(int ifr=0;ifr<Nf;ifr++)
  {
    double fr = frmin*exp(ifr*df); // running logarithmically from 20 Hz to 20000 kHz
    cout << fr << "\t" << AFuncWater(&fr,waterpars) << "\t" << water_absorption(&fr,waterpars) << endl;
  }
  cout << endl;
  delete waterpars;
  
  //Chosing a material
  double nmat = 5; //"Smooth brickwork with flush pointing"
  //Printing results
  cout << "Printing " << materialnames[(int)floor(nmat)] << " absorption data, original data" << endl;
  for(int ifr=0;ifr<NAI;ifr++)
  {
    double fr = frequenciesOtherMat[ifr];
    cout << fr << "\t" << aivalues[(int)floor(nmat)][ifr] << endl;
  }
  cout << "Printing " << materialnames[(int)floor(nmat)] << " absorption data, simplified" << endl;
  for(int ifr=0;ifr<Nf;ifr++)
  {
    double fr = frmin*exp(ifr*df); // running logarithmically from 20 Hz to 20000 kHz
    cout << fr << "\t" << AlphaOtherMaterial(&fr,&nmat) << endl;
  }
  cout << endl;

  return 0;
}