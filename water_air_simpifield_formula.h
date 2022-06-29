#include "base_formulas_original.h"

using namespace std;

//Global water parameters
double flim1w, flim2w, flim3w, flim4w;
double n1w, n2w, n3w, n4w, n5w;
double a1w, a2w, a3w, a4w, a5w;
const int NParsw = 4;
bool water_setup = false;

//Global air parameters
double flim1a, flim2a, flim3a, flim4a;
double n1a, n2a, n3a, n4a, n5a;
double a1a, a2a, a3a, a4a, a5a;
const int NParsa = 3;
bool air_setup = false;

void setup_parameters_water(double TCelsius, double Salinity, double Depth, double pH)
{
  //Basic parameters for the first estimate
  double fW = FrW(TCelsius, Salinity, Depth, pH)*1000; //in Hz
  double fB = FrB(TCelsius, Salinity, Depth, pH)*1000; //in Hz
  double fM = FrM(TCelsius, Salinity, Depth, pH)*1000; //in Hz

  double pars[NParsw]= {TCelsius,Salinity,Depth,pH};
  double cW =   WFuncWater(&fW,&pars[0]); //c=  Func(f0), as then f^2/f0^2       = 1
  double cB = 2*BFuncWater(&fB,&pars[0]); //c=2*Func(f0), as then f^2/(f^2+f0^2) = 1/2
  double cM = 2*MFuncWater(&fM,&pars[0]); //c=2*Func(f0), as then f^2/(f^2+f0^2) = 1/2
  
  double fiBM = sqrt(fB*fM)*pow(cB/cM+fB*fB/fM/fM,0.25); //f_i for transition from B to M
  double fiMW = sqrt(fM*fW)*pow(cM/cW+fM*fM/fW/fW,0.25); //f_i for transition from M to W
  
  double nBM = 4*cM*fiBM*fiBM / (cM*fiBM*fiBM + cB*fM*fM + cM*fB*fB);
  double nMW = 4*cW*fiMW*fiMW / (cW*fiMW*fiMW + cM*fW*fW + cW*fM*fM);
  
  double ciBM = AFuncWater(&fiBM,&pars[0]);
  double ciMW = AFuncWater(&fiMW,&pars[0]);
  
  //Estimated parameters  
  flim1w = fiBM*pow(cB/ciBM*fiBM*fiBM/fB/fB,1.0/(nBM-2));
  flim2w = fiBM*pow(cM/ciBM*fiBM*fiBM/fM/fM,1.0/(nBM-2));
  flim3w = fiMW*pow(cM/ciMW*fiMW*fiMW/fM/fM,1.0/(nMW-2));
  flim4w = fiMW*pow(cW/ciMW*fiMW*fiMW/fW/fW,1.0/(nMW-2));
  
  n1w = 2;
  n2w = nBM;
  n3w = 2;
  n4w = nMW;
  n5w = 2;

  //Estimates based on Excel calculation
  flim1w = 0.8282*flim1w + 25.818;
  flim2w = 0.7481*flim2w + 647.67;
  //flim3w and flim4w are not relevant for the interval [20 Hz, 20 kHz], hence no correction is introduced for them
  
  n1w = 1.9093 + 0.000968093*TCelsius + 0.000928211*Salinity + 4.61102e-7*TCelsius*Salinity -  9.12418e-7*Depth - 4.34584e-10*Salinity*Depth;
  n2w = 0.8716*n2w + 0.2902;
  n3w = 1.51232 - 0.000024935*Depth + 0.00661563*Salinity +  0.00219034*TCelsius - 0.000312489*TCelsius*TCelsius;
  //n4w and n5w are not relevant for the interval [20 Hz, 20 kHz], hence no correction is introduced for them
  
  a1w = 0.0000002651151*(1 - 0.040785676*TCelsius+0.0006424371*TCelsius*TCelsius)*(1-0.008556518*Salinity)*(1-0.000012658*Depth);
  
  //Formulas for the constant prefactors ensuring continuity
  a2w = a1w*pow(flim1w,n1w-n2w);
  a3w = a2w*pow(flim2w,n2w-n3w);
  a4w = a3w*pow(flim3w,n3w-n4w);
  a5w = a4w*pow(flim4w,n4w-n5w);
  
  water_setup = true;
}

void setup_parameters_air(double Pressure, double TCelsius, double Humidity)
{
  double PPascal = Pressure * Patm;
  double TKelvin = TCelsius + T01;
  double hvalue = h(Humidity,TKelvin,PPascal);
  
  //Basic parameters for the first estimate
  double fN = FrN(TKelvin,hvalue,PPascal);
  double fO = FrO(hvalue,PPascal);
  double fC = 1; //fC does not appear in the formulas, but can be understood as cC -> cC/fC^2, for symmetry reasons
  
  double pars[NParsa]= {Pressure,TCelsius,Humidity};
  
  double cN = 2*NFuncAir(&fN,&pars[0]); //cN=2*NFunc(fN), as then f^2/(f^2+fN^2) = 1/2
  double cO = 2*OFuncAir(&fO,&pars[0]); //cO=2*OFunc(fO), as then f^2/(f^2+fO^2) = 1/2
  double cC =   CFuncAir(&fC,&pars[0]); //cO=CFunc(fC), as then f^2=1
  
  double fiNO = sqrt(fN*fO)*pow(cN/cO+fN*fN/fO/fO,0.25); //f_i for transition from N to O
  double fiOC = sqrt(fO*fC)*pow(cO/cC+fO*fO/fC/fC,0.25); //f_i for transition from O to classical
  
  double nNO = 4*cO*fiNO*fiNO / (cO*fiNO*fiNO + cN*fO*fO + cO*fN*fN);
  double nOC = 4*cC*fiOC*fiOC / (cC*fiOC*fiOC + cO*fC*fC + cC*fO*fO);
  
  double ciNO = AFuncAir(&fiNO,&pars[0]);
  double ciOC = AFuncAir(&fiOC,&pars[0]);
  
  //Estimated parameters  
  flim1a = fiNO*pow(cN/ciNO*fiNO*fiNO/fN/fN,1.0/(nNO-2));
  flim2a = fiNO*pow(cO/ciNO*fiNO*fiNO/fO/fO,1.0/(nNO-2));
  flim3a = fiOC*pow(cO/ciOC*fiOC*fiOC/fO/fO,1.0/(nOC-2));
  flim4a = fiOC*pow(cC/ciOC*fiOC*fiOC/fC/fC,1.0/(nOC-2));
  
  n1a = 2;
  n2a = nNO;
  n3a = 2;
  n4a = nOC;
  n5a = 2;
  
  //Corrections based on Excel calculation
  flim1a = 0.92306*flim1a - 8.2963E-05*flim1a*flim1a; 
  flim2a = 0.90727*flim2a - 1.5649E-05*flim2a*flim2a;
  if(flim3a<0 || flim3a>20000) flim3a = 22222; //This means that if original flim3 is out of range, then do not use it
  else flim3a = 0.89696*flim3a + 5.3452E-07*flim3a*flim3a;
  //flim4a is not relevant for the interval [20 Hz, 20 kHz], so no correction is introduced for it
  
  n1a = 1.763454532 + 0.003249811*TCelsius + 0.001091377*Humidity;
  n2a = 1.6968*n2a - 0.537*n2a*n2a;
  n3a = 1.671471966 + 0.000832291*TCelsius + 0.001112712*Humidity;
  n4a = 1.7413*n4a;
  //n5a is not relevant for the interval [20 Hz, 20 kHz], so no correction is introduced for it
  
  double a1temp = 1.96005E-05*exp(-0.0159*TCelsius)*pow(Humidity,-0.4);
  a1a = -0.83774*a1temp + 7.52E+05*a1temp*a1temp + 2.35E+10*a1temp*a1temp*a1temp;
  
  //Formulas for the constant prefactors ensuring continuity
  a2a = a1a*pow(flim1a,n1a-n2a);
  a3a = a2a*pow(flim2a,n2a-n3a);
  a4a = a3a*pow(flim3a,n3a-n4a);
  a5a = a4a*pow(flim4a,n4a-n5a);
  
  air_setup = true;
}

double water_absorption(const double *x, const double *pars)
{
  double f = x[0]; //in Hz
  double TCelsius = pars[0];
  double Salinity = pars[1];
  double Depth    = pars[2];
  double pH       = pars[3];
  
  if(!water_setup) setup_parameters_water(TCelsius, Salinity, Depth, pH);
  
  if(f<flim1w)      return a1w*pow(f,n1w);
  else if(f<flim2w) return a2w*pow(f,n2w);
  else if(f<flim3w) return a3w*pow(f,n3w);
  else if(f<flim4w) return a4w*pow(f,n4w);
  else              return a5w*pow(f,n5w);
}

double air_absorption(const double *x, const double *pars)
{
  double f = x[0]; //in Hz
  double Pressure = pars[0];
  double TCelsius = pars[1];
  double Humidity = pars[2];
  
  if(!air_setup) setup_parameters_air(Pressure, TCelsius, Humidity);
  
  if(f<flim1a)      return a1a*pow(f,n1a);
  else if(f<flim2a) return a2a*pow(f,n2a);
  else if(f<flim3a) return a3a*pow(f,n3a);
  else if(f<flim4a) return a4a*pow(f,n4a);
  else              return a5a*pow(f,n5a);
}