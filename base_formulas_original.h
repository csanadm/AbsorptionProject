const double T0 = 293.15;     //temperature scale for the calibrated functions
const double Patm = 101325.0; //atmospheric pressure
const double T01 = 273.16;    //triple point temperature

const double KelvinConversion = 273.15;

//classical absorption
double classical(double TKelvin, double f)
{
  return 868.6*1.84e-11*pow(TKelvin/T0,1./2.)*f*f;
}

//Nitrogen relaxation frequency
double FrN(double TKelvin, double h, double P)
{ 
  return P/Patm*(9 + 280*h*exp(-4.17*(pow(TKelvin/T0,-1./3.) - 1)))*pow(TKelvin/T0,-1./2.);
}

//Oxygen relaxation frequency
double FrO(double h, double P)
{
  return P/Patm*(24 + 4.04e4*h*(0.02+h)/(0.391+h));
}

//Water vapour saturation pressure
double psat(double TKelvin)
{
  return Patm*pow(10,-6.8346*pow(T01/TKelvin,1.261) + 4.6151);
}

//Absolute humidity
double habs(double hrel, double TKelvin, double P)
{
  return hrel*psat(TKelvin)/P;
}

//Nitrogen absorption
double Nit(double TKelvin, double h, double P, double f)
{
  double frn = FrN(TKelvin,h,P);
  return  868.6*0.1068*exp(-3352./TKelvin)*pow(TKelvin/T0,-5./2.)*frn * f*f/(frn*frn + f*f);
}

//Oxygen absorption
double Oxyg(double TKelvin, double h, double P, double f)
{
  double  fro = FrO(h,P);
  return  868.6*0.01275*exp(-2239.1/TKelvin)*pow(TKelvin/T0,-5./2.)*fro * f*f/(fro*fro+f*f);
}

/* //Total air absorption, this function is not needed actually
double AlphaAir(double TKelvin, double h, double P, double f)
{
  return classical(TKelvin,f) + Nit(TKelvin,h,P,f) + Oxyg(TKelvin,h,P,f);
}
*/

//Classical absorption, ROOT style c++ function
double CFuncAir(const double *x, const double *pars)
{
  double f = x[0];
  double TKelvin = pars[0] + T01;  //convert from celsius to kelvin
  
  return classical(TKelvin, f);
}

//Nitrogen absorption, ROOT style c++ function
double NFuncAir(const double *x, const double *pars)
{
  double f = x[0];
  double TKelvin = pars[0] + T01;  //convert from celsius to kelvin
  double HRel    = pars[1];        //Relative humidity
  double P       = pars[2] * Patm; //convert from atmosphere to pascal
  
  return Nit(TKelvin, habs(HRel,TKelvin,P), P, f);
}

//Oxygen absorption, ROOT style c++ function
double OFuncAir(const double *x, const double *pars)
{
  double f = x[0];
  double TKelvin = pars[0] + T01;  //convert from celsius to kelvin
  double HRel    = pars[1];        //Relative humidity
  double P       = pars[2] * Patm; //convert from atmosphere to pascal
  
  return Oxyg(TKelvin, habs(HRel,TKelvin,P), P, f);
}

//Total air absorption, ROOT style c++ function
double AFuncAir(const double *x, const double *pars)
{
  return CFuncAir(x,pars) + NFuncAir(x,pars) + OFuncAir(x,pars);
}

//Four-parameter approximation of air absorption
double Approx4Air(const double *x, const double *pars)
{
  double a1 = pars[3];
  double n1 = pars[4];
  double n2 = pars[5];
  double n3 = pars[6];
  double n4 = pars[7];
  double flim1 = pars[8];
  double flim2 = pars[9];
  double flim3 = pars[10];
  //if(flim1>flim2) swap(flim1,flim2);
  //if(flim2>flim3) swap(flim2,flim3);
  //if(flim1>flim2) swap(flim1,flim2);
  double a2 = a1*pow(flim1,n1-n2);
  double a3 = a2*pow(flim2,n2-n3);
  double a4 = a3*pow(flim3,n3-n4);
  double f = x[0];
  if(f<flim1)      return a1*pow(f,n1);
  else if(f<flim2) return a2*pow(f,n2);
  else if(f<flim3) return a3*pow(f,n3);
  else             return a4*pow(f,n4);
}

//Five-parameter approximation of air absorption
double Approx5Air(const double *x, const double *pars)
{
  double a1 = pars[3];
  double n1 = pars[4];
  double n2 = pars[5];
  double n3 = pars[6];
  double n4 = pars[7];
  double n5 = pars[8];
  double flim1 = pars[9];
  double flim2 = pars[10];
  double flim3 = pars[11];
  double flim4 = pars[12];
  double a2 = a1*pow(flim1,n1-n2);
  double a3 = a2*pow(flim2,n2-n3);
  double a4 = a3*pow(flim3,n3-n4);
  double a5 = a4*pow(flim4,n4-n5);
  double f = x[0];
  if(f<flim1)      return a1*pow(f,n1);
  else if(f<flim2) return a2*pow(f,n2);
  else if(f<flim3) return a3*pow(f,n3);
  else if(f<flim4) return a4*pow(f,n4);
  else             return a5*pow(f,n5);
}

//Sound speed in water
double csound(double TCelsius, double Salinity, double Depth, double pH)
{
  return 1412 + 3.21*TCelsius + 1.19*Salinity + 0.0167*Depth;
}

//Water absorption frequency scale, to make formulas simpler
double FrW(double TCelsius, double Salinity, double Depth, double pH) //frequencies in kHz here!
{
  return 1/(1 - 3.83e-5*Depth + 4.9e-10*Depth*Depth);
}

//Boric acid relaxation frequency
double FrB(double TCelsius, double Salinity, double Depth, double pH) //frequencies in kHz here!
{ 
  double TKelvin = TCelsius + KelvinConversion;
  return 2.8*sqrt(Salinity/35)*pow(10,4-1245/TKelvin);
}

//Magnesium sulphate relaxation frequency
double FrM(double TCelsius, double Salinity, double Depth, double pH) //frequencies in kHz here!
{
  double TKelvin = TCelsius + KelvinConversion;
  return 8.17*pow(10,8-1990/TKelvin)/(1 + 0.0018*(Salinity - 35));
}

//Classical absorption in water
double Water(double TCelsius, double Salinity, double Depth, double pH, double f) //frequencies in kHz here!
{
  double frw = FrW(TCelsius, Salinity, Depth, pH);
  double cw = -1;
  if(TCelsius<=20) cw = 4.937e-4 - 2.590e-5*TCelsius + 9.11e-7*TCelsius*TCelsius - 1.50e-8*TCelsius*TCelsius*TCelsius;
  else             cw = 3.964e-4 - 1.146e-5*TCelsius + 1.45e-7*TCelsius*TCelsius - 6.50e-10*TCelsius*TCelsius*TCelsius;
  return cw*f*f/(frw*frw);
}

//Magnesium sulphate absorption
double Magsulf(double TCelsius, double Salinity, double Depth, double pH, double f) //frequencies in kHz here!
{
  double frm = FrM(TCelsius, Salinity, Depth, pH);
  double c = csound(TCelsius, Salinity, Depth, pH);
  return 21.44*Salinity/c*(1+0.025*TCelsius)*(1-1.37e-4*Depth+ 6.2e-9*Depth*Depth)*frm*f*f/(frm*frm+f*f);
}

//Boric acid absorption
double Boric(double TCelsius, double Salinity, double Depth, double pH, double f) //frequencies in kHz here!
{
  double frb = FrB(TCelsius, Salinity, Depth, pH);
  double c = csound(TCelsius, Salinity, Depth, pH);
  return (8.86/c)*pow(10,0.78*pH-5)*frb*f*f/(frb*frb+f*f);
}

/* //Total water absorption, this function is not needed actually
double AlphaWater(double TCelsius, double Salinity, double Depth, double pH, double f) //frequencies in kHz here!
{
  return Water(TCelsius,Salinity,Depth,pH,f) + Magsulf(TCelsius,Salinity,Depth,pH,f) + Boric(TCelsius,Salinity,Depth,pH,f);
}
*/

//Classical absorption in water, ROOT style c++ function
double WFuncWater(const double *x, const double *pars)
{
  double f = x[0]; //in Hz
  double TCelsius = pars[0];
  double Salinity = pars[1];
  double Depth    = pars[2];
  double pH       = pars[3];
  
  return Water(TCelsius,Salinity,Depth,pH,f/1000);
}

//Boric acid absorption in water, ROOT style c++ function
double BFuncWater(const double *x, const double *pars)
{
  double f = x[0]; //in Hz
  double TCelsius = pars[0];
  double Salinity = pars[1];
  double Depth    = pars[2];
  double pH       = pars[3];
  
  return Boric(TCelsius,Salinity,Depth,pH,f/1000);
}

//Magnesium sulphate absorption in water, ROOT style c++ function
double MFuncWater(const double *x, const double *pars)
{
  double f = x[0]; //in Hz
  double TCelsius = pars[0];
  double Salinity = pars[1];
  double Depth    = pars[2];
  double pH       = pars[3];
  
  return Magsulf(TCelsius,Salinity,Depth,pH,f/1000);
}

//Total absorption in water, ROOT style c++ function
double AFuncWater(const double *x, const double *pars)
{
  return WFuncWater(x,pars) + BFuncWater(x,pars) + MFuncWater(x,pars);
}

//Water absorption approximation
double ApproxWater(const double *x, const double *pars)
{
  double a1 = pars[4];
  double n1 = pars[5];
  double n2 = pars[6];
  double n3 = pars[7];
  double n4 = pars[8];
  double n5 = pars[9];
  double flim1 = pars[10];
  double flim2 = pars[11];
  double flim3 = pars[12];
  double flim4 = pars[13];
  double a2 = a1*pow(flim1,n1-n2);
  double a3 = a2*pow(flim2,n2-n3);
  double a4 = a3*pow(flim3,n3-n4);
  double a5 = a4*pow(flim4,n4-n5);
  double f = x[0];
  if(f<flim1)      return a1*pow(f,n1);
  else if(f<flim2) return a2*pow(f,n2);
  else if(f<flim3) return a3*pow(f,n3);
  else if(f<flim4) return a4*pow(f,n4);
  else             return a5*pow(f,n5);
}