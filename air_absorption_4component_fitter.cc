#include <iostream>
#include <fstream>
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TLatex.h"
#include "TStyle.h"
#include "base_formulas_original.h"

using namespace std;

const int NPARS = 14;

const double frmin = 20.;      //min frequency in Hz
const double frmax = 20000.;  //max frequency in Hz
const int Nf = 100;            //number of frequency steps

double AFunc(const double *x, const double *pars)
{
  return AFuncAir(x,pars);
}

double Approx(const double *x, const double *pars)
{
  return Approx4Air(x,pars);
}

double ApproxDiff(const double *x, const double *pars)
{
  return (AFunc(x,pars)-Approx(x,pars)) / AFunc(x,pars);
}

double FuncToMin(const double *pars)
{
  double d = log(frmax/frmin)/Nf;        //frequency step inverse exponent
  double chi2 = 0;
  for(int i=0;i<Nf;i++)
  {
    double f = frmin*exp(i*d); // running logarithmically from 20 Hz to 20000 kHz
    double fullvalue = AFunc(&f,pars);
    double approx = Approx(&f,pars);
    double chi = (fullvalue-approx)/fullvalue;
    chi2 += chi*chi;
  }
  return chi2;
}


int main(int argc, char *argv[])
{
  double TCelsius = 20;
  double Humidity = 30;
  double Pressure = 1;
  if(argc>1) TCelsius = atof(argv[1]);
  if(argc>2) Humidity = atof(argv[2]);
  // Choose method upon creation between:
  // kMigrad, kSimplex, kCombined, 
  // kScan, kFumili
  ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
  
  cout << endl;
  cout << "***************************************************************************" << endl;
  cout << "*********** Running minimization for T = " << TCelsius << " C and H = " << Humidity << "% *******************" << endl;
  cout << "***************************************************************************" << endl;

  min.SetMaxFunctionCalls(1000000);
  min.SetMaxIterations(100000);
  min.SetTolerance(0.001);

  ROOT::Math::Functor ftr(&FuncToMin,11); 

  min.SetFunction(ftr);

  // Set the variables to be minimized
  min.SetFixedVariable(0, "TCelsius", TCelsius); //T in Celsius
  min.SetFixedVariable(1, "Humidity", Humidity); //H in percent
  min.SetFixedVariable(2, "Pressure", Pressure);  //P in atm
  
  double PPascal = Pressure * Patm;
  double TKelvin = TCelsius + T01;
  double hvalue = habs(Humidity,TKelvin,PPascal);
  double fN = FrN(TKelvin,hvalue,PPascal);
  double fO = FrO(hvalue,PPascal);
  double fC = 1; //fC does not appear in the formulas, but can be understood as cC -> cC/fC^2, for symmetry reasons
  cerr << "N and O frequencies: " << fN << "\t" << fO << endl;
  //double flim1 = fN/2; //4.2*fN*hvalue; //50?
  //double flim3 = fO/2; //0.57*fO*hvalue; //2000?
  //double flim2 = fN*2; //0.12*sqrt(flim1*flim1+flim3*flim3); //200?
  
  //IDEA: TRY TO FIX flim1 AND flim2 BASED ON THE FORMULA (THE INFLEXION POINT)
  double cN = 2*NFuncAir(&fN,min.X()); //cN=2*NFunc(fN), as then f^2/(f^2+fN^2) = 1/2
  double cO = 2*OFuncAir(&fO,min.X()); //cO=2*OFunc(fO), as then f^2/(f^2+fO^2) = 1/2
  double cC =   CFuncAir(&fC,min.X()); //cO=CFunc(fC), as then f^2=1
  
  double fiNO = sqrt(fN*fO)*pow(cN/cO+fN*fN/fO/fO,0.25); //f_i for transition from N to O
  double fiOC = sqrt(fO*fC)*pow(cO/cC+fO*fO/fC/fC,0.25); //f_i for transition from O to classical
  
  double nNO = 4*cO*fiNO*fiNO / (cO*fiNO*fiNO + cN*fO*fO + cO*fN*fN);
  double nOC = 4*cC*fiOC*fiOC / (cC*fiOC*fiOC + cO*fC*fC + cC*fO*fO);
  
  double ciNO = AFunc(&fiNO,min.X());
  double ciOC = AFunc(&fiOC,min.X());
  
  double flim1 = fiNO*pow(cN/ciNO*fiNO*fiNO/fN/fN,1.0/(nNO-2));
  double flim2 = fiNO*pow(cO/ciNO*fiNO*fiNO/fO/fO,1.0/(nNO-2));
  double flim3 = fiOC*pow(cO/ciOC*fiOC*fiOC/fO/fO,1.0/(nOC-2));
  double flim4 = fiOC*pow(cC/ciOC*fiOC*fiOC/fC/fC,1.0/(nOC-2)); //not needed
  
  //Correction based on Excel calculation
  flim1 = 0.92306*flim1 - 8.2963E-05*flim1*flim1; 
  flim2 = 0.90727*flim2 - 1.5649E-05*flim2*flim2;
  if(flim3<0 || flim3>20000) flim3 = 22222;
  else flim3 = 0.89696*flim3 + 5.3452E-07*flim3*flim3;
  //former values from previous Excel calculation
  //flim1 = 0.92687*flim1 - 1.1196e-04*flim1*flim1; 
  //flim2 = 0.90990*flim2 - 1.7660e-05*flim2*flim2;
  //if(flim3<0 || flim3>20000) flim3 = 22222;
  //else flim3 = 1.09470*flim3 - 1.9820e-05*flim3*flim3;
  
  cerr << TCelsius << "\t" << Humidity << "\t cN, cO, cC: \t" << cN << "\t" << cO << "\t" << cC << endl;
  cerr << TCelsius << "\t" << Humidity << "\t fN, fO, fC: \t" << fN << "\t" << fO << "\t" << fC << endl;
  cerr << TCelsius << "\t" << Humidity << "\t nNO, nOC: \t" << nNO << "\t" << nOC << endl;
  cerr << TCelsius << "\t" << Humidity << "\t fiNO, fiOC: \t" << fiNO << "\t" << fiOC << endl;
  cerr << TCelsius << "\t" << Humidity << "\t flim1, flim2, flim3, flim4: \t" << flim1 << "\t" << flim2 << "\t" << flim3 << "\t" << flim4 << endl;
  
  double n1 = 2;
  double n2 = nNO;
  double n3 = 2;
  double n4 = nOC;
  
  //Correction based on Excel calculation
  n1 = 1.763454532 + 0.003249811*TCelsius + 0.001091377*Humidity;
  n2 = 1.6968*n2 - 0.537*n2*n2;
  n3 = 1.671471966 + 0.000832291*TCelsius + 0.001112712*Humidity;
  n4 = 1.7413*n4;
  
  //double a1 = cN/pow(fN,n1); //this was the original
  double a1temp = 1.96005E-05*exp(-0.0159*TCelsius)*pow(Humidity,-0.4);
  double a1 = -0.83774*a1temp + 7.52E+05*a1temp*a1temp + 2.35E+10*a1temp*a1temp*a1temp;
  
  min.SetLowerLimitedVariable(3, "a1",    a1,     1e-8, 0);
  min.SetLimitedVariable(     4, "n1",    n1,      0.1,  0.5, 2.5);
  min.SetLimitedVariable(     5, "n2",    n2,      0.1,  0.5, 2.5);
  min.SetLimitedVariable(     6, "n3",    n3,      0.1,  0.5, 2.5);
  min.SetLimitedVariable(     7, "n4",    n4 ,     0.1,  0.0, 2.5);
  min.SetLimitedVariable(     8, "flim1", flim1,   1.0,  flim1/5, flim1*5); //was 140.169097 
  min.SetLimitedVariable(     9, "flim2", flim2,   5.0,  flim2/5, flim2*5); //was 1219.436748
  min.SetLimitedVariable(     10,"flim3", flim3,   10.,  flim3/5, flim3*5); //limits were fO/5, fO*2;
  
  min.FixVariable(3); //a1
  min.FixVariable(4); //n1
  min.FixVariable(5); //n2
  min.FixVariable(6); //n3
  min.FixVariable(7); //n4
  min.FixVariable(8); //flim1
  min.FixVariable(9); //flim2
  min.FixVariable(10); //flim3
  if(flim3>20000) min.FixVariable(10);
  cout << "Function value before minimization: " << FuncToMin(min.X()) << endl;
/*
  min.Minimize(); 
  
  if(min.X()[10]>=frmax)
  {
    cout << "Redoing minimization with fixed flim3 and n4" << endl;
    min.SetVariableValue(7,1.0);
    min.FixVariable(7);
    min.SetVariableValue(10,frmax*1.1);
    min.FixVariable(10);
    min.Minimize();
  }
  
  min.PrintResults();
  min.Hesse();
*/
  //PRINTING "BY HAND" IF PrintResults DOES NOT WORK
  for(int ipar=0;ipar<11;ipar++) cerr << min.VariableName(ipar) << "\t" << min.X()[ipar] << "\t +- \t" << min.Errors()[ipar] << endl;

  TCanvas *c = new TCanvas("c","c",1024,768);
  gStyle->SetOptStat(0);
  c->SetGrid(1,1);
  TH1* frame = new TH1F("frame",Form("#alpha versus f for T=%.0f C, H = %.0f%%, P = %.1f atm",TCelsius,Humidity,Pressure),100,frmin,frmax);
  frame->GetXaxis()->SetTitle("f [Hz]");
  frame->GetYaxis()->SetTitle("#alpha [dB/100 m]");
  frame->SetMinimum(0.001);
  frame->SetMaximum(100);
  frame->Draw();
  
  TF1* alphafunc = new TF1("alphafunc",AFunc,frmin,frmax,11);
  alphafunc->SetParameters(min.X());
  alphafunc->SetLineColor(8);
  alphafunc->SetLineWidth(2);
  alphafunc->Draw("SAME");
  
  TF1* approxfunc = new TF1("approxfunc",Approx,frmin,frmax,11);
  approxfunc->SetParameters(min.X());
  approxfunc->SetLineColor(2);
  approxfunc->SetLineStyle(7);
  approxfunc->Draw("SAME");
  
  c->SetLogx(1);
  c->SetLogy(1);
  c->Modified();
  //c->Print(Form("pngplots/4comp_alphaapproxfit_T%.0f_H%.0f.png",TCelsius,Humidity));
  
  delete alphafunc;
  delete approxfunc;
  
  c->SetLogy(0);
  frame->GetYaxis()->SetTitle("#alpha relative difference");
  frame->SetMinimum(-0.15);
  frame->SetMaximum(0.15);
  frame->Draw();
  TF1* approxdiff = new TF1("approxdiff",ApproxDiff,frmin,frmax,11);
  approxdiff->SetParameters(min.X());
  approxdiff->SetLineColor(2);
  approxdiff->SetLineStyle(7);
  approxdiff->Draw("SAME");
  //c->Print(Form("pngplots/4comp_alphaapproxfitdiff_T%.0f_H%.0f.png",TCelsius,Humidity));
  
  //double avgdiff = sqrt(min.MinValue()/Nf);
  double avgdiff = sqrt(FuncToMin(min.X())/Nf);
  double maxdiff = approxdiff->GetMaximum(frmin,frmax);
  double mindiff = approxdiff->GetMinimum(frmin,frmax);
  
  cout << "Avg diff: " << avgdiff << ", max diff: +" << maxdiff << " -" << mindiff << endl;

/*
  ofstream outfile("air_absorption_4component_fitter.allfixed.out",ofstream::app);
  outfile << TCelsius << "\t" << Humidity << "\t" << avgdiff << "\t" << maxdiff << "\t" << mindiff << "\t" << min.Edm() << "\t" << min.NCalls() << "\t" << min.Status() << "\t"
        << min.X()[3] << "\t" << min.Errors()[3] << "\t"
        << min.X()[4] << "\t" << min.Errors()[4] << "\t"
        << min.X()[5] << "\t" << min.Errors()[5] << "\t"
        << min.X()[6] << "\t" << min.Errors()[6] << "\t"
        << min.X()[7] << "\t" << min.Errors()[7] << "\t"
        << min.X()[8] << "\t" << min.Errors()[8]<< "\t"
        << min.X()[9] << "\t" << min.Errors()[9]<< "\t"
        << min.X()[10] << "\t" << min.Errors()[10]
          << endl;
*/
  delete frame;
  delete approxdiff;
  delete c;

  return 0;
} 
