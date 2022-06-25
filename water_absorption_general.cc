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
  return AFuncWater(x,pars);
}

double Approx(const double *x, const double *pars)
{
  return ApproxWater(x,pars);
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
  double TCelsius = 8;
  double Salinity = 35;
  double Depth    = 50;
  double pH       = 8;
  if(argc>1) TCelsius = atof(argv[1]);
  if(argc>2) Salinity = atof(argv[2]);
  if(argc>3) Depth    = atof(argv[3]);
  if(argc>4) pH       = atof(argv[4]);
  // Choose method upon creation between:
  // kMigrad, kSimplex, kCombined, 
  // kScan, kFumili
  ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
  
  cout << endl;
  cout << "***************************************************************************" << endl;
  cout << "*********** Running minimization for T = " << TCelsius << " C and S = " << Salinity << " ppt *****************" << endl;
  cout << "***************************************************************************" << endl;

  min.SetMaxFunctionCalls(1000000);
  min.SetMaxIterations(100000);
  min.SetTolerance(0.001);

  ROOT::Math::Functor ftr(&FuncToMin,NPARS); 

  min.SetFunction(ftr);

  // Set the variables to be minimized
  min.SetFixedVariable(0, "TCelsius", TCelsius); //T in Celsius
  min.SetFixedVariable(1, "Salinity", Salinity); //S in ppt
  min.SetFixedVariable(2, "Depth",    Depth);    //D in m
  min.SetFixedVariable(3, "pH",       pH);       //D in m
  
  double fW = FrW(TCelsius, Salinity, Depth, pH)*1000; //in Hz
  double fB = FrB(TCelsius, Salinity, Depth, pH)*1000; //in Hz
  double fM = FrM(TCelsius, Salinity, Depth, pH)*1000; //in Hz
  cerr << "W, B, M frequencies (in kHz): " << fW << "\t" << fB << "\t" << fM << endl;

  
  //IDEA: TRY TO FIX flim1 AND flim2 BASED ON THE FORMULA (THE INFLEXION POINT)
  double cW =   WFuncWater(&fW,min.X()); //c=  Func(f0), as then f^2/f0^2       = 1
  double cB = 2*BFuncWater(&fB,min.X()); //c=2*Func(f0), as then f^2/(f^2+f0^2) = 1/2
  double cM = 2*MFuncWater(&fM,min.X()); //c=2*Func(f0), as then f^2/(f^2+f0^2) = 1/2
  
  double fiBM = sqrt(fB*fM)*pow(cB/cM+fB*fB/fM/fM,0.25); //f_i for transition from B to M
  double fiMW = sqrt(fM*fW)*pow(cM/cW+fM*fM/fW/fW,0.25); //f_i for transition from M to W
  
  double nBM = 4*cM*fiBM*fiBM / (cM*fiBM*fiBM + cB*fM*fM + cM*fB*fB);
  double nMW = 4*cW*fiMW*fiMW / (cW*fiMW*fiMW + cM*fW*fW + cW*fM*fM);
  
  double ciBM = AFunc(&fiBM,min.X());
  double ciMW = AFunc(&fiMW,min.X());
  
  double flim1 = fiBM*pow(cB/ciBM*fiBM*fiBM/fB/fB,1.0/(nBM-2));
  double flim2 = fiBM*pow(cM/ciBM*fiBM*fiBM/fM/fM,1.0/(nBM-2));
  double flim3 = fiMW*pow(cM/ciMW*fiMW*fiMW/fM/fM,1.0/(nMW-2));
  double flim4 = fiMW*pow(cW/ciMW*fiMW*fiMW/fW/fW,1.0/(nMW-2));
  
  cerr << TCelsius << "\t" << Salinity << "\t" << Depth << "\t cB, cM, cW: \t" << cB << "\t" << cM << "\t" << cW << endl;
  cerr << TCelsius << "\t" << Salinity << "\t" << Depth << "\t fB, fM, fW: \t" << fB << "\t" << fM << "\t" << fW << endl;
  cerr << TCelsius << "\t" << Salinity << "\t" << Depth << "\t nBM, nBW: \t" << nBM << "\t" << nMW << endl;
  cerr << TCelsius << "\t" << Salinity << "\t" << Depth << "\t fiBM, fiMW: \t" << fiBM << "\t" << fiMW << endl;
  cerr << TCelsius << "\t" << Salinity << "\t" << Depth << "\t flim1, flim2, flim3, flim4: \t" << flim1 << "\t" << flim2 << "\t" << flim3 << "\t" << flim4 << endl;
  
  //Correction based on Excel calculation
  flim1 = 0.8282*flim1 + 25.818;
  flim2 = 0.7481*flim2 + 647.67;
  
  double n1 = 2;
  double n2 = nBM;
  double n3 = 2;
  double n4 = nMW;
  double n5 = 2;
  
  //Correction based on Excel calculation
  n1 = 1.9093 + 0.000968093*TCelsius + 0.000928211*Salinity + 4.61102e-7*TCelsius*Salinity -  9.12418e-7*Depth - 4.34584e-10*Salinity*Depth;
  n2 = 0.8716*n2 + 0.2902;
  n3 = 1.51232 - 0.000024935*Depth + 0.00661563*Salinity +  0.00219034*TCelsius - 0.000312489*TCelsius*TCelsius;
  
  double a1 = cB/pow(fB,n1);
  //Correction based on Excel calculation
  a1 = 0.0000002651151*(1 - 0.040785676*TCelsius+0.0006424371*TCelsius*TCelsius)*(1-0.008556518*Salinity)*(1-0.000012658*Depth);
  
  min.SetLowerLimitedVariable(4, "a1",    a1,     1e-8, 0);
  min.SetLimitedVariable(     5, "n1",    n1,      0.1,  0.5, 2.5);
  min.SetLimitedVariable(     6, "n2",    n2,      0.1,  0.5, 2.5);
  min.SetLimitedVariable(     7, "n3",    n3,      0.1,  0.5, 2.5);
  min.SetLimitedVariable(     8, "n4",    n4 ,     0.1,  0.0, 2.5);
  min.SetLimitedVariable(     9, "n5",    n5 ,     0.1,  0.0, 2.5);
  min.SetLimitedVariable(     10, "flim1", flim1,   1.0,  flim1/5, flim1*5);
  min.SetLimitedVariable(     11, "flim2", flim2,   5.0,  flim2/5, flim2*5);
  min.SetLimitedVariable(     12, "flim3", flim3,   10.,  flim3/5, flim3*5);
  min.SetLimitedVariable(     13, "flim4", flim4,   10.,  flim4/5, flim4*5); 
  
  min.FixVariable(4); //a1
  min.FixVariable(5); //n1 -> now fixed to corrected value above
  min.FixVariable(6); //n2 -> now fixed to corrected value above
  min.FixVariable(7); //n3 -> now fixed to corrected value above
  min.FixVariable(8); //n4 -> plays no role if flim3,4 above 20 kHz
  min.FixVariable(9); //n5 -> plays no role if flim3,4 above 20 kHz
  min.FixVariable(10); //flim1 -> now fixed to corrected value above
  min.FixVariable(11); //flim2 -> now fixed to corrected value above
  min.FixVariable(12); //flim3 -> plays no role if flim3,4 above 20 kHz
  min.FixVariable(13); //flim4 -> plays no role if flim3,4 above 20 kHz
  if(flim1<20) min.FixVariable(10);
  if(flim3>20000) min.FixVariable(12);
  if(flim4>20000) min.FixVariable(13);
  cout << "Function value before minimization: " << FuncToMin(min.X()) << endl;

  //min.Minimize(); 
  //min.PrintResults();
  //min.Hesse();

  //PRINTING "BY HAND" IF PrintResults DOES NOT WORK
  for(int ipar=0;ipar<NPARS;ipar++) cerr << min.VariableName(ipar) << "\t" << min.X()[ipar] << "\t +- \t" << min.Errors()[ipar] << endl;

  TCanvas *c = new TCanvas("c","c",1024,768);
  gStyle->SetOptStat(0);
  c->SetGrid(1,1);
  TH1* frame = new TH1F("frame",Form("#alpha versus f for T=%.0f C, S = %.0f ppt, D = %.0f m, pH = %.0f",TCelsius,Salinity,Depth,pH),100,frmin,frmax);
  frame->GetXaxis()->SetTitle("f [Hz]");
  frame->GetYaxis()->SetTitle("#alpha [dB/km]");
  frame->SetMinimum(0.00001);
  frame->SetMaximum(10);
  frame->Draw();
  
  TF1* alphafunc = new TF1("alphafunc",AFunc,frmin,frmax,NPARS);
  alphafunc->SetParameters(min.X());
  alphafunc->SetLineColor(8);
  alphafunc->SetLineWidth(2);
  alphafunc->Draw("SAME");
  
  TF1* approxfunc = new TF1("approxfunc",Approx,frmin,frmax,NPARS);
  approxfunc->SetParameters(min.X());
  approxfunc->SetLineColor(2);
  approxfunc->SetLineStyle(7);
  approxfunc->Draw("SAME");
  
  c->SetLogx(1);
  c->SetLogy(1);
  c->Modified();
  c->Print(Form("pngplots/water_absorption_fit_T%.0f_S%.0f_D%.0f.png",TCelsius,Salinity,Depth));
  
  delete alphafunc;
  delete approxfunc;
  
  c->SetLogy(0);
  frame->GetYaxis()->SetTitle("#alpha relative difference");
  frame->SetMinimum(-0.5);
  frame->SetMaximum(0.5);
  frame->Draw();
  TF1* approxdiff = new TF1("approxdiff",ApproxDiff,frmin,frmax,NPARS);
  approxdiff->SetParameters(min.X());
  approxdiff->SetLineColor(2);
  approxdiff->SetLineStyle(7);
  approxdiff->Draw("SAME");
  c->Print(Form("pngplots/water_absorption_fitdiff_T%.0f_S%.0f_D%.0f.png",TCelsius,Salinity,Depth));
  
  //double avgdiff = sqrt(min.MinValue()/Nf);
  double avgdiff = sqrt(FuncToMin(min.X())/Nf);
  double maxdiff = approxdiff->GetMaximum(frmin,frmax);
  double mindiff = approxdiff->GetMinimum(frmin,frmax);
  
  cout << "Avg diff: " << avgdiff << ", max diff: +" << maxdiff << " -" << mindiff << endl;

  ofstream outfile("water_absorption_general.out",ofstream::app);
  outfile << TCelsius << "\t" << Salinity << "\t" << Depth << "\t" << avgdiff << "\t" << maxdiff << "\t" << mindiff << "\t" << min.Edm() << "\t" << min.NCalls() << "\t" << min.Status() << "\t"
        << min.X()[4] << "\t" << min.Errors()[4] << "\t"
        << min.X()[5] << "\t" << min.Errors()[5] << "\t"
        << min.X()[6] << "\t" << min.Errors()[6] << "\t"
        << min.X()[7] << "\t" << min.Errors()[7] << "\t"
        << min.X()[8] << "\t" << min.Errors()[8] << "\t"
        << min.X()[9] << "\t" << min.Errors()[9] << "\t"
        << min.X()[10] << "\t" << min.Errors()[10] << "\t"
        << min.X()[11] << "\t" << min.Errors()[11] << "\t"
        << min.X()[12] << "\t" << min.Errors()[12] << "\t"
        << min.X()[13] << "\t" << min.Errors()[13]
          << endl;
 
  delete frame;
  delete approxdiff;
  delete c;

  return 0;
} 
