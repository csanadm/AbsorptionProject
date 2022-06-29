#include <iostream>
#include <fstream>
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGraph.h"
#include "water_air_simpifield_formula.h"
#include "other_materials_ai.h"

using namespace std;

const double frmin = 20.;      //min frequency in Hz
const double frmax = 20000.;  //max frequency in Hz

int main()
{
  TCanvas *c = new TCanvas("c","c",1024,768);
  gStyle->SetOptStat(0);
  c->SetGrid(1,1);
  TH1* frame = new TH1F("frame","frame",100,frmin,frmax);
  frame->GetXaxis()->SetTitle("f [Hz]");
  
  double TCelsiusAir = 20;
  double Humidity = 30;
  double Pressure = 1;
  frame->SetTitle(Form("#alpha versus f for T=%.0f C, H = %.0f%%, P = %.1f atm",TCelsiusAir,Humidity,Pressure));  
  frame->GetYaxis()->SetTitle("#alpha [dB/100 m]");
  frame->SetMinimum(0.001);
  frame->SetMaximum(100);
  frame->Draw();
  TF1* air_absorption_orig = new TF1("air_absorption_orig",AFuncAir,frmin,frmax,NParsa);
  air_absorption_orig->SetParameters(TCelsiusAir,Humidity,Pressure);
  air_absorption_orig->SetLineColor(8);
  air_absorption_orig->SetLineWidth(2);
  air_absorption_orig->Draw("SAME");
  TF1* air_absorption_func = new TF1("air_absorption_func",air_absorption,frmin,frmax,NParsa);
  air_absorption_func->SetParameters(TCelsiusAir,Humidity,Pressure);
  air_absorption_func->SetLineColor(2);
  air_absorption_func->SetLineStyle(7);
  air_absorption_func->Draw("SAME");
  c->SetLogx(1);
  c->SetLogy(1);
  c->Modified();
  c->Print(Form("air_approximation_T%.0f_H%.0f_P%.0f.png",TCelsiusAir,Humidity,Pressure));
  
  double TCelsiusWater = 8;
  double Salinity = 35;
  double Depth    = 500;
  double pH       = 8;
  frame->SetTitle(Form("#alpha versus f for T=%.0f C, S = %.0f ppt, D = %.0f m, pH = %.0f",TCelsiusWater,Salinity,Depth,pH));
  frame->GetYaxis()->SetTitle("#alpha [dB/km]");
  frame->SetMinimum(0.00001);
  frame->SetMaximum(10);
  frame->Draw();
  TF1* water_absorption_orig = new TF1("water_absorption_orig",AFuncWater,frmin,frmax,NParsw);
  water_absorption_orig->SetParameters(TCelsiusWater,Salinity,Depth,pH);
  water_absorption_orig->SetLineColor(8);
  water_absorption_orig->SetLineWidth(2);
  water_absorption_orig->Draw("SAME");
  TF1* water_absorption_func = new TF1("water_absorption_func",water_absorption,frmin,frmax,NParsw); 
  water_absorption_func->SetParameters(TCelsiusWater,Salinity,Depth,pH);
  water_absorption_func->SetLineColor(2);
  water_absorption_func->SetLineStyle(7);
  water_absorption_func->Draw("SAME");
  c->SetLogx(1);
  c->SetLogy(1);
  c->Modified();
  c->Print(Form("water_approximation_T%.0f_S%.0f_D%.0f_pH%.0f.png",TCelsiusWater,Salinity,Depth,pH));
  
  int nmat = 5; //"Smooth brickwork with flush pointing"
  frame->SetTitle(Form("#alpha versus f for %s",materialnames[nmat]));
  frame->GetYaxis()->SetTitle("#alpha for full material width");
  TGraph* othermat_data = new TGraph(NAI,frequenciesOtherMat,aivalues[nmat]);
  frame->SetMinimum(0.002);
  frame->SetMaximum(2);
  frame->Draw();
  othermat_data->SetMarkerStyle(20);
  othermat_data->SetMarkerSize(1);
  othermat_data->Draw("P");
  TF1* othermat_absorption_func = new TF1("othermat_absorption_func",AlphaOtherMaterial,frmin,frmax,1); 
  othermat_absorption_func->SetParameter(0,nmat);
  othermat_absorption_func->SetLineColor(2);
  othermat_absorption_func->SetLineStyle(7);
  othermat_absorption_func->Draw("SAME");
  c->SetLogx(1);
  c->SetLogy(1);
  c->Modified();
  c->Print(Form("othermat_approximation_nmat%d.png",nmat));
  
  delete frame;
  delete water_absorption_func;
  delete water_absorption_orig;
  delete air_absorption_func;
  delete air_absorption_orig;
  delete c;

  return 0;
}