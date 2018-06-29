
#include <iostream>
#include <fstream>

void PlotAngularDist(int npt = 180, double n = 2.80, std::string filename = "AngularDist.dat")
{
    double lo = 0; // radians
    double hi = TMath::Pi()/2.; // radians
    TF1 * angdist = new TF1("angdist","[0]*TMath::Power(TMath::Cos(x),[1]+1.)*TMath::Sin(x)",lo,hi);
    //TF1 * angdist = new TF1("angdist","[0]*TMath::Power(TMath::Cos(x),[1])",lo,hi);
    double I_0 = 1.0;
    angdist->SetParameters(I_0,n);
    angdist->SetNpx(10000);

    angdist->Draw();

    TGraph * angdistgr = new TGraph;
    std::ofstream outfile;
    outfile.open(filename);
    double inc = (hi-lo)/double(npt);
    for(double xx=lo; xx<=hi; xx+=inc, npt++) {
        outfile << xx << "\t" << angdist->Eval(xx) << "\n";
        angdistgr->SetPoint(npt,xx,angdist->Eval(xx));
    }
    outfile.close();

    new TCanvas();
    angdistgr->Draw("a*");
}

