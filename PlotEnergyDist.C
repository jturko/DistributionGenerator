
#include <iostream>
#include <fstream>

// taken from Eq 8 in the paper "Muon-Induced Background Study 
// for Underground Laboratories" by D. Mei and A. Hime
void PlotEnergyDist(int dx = 10)
{
    double elo = 0;  // GeV
    double ehi = 1e4;   // GeV
    TF1 * edist = new TF1("edist","[0]*TMath::Exp(-[1]*[2]*([3]-1.)) * TMath::Power( x+[4]*(1-TMath::Exp(-[1]*[2])) , -[3] )", elo, ehi);
    double A = 1;               // 0 - A
    double b = 0.4;             // 1 - b
    double h = 1.585;           // 2 - h
    double gamma_mu = 3.77;     // 3 - gamma_mu
    double epsilon_mu = 693;    // 4 - epsilon_mu
    edist->SetParameters(A, b, h, gamma_mu, epsilon_mu);
    edist->SetNpx(10000);

    edist->GetXaxis()->SetRangeUser(1,1e3);
    edist->GetYaxis()->SetRangeUser(1e-13,1e-10);
    edist->Draw();
    //gPad->SetLogx(); gPad->SetLogy();

    TGraph * edistgr = new TGraph;
    std::ofstream outfile;
    outfile.open("EnergyDist.dat");
    double inc = dx; int npt = 0;
    for(double xx=elo; xx<=ehi; xx+=inc, npt++) {
        outfile << xx << "\t" << edist->Eval(xx) << "\n";
        edistgr->SetPoint(npt,xx,edist->Eval(xx));
    }
    outfile.close();

    new TCanvas();
    edistgr->Draw("a*");
}

// taken from Eq 1 in the paper "The Cosmic Ray Muon Flux 
// at WIPP" by E, Esch et al
void EschEnergyDist()
{
    double a = 3.09e-3; // cm^-2 s^-1 sr^-1
    double b = 0.5483;
    double c = 0.3977;

    TF1 * pfunc = new TF1("pfunc","[0]*TMath::Power(x,-[1]-[2]*TMath::Log(x))",0.1,1e4);
    pfunc->SetParameter(0,a);
    pfunc->SetParameter(1,b);
    pfunc->SetParameter(2,c);
    pfunc->GetXaxis()->SetTitle("Muon momentum [GeV/c]");
    pfunc->GetYaxis()->SetTitle("Relative intensity [arb. units]");
    pfunc->SetNpx(10000);
    
    new TCanvas();
    pfunc->Draw("");    
    gPad->SetLogy();
    gPad->SetLogx();
    
    // ([3]*TMath::Sqrt(TMath::Power((x/[3])+1,2)-1)) // this is the p(T) function
    double m = 0.10566; // GeV c^-2
    TF1 * tfunc = new TF1("tfunc","[0]*TMath::Power( ([3]*TMath::Sqrt(TMath::Power((x/[3])+1.,2.)-1.)) , -[1]-[2]*TMath::Log( ([3]*TMath::Sqrt(TMath::Power((x/[3])+1.,2.)-1.)) ) )",0.1,1e4);
    tfunc->SetParameter(0,a);
    tfunc->SetParameter(1,b);
    tfunc->SetParameter(2,c);
    tfunc->SetParameter(3,m);
    tfunc->GetXaxis()->SetTitle("Muon energy [GeV]");
    tfunc->GetYaxis()->SetTitle("Relative intensity [arb. units]");
    tfunc->SetNpx(10000);

    //new TCanvas();
    tfunc->SetLineColor(kBlack);
    tfunc->Draw("");    
    gPad->SetLogy();
    gPad->SetLogx();

    TF1 * T2p = new TF1("T2p","([0]*TMath::Sqrt(TMath::Power((x/[0])+1.,2.)-1.))",0,1e4);
    T2p->SetParameter(0,m);
    
}
