
#include "/Users/JTurko/root_scripts/relativistic_tools.C"

#include <iostream>
#include <fstream>

// taken from Eq 8 in the paper "Muon-Induced Background Study 
// for Underground Laboratories" by D. Mei and A. Hime
TF1 * PlotEnergyDist(double dx = 10, double h = 1.585, std::string filename = "EnergyDist.dat")
{
    double elo = 1;  // GeV
    double ehi = 1e3;   // GeV
    TF1 * edist = new TF1("edist","[0]*TMath::Exp(-[1]*[2]*([3]-1.)) * TMath::Power( x+[4]*(1-TMath::Exp(-[1]*[2])) , -[3] )", elo, ehi);
    double A = 1;               // 0 - A
    if(h==1.585) A/=5.76e-11;
    double b = 0.4;             // 1 - b
                                // 2 - h
    double gamma_mu = 3.77;     // 3 - gamma_mu
    double epsilon_mu = 693;    // 4 - epsilon_mu
    edist->SetParameters(A, b, h, gamma_mu, epsilon_mu);
    edist->SetNpx(10000);

    edist->GetXaxis()->SetRangeUser(1,1e3);
    edist->GetYaxis()->SetRangeUser(1e-13,1e-10);
    new TCanvas();
    edist->Draw();
    //gPad->SetLogx(); gPad->SetLogy();

    TGraph * edistgr = new TGraph;
    std::ofstream outfile;
    outfile.open(filename);
    double inc = dx; int npt = 0;
    for(double xx=elo; xx<=ehi; xx+=inc, npt++) {
        outfile << xx << "\t" << edist->Eval(xx) << "\n";
        edistgr->SetPoint(npt,xx,edist->Eval(xx));
    }
    outfile.close();

    //new TCanvas();
    //edistgr->Draw("a*");

    return edist;
}

// eqn 3.27 from Greider2001 but momentum converted to energy 
double EnergyDistGrieder327(double * x, double * par)
{
    double energy = x[0];
    double mu_mass = 0.1056; // GeV/c^2
    double momentum = T2P(energy,mu_mass);
    double a = 2.47e-3, b = 0.4854, c = 0.3406;
    double y = a*std::pow(momentum,-b-c*std::log(momentum));
    return y/2.35e-3;
}

// eqn 3.28 from Greider2001 but momentum converted to energy 
double EnergyDistGrieder328(double * x, double * par)
{
    double energy = x[0];
    double mu_mass = 0.1056; // GeV/c^2
    double momentum = T2P(energy,mu_mass);
    double a = 3.09e-3, b = 0.5483, c = 0.3977;
    double y = a*std::pow(momentum,-b-c*std::log(momentum));
    return y/2.35e-3;
}

// eqn 3.27 from Greider2001, for the momentum though
TF1 * PlotMomentumEnergyDistGrieder327(double dx=10, std::string filename="EnergyDistGrieder.dat")
{
    double plo = 1;  // GeV/c
    double phi = 1e3;   // GeV/c
    TF1 * func = new TF1("IvsP","[0]*TMath::Power(x,-[1]-[2]*TMath::Log(x))",plo,phi);
    func->SetParameters(2.47e-3, 0.4854, 0.3406);
    func->SetNpx(10000);
 
    //new TCanvas();
    func->GetXaxis()->SetTitle("Momentum [GeV/c]");
    func->GetYaxis()->SetTitle("Intensity [s^{-1} sr^{-1} m^{-2} (GeV/c)^{-1}]");   
    //func->Draw();

    double elo = 1; // GeV
    double ehi = 1e3; // GeV
    TF1 * func2 = new TF1("IvsE",EnergyDistGrieder327,elo,ehi,0);
    func2->SetNpx(10000);

    new TCanvas();
    func2->GetXaxis()->SetTitle("Energy [GeV]");
    func2->GetYaxis()->SetTitle("Intensity [s^{-1} sr^{-1} m^{-2} (GeV/c)^{-1}]");   
    func2->Draw();    
    
    std::ofstream outfile;
    outfile.open(filename);
    double inc = dx;
    for(double xx=elo; xx<=ehi; xx+=inc) {
        outfile << xx << "\t" << func2->Eval(xx) << "\n";
    }
    outfile.close();

    return func2;
}

// eqn 3.28 from Greider2001, for the momentum though
TF1 * PlotMomentumEnergyDistGrieder328(double dx=10, std::string filename="EnergyDistGrieder.dat")
{
    double plo = 1;  // GeV/c
    double phi = 1e3;   // GeV/c
    TF1 * func = new TF1("IvsP","[0]*TMath::Power(x,-[1]-[2]*TMath::Log(x))",plo,phi);
    func->SetParameters(2.47e-3, 0.4854, 0.3406);
    func->SetNpx(10000);
 
    //new TCanvas();
    func->GetXaxis()->SetTitle("Momentum [GeV/c]");
    func->GetYaxis()->SetTitle("Intensity [s^{-1} sr^{-1} m^{-2} (GeV/c)^{-1}]");   
    //func->Draw();

    double elo = 1; // GeV
    double ehi = 1e3; // GeV
    TF1 * func2 = new TF1("IvsE",EnergyDistGrieder328,elo,ehi,0);
    func2->SetNpx(10000);

    new TCanvas();
    func2->GetXaxis()->SetTitle("Energy [GeV]");
    func2->GetYaxis()->SetTitle("Intensity [s^{-1} sr^{-1} m^{-2} (GeV/c)^{-1}]");   
    func2->Draw();    
    
    std::ofstream outfile;
    outfile.open(filename);
    double inc = dx;
    for(double xx=elo; xx<=ehi; xx+=inc) {
        outfile << xx << "\t" << func2->Eval(xx) << "\n";
    }
    outfile.close();

    return func2;
}
