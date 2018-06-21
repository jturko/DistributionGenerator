
#include <fstream>
#include <iostream>
#include <string>

class DistributionGenerator
{
  public:
    DistributionGenerator(std::string filename = "EnergyDist.dat");
    ~DistributionGenerator();

    double GenerateValue();
    TH1F * GenerateDistribution(int samples = 1000, int nbins = 100);

  private:
    std::string fFileName;
    std::vector<double> fX;
    std::vector<double> fY; 
    std::vector<double> fSlp; 
    int fN;
    double fYmax;

    TRandom3 fRandom;
};

DistributionGenerator::DistributionGenerator(std::string filename)
{
    fFileName = filename;
    // load the distribution function
    double xin, yin;
    fYmax = -1;
    std::ifstream infile(fFileName);
    if(!infile.is_open()) { std::cout << "could not open file: " << fFileName << std::endl; return; }
    while(!infile.eof()) {
        infile >> xin >> yin;
        //std::cout << " x = " << xin << " y = " << yin << std::endl;
        fX.push_back(xin);
        fY.push_back(yin);
        if(yin>fYmax) fYmax = yin;
    }
    fN = (int)fX.size();

    // calculate slope
    double dx = fX[1]-fX[0];
    fSlp.resize(fN);
    for(int i=0; i<fN-1; i++) {
        fSlp[i] = (fY[i+1]-fY[i])/dx;
    }
    
}

DistributionGenerator::~DistributionGenerator() { }

double DistributionGenerator::GenerateValue()
{
    double xrndm = 0.;
    double yrndm = 0.;
    double yinter = -1.;
    while(yrndm>yinter) {
        // randomly choose point
        xrndm = fX[0] + fRandom.Rndm()*(fX[fN-1]-fX[0]);
        yrndm = fRandom.Rndm()*fYmax;
        // find bin
        int j = fN-2;
        while( ( fX[j] > xrndm ) && (j > 0) ) j--;
        yinter = fY[j] + fSlp[j]*(xrndm-fX[j]);
    }
    return xrndm;
}

TH1F * DistributionGenerator::GenerateDistribution(int samples, int nbins) 
{
    TH1F * histo = new TH1F(Form("%s",fFileName.c_str()),Form("%s",fFileName.c_str()),nbins,fX[0],fX[fN-1]);
    for(int i=0; i<samples; i++) histo->Fill(GenerateValue());

    return histo;

}
