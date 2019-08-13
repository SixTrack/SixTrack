#include <iostream>

#include "RootEDepHist.h"
#include "TH2.h"
#include "TH3.h"

/*
name
x bin size
y bin size
z bin size
n x-bin
n y-bin
z bins will be collimator length / bin size
*/
//RootEDepHist::RootEDepHist(std::string name, int type, double CollimatorLength, double halfgap, double x_bin, double y_bin, double z_bin, size_t nx, size_t ny) : CollimatorName(name)
RootEDepHist::RootEDepHist(std::string name, int type, double halfgap, double xmax, double ymax, double zmax, double xbin, double ybin, double zbin, double xbigmax, double ybigmax, double zbigmax, double xbigstep, double ybigstep, double zbigstep) : CollimatorName(name), HistogramType(type)
{
//CollimatorLength should be in mm
	size_t nz = ceil(zmax / zbin);
	size_t ny = ceil(ymax / ybin);
	size_t nx = ceil(xmax / xbin);

	//extra +1 for boundaries, x has to have the "jaw gap" bin

	std::vector<double> zbins;
	GenerateBinVector(false, zbins, zbin, nz, zbigmax, zbigstep, 0.0);
	nz = zbins.size()-1;

	std::vector<double> ybins;
	GenerateBinVector(true, ybins, ybin, ny, ybigmax, ybigstep, 0.0);
	ny = ybins.size()-1;

	std::vector<double> xbins;
	GenerateBinVector(true, xbins, xbin, nx, xbigmax, xbigstep, halfgap-1e-8);
	nx = xbins.size()-1;

	std::cout << "GEANT4> Adding " << type << "D energy deposition histogram for \"" << name << "\" with n bins (x,y,z) (" << nx << "," << ny << "," << nz << ")" << std::endl;
	std::cout << "GEANT4> x: " << xbins.at(0) << " -> " << xbins.at(xbins.size()-1) << std::endl;
	std::cout << "GEANT4> y: " << ybins.at(0) << " -> " << ybins.at(ybins.size()-1) << std::endl;
	std::cout << "GEANT4> z: " << zbins.at(0) << " -> " << zbins.at(zbins.size()-1) << std::endl;

	if(type == 2)
	{
		std::string n_xy = "edep_" + name + "_xy";
		std::string n_zx = "edep_" + name + "_zx";
		std::string n_zy = "edep_" + name + "_zy";

		Hist_xy = new TH2F(n_xy.c_str(), name.c_str(), nx, xbins.data(), ny, ybins.data());
		Hist_zx = new TH2F(n_zx.c_str(), name.c_str(), nz, zbins.data(), nx, xbins.data());
		Hist_zy = new TH2F(n_zy.c_str(), name.c_str(), nz, zbins.data(), ny, ybins.data());
	}
	else if(type == 3)
	{
		std::string n_xyz = "edep_" + name + "_xyz";
		Hist_xyz = new TH3F(n_xyz.c_str(), name.c_str(), nx, xbins.data(), ny, ybins.data(), nz, zbins.data());
	}
	else
	{
		//Bad input?
		std::cerr << "GEANT4> ERROR: Neither a 2d or 3d histogram was requested for \"" << name << "\" - exiting!" << std::endl;
		exit(EXIT_FAILURE);
	}
}

void RootEDepHist::SetName(std::string in)
{
	CollimatorName = in;
}

std::string RootEDepHist::GetName() const
{
	return CollimatorName;
}

TH2F* RootEDepHist::GetHist_xy() const
{
	return Hist_xy;
}

TH2F* RootEDepHist::GetHist_zx() const
{
	return Hist_zx;
}

TH2F* RootEDepHist::GetHist_zy() const
{
	return Hist_zy;
}

TH3F* RootEDepHist::GetHist3D() const
{
	return Hist_xyz;
}

void RootEDepHist::GenerateBinVector(bool NegativeBins, std::vector<double>& bins, double bin_size, size_t nbins, double LargeMaximum, double LargeStep, double gap)
{
	std::vector<double> prebin;

	//Fill the core region of interest
	for(size_t n=0; n <= nbins; n++)
	{
		prebin.push_back((n*bin_size) + gap);
	}

	//Fill up to a defined maximum at lower resolution (LargeMaximum/LargeStep)
	//Jump up to the next integer in mm
	double binStart = ceil(prebin.at(prebin.size()-1));
	//Check we are not adding the same bin twice
	if(binStart == prebin.at(prebin.size()-1))
	{
		binStart+=LargeStep;
	}

	double ThisBin = binStart;
	while(ThisBin <= LargeMaximum)
	{
		prebin.push_back(ThisBin);
		ThisBin+=LargeStep;
	}

	if(NegativeBins)
	{
		for(size_t n=prebin.size(); n-->0;)
		{
			bins.push_back(-prebin.at(n));
		}
	}

	for(size_t n=0; n<prebin.size(); n++)
	{
		bins.push_back(prebin.at(n));
	}
}

int RootEDepHist::GetHistogramType()
{
	return  HistogramType;
}
