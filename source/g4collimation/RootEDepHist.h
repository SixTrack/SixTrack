#ifndef _RootEDepHist_h
#define _RootEDepHist_h

#include <string>
#include <vector>

#include "TH2.h"
#include "TH3.h"

class RootEDepHist
{
public:

/*
name
x bin size
y bin size
z bin size
n x-bin
n y-bin
z bins will be collimator length / bin size
*/
	RootEDepHist(std::string, int, double, double, double, double, double, double, double, double, double, double, double, double, double);

	void SetName(std::string);
	std::string GetName() const;
	TH2F* GetHist_xy() const;
	TH2F* GetHist_zx() const;
	TH2F* GetHist_zy() const;

	TH3F* GetHist3D() const;

	int GetHistogramType();

	void GenerateBinVector(bool, std::vector<double>& bins, double bin_size, size_t nbins, double LargeMaximum, double LargeStep, double gap);
private:

	std::string CollimatorName;
	int HistogramType;

	TH2F* Hist_xy;
	TH2F* Hist_zx;
	TH2F* Hist_zy;
	TH3F* Hist_xyz;
};

#endif

