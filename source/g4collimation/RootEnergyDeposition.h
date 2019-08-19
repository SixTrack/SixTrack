#ifndef _RootEnergyDeposition_h
#define _RootEnergyDeposition_h

#include <map>

#include "CollimationJawHit.h"
#include "CollimationStorage.h"
#include "RootEnergyDeposition.h"
#include "RootEDepHist.h"

#include "TH2.h"

class RootEnergyDeposition
{
public:

	RootEnergyDeposition(std::vector<CollimationEnergyDeposition>);

	void SetCollimator(std::string, double, double);

	void Process(CollimationJawHitsCollection*);

	void WriteHistogramsToFile();
	void WriteHist();

	void FillHit(TH2F* hist, double xs, double xe, double ys, double ye, double edep);
	void FillHit(TH3F* hist, double xs, double xe, double ys, double ye, double zs, double ze, double edep);

private:

	RootEDepHist* ThisEntry;
	std::map<std::string,RootEDepHist*> HistogramMap;

	std::vector<CollimationEnergyDeposition> EnergyDepositionConfiguration;

};

#endif

