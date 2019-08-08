#include <algorithm>
#include <map>

#include "RootEnergyDeposition.h"
#include "CollimationJawHit.h"

#include "TROOT.h"
#include "TCanvas.h"
#include "TH2.h"

RootEnergyDeposition::RootEnergyDeposition(std::vector<CollimationEnergyDeposition> EnergyDepositionConfiguration_in)
{
	EnergyDepositionConfiguration = EnergyDepositionConfiguration_in;
}

void RootEnergyDeposition::SetCollimator(std::string name, double length, double hgap)
{
	std::map<std::string,RootEDepHist*>::iterator itr = HistogramMap.begin();

	double xmax,ymax,zmax;
	double xstep,ystep,zstep;
	double xbigmax,ybigmax,zbigmax;
	double xbigstep,ybigstep,zbigstep;
	int type;

	itr = HistogramMap.find(name);
	if(itr != HistogramMap.end())
	{
		ThisEntry = itr->second;
	}
	else
	{
		//Create
		//Find settings for this collimatior histogram
		std::vector<CollimationEnergyDeposition>::iterator itf = std::find_if(EnergyDepositionConfiguration.begin(), EnergyDepositionConfiguration.end(), find_edep_name(name));
		if(itf != EnergyDepositionConfiguration.end())
		{
			//Found the collimator
			//Extract values
			type=itf->type;

			xstep=itf->xstep;
			ystep=itf->ystep;
			zstep=itf->zstep;

			xmax=itf->xmax;
			ymax=itf->ymax;
			zmax=itf->zmax;

			xbigstep=itf->xbigstep;
			ybigstep=itf->ybigstep;
			zbigstep=itf->zbigstep;

			xbigmax=itf->xbigmax;
			ybigmax=itf->ybigmax;
			zbigmax=itf->zbigmax;

			if(itf->zmax == -1) zmax=length;
			if(itf->zbigmax == -1) zbigmax=length;
		}
		else
		{
			//Not found, see if the ALL entry exists
			itf = std::find_if(EnergyDepositionConfiguration.begin(), EnergyDepositionConfiguration.end(), find_edep_name("ALL"));
			if(itf!=EnergyDepositionConfiguration.end())
			{
				//Found global/fallback values
				type=itf->type;

				xstep=itf->xstep;
				ystep=itf->ystep;
				zstep=itf->zstep;

				xmax=itf->xmax;
				ymax=itf->ymax;
				zmax=itf->zmax;

				xbigstep=itf->xbigstep;
				ybigstep=itf->ybigstep;
				zbigstep=itf->zbigstep;

				xbigmax=itf->xbigmax;
				ybigmax=itf->ybigmax;
				zbigmax=itf->zbigmax;

				if(itf->zmax == -1) zmax=length;
				if(itf->zbigmax == -1) zbigmax=length;
			}
			else
			{
				//Did not find anything?
				//Skip this collimator
				ThisEntry = nullptr;
				return;
				//std::cerr << "GEANT4> ERROR: Collimation energy deposition was requested but no histogram settings were found!" << std::endl;
				//exit(EXIT_FAILURE);
			}
		}

		//Geant4 units are in mm
		RootEDepHist* ThisHist = new RootEDepHist(name, type, hgap, xmax, ymax, zmax, xstep, ystep, zstep, xbigmax, ybigmax, zbigmax, xbigstep, ybigstep, zbigstep);

		//Insert with a safety check
		std::pair<std::map< std::string,RootEDepHist* >::const_iterator,bool > check;
		check = HistogramMap.insert(std::pair<std::string,RootEDepHist*>(name, ThisHist));

		if(check.second == false)
		{
			std::cerr << "GEANT4> ERROR: Multiple definitions of collimator: \"" << name << "\"" << std::endl;
			exit(EXIT_FAILURE);
		}

		ThisEntry = ThisHist;
	}
}

void RootEnergyDeposition::Process(CollimationJawHitsCollection* hits)
{

	if(!hits)
	{
		std::cerr << "GEANT4> ERROR: Trying to process empty CollimationJawHitsCollection pointer!" << std::endl;
		exit(EXIT_FAILURE);
	}

	if(!ThisEntry)
	{
		//No processing to be done here
		return;
	}


	if(ThisEntry->GetHistogramType() == 2)
	{
		TH2F* xy = ThisEntry->GetHist_xy();
		TH2F* zx = ThisEntry->GetHist_zx();
		TH2F* zy = ThisEntry->GetHist_zy();
	
		if(!xy || !zx || !zy)
		{
			std::cerr << "GEANT4> ERROR: Trying to process empty histogram pointers!" << std::endl;
			exit(EXIT_FAILURE);
		}
	
		size_t entries = hits->entries();
		for(size_t n=0; n < entries; n++)
		{
			double edep = (*hits)[n]->GetEdep();

/*
			G4ThreeVector spos = (*hits)[n]->GetStartPosition();
			G4ThreeVector epos = (*hits)[n]->GetEndPosition();

			double sx = spos.x();
			double sy = spos.y();
			double sz = spos.z();

			double ex = epos.x();
			double ey = epos.y();
			double ez = epos.z();
*/

			double sx = (*hits)[n]->GetStartX();
			double sy = (*hits)[n]->GetStartY();
			double sz = (*hits)[n]->GetStartZ();

			double ex = (*hits)[n]->GetEndX();
			double ey = (*hits)[n]->GetEndY();
			double ez = (*hits)[n]->GetEndZ();

			FillHit(xy, sx, ex, sy, ey, edep);
			FillHit(zx, sz, ez, sx, ex, edep);
			FillHit(zy, sz, ez, sy, ey, edep);
		}
	}
	else if(ThisEntry->GetHistogramType() == 3)
	{
		TH3F* xyz = ThisEntry->GetHist3D();

		if(!xyz)
		{
			std::cerr << "GEANT4> ERROR: Trying to process empty histogram pointers!" << std::endl;
			exit(EXIT_FAILURE);
		}

		size_t entries = hits->entries();
		for(size_t n=0; n < entries; n++)
		{
			double edep = (*hits)[n]->GetEdep();
/*
			G4ThreeVector spos = (*hits)[n]->GetStartPosition();
			G4ThreeVector epos = (*hits)[n]->GetEndPosition();

			double sx = spos.x();
			double sy = spos.y();
			double sz = spos.z();

			double ex = epos.x();
			double ey = epos.y();
			double ez = epos.z();
*/

			double sx = (*hits)[n]->GetStartX();
			double sy = (*hits)[n]->GetStartY();
			double sz = (*hits)[n]->GetStartZ();

			double ex = (*hits)[n]->GetEndX();
			double ey = (*hits)[n]->GetEndY();
			double ez = (*hits)[n]->GetEndZ();
			FillHit(xyz, sx, ex, sy, ey, sz, ez, edep);
		}
	}
}

/*
 * If a forced flush is needed
 * */
void RootEnergyDeposition::WriteHistogramsToFile()
{

	std::map<std::string,RootEDepHist*>::iterator itr = HistogramMap.begin();
	while(itr != HistogramMap.end())
	{
		TH2F* xy = itr->second->GetHist_xy();
		TH2F* zx = itr->second->GetHist_zx();
		TH2F* zy = itr->second->GetHist_zy();
		xy->Write();
		zx->Write();
		zy->Write();
		itr++;
	}
}

void RootEnergyDeposition::FillHit(TH2F* hist, double xs, double xe, double ys, double ye, double edep)
{
	int n_slices = 50;

	//A very crude solution
	double xstep = (xe - xs) / n_slices;
	double ystep = (ye - ys) / n_slices;

	double estep = edep / n_slices;
	double x,y;

	for(int n=0; n < n_slices; n++)
	{
		x = xs + (n*xstep);
		y = ys + (n*ystep);
		hist->Fill(x,y, estep);
	}
}

void RootEnergyDeposition::FillHit(TH3F* hist, double xs, double xe, double ys, double ye, double zs, double ze, double edep)
{
	int n_slices = 100;

	//A very crude solution
	double xstep = (xe - xs) / n_slices;
	double ystep = (ye - ys) / n_slices;
	double zstep = (ze - zs) / n_slices;

	double estep = edep / n_slices;
	double x,y,z;

	for(int n=0; n < n_slices; n++)
	{
		x = xs + (n*xstep);
		y = ys + (n*ystep);
		z = zs + (n*zstep);
		hist->Fill(x,y,z,estep);
	}
}
void RootEnergyDeposition::WriteHist()
{
}
