#include "BunchDump_root.h"

#include "TTree.h"
#include "TROOT.h"

#include <iostream>
#include <string>
#include <algorithm>
#include <cctype>
#include <cstring>

BunchDumpRootOutput* RootBunchDumpOutput = nullptr;

extern "C" void root_BunchDumpInit()
{
	if(!RootBunchDumpOutput)
	RootBunchDumpOutput = new BunchDumpRootOutput();
}

extern "C" void root_DumpBunch(char* name_in, int name_len, int i_in, int ix_in, int turn_in, int particleID_in, int parentID_in, int pdgID_in, int16_t q_in, double weight_in, double s_in, double x_in, double xp_in, double y_in, double yp_in, double z_in, double dp_in, double sx_in, double sy_in, double sz_in, double mass_in)
{
	if(RootBunchDumpOutput)
	{
		RootBunchDumpOutput->DumpBunch(name_in, name_len, i_in, ix_in, turn_in, particleID_in, parentID_in, pdgID_in, q_in, weight_in, s_in, x_in, xp_in, y_in, yp_in, z_in, dp_in, sx_in, sy_in, sz_in, mass_in);
	}
	else
	{
		std::cerr << "RootBunchDumpOutput is not yet created" << std::endl;
		exit(EXIT_FAILURE);
	}
}


//Class functions
BunchDumpRootOutput::BunchDumpRootOutput()
{
	//Tree stuff
	BunchDumpTree = new TTree("BunchDump","BunchDumpTree");
	BunchDumpTree->Branch("name",		name,		"name[49]/C");
	BunchDumpTree->Branch("i",			&i,			"i/I");
	BunchDumpTree->Branch("ix",			&ix,		"ix/I");
	BunchDumpTree->Branch("turn",		&turn,		"turn/I");
	BunchDumpTree->Branch("particleID",	&particleID,"particleID/I");
	BunchDumpTree->Branch("parentID",	&parentID,	"parentID/I");
	BunchDumpTree->Branch("pdgID",		&pdgID,		"pdgID/I");
	BunchDumpTree->Branch("q",			&q,			"q/S");
	BunchDumpTree->Branch("weight",		&weight,	"weight/D");
	BunchDumpTree->Branch("s",			&s,			"s/D");
	BunchDumpTree->Branch("x",			&x,			"x/D");
	BunchDumpTree->Branch("xp",			&xp,		"xp/D");
	BunchDumpTree->Branch("y",			&y,			"y/D");
	BunchDumpTree->Branch("yp",			&yp,		"yp/D");
	BunchDumpTree->Branch("z",			&z,			"z/D");
	BunchDumpTree->Branch("dp",			&dp,		"dp/D");
	BunchDumpTree->Branch("sx",			&sx,		"sx/D");
	BunchDumpTree->Branch("sy",			&sy,		"sy/D");
	BunchDumpTree->Branch("sz",			&sz,		"sz/D");
	BunchDumpTree->Branch("m",			&mass,		"m/D");
}

void BunchDumpRootOutput::DumpBunch(char* name_in, int name_len, int i_in, int ix_in, int turn_in, int particleID_in, int parentID_in, int pdgID_in, int16_t q_in, double weight_in, double s_in, double x_in, double xp_in, double y_in, double yp_in, double z_in, double dp_in, double sx_in, double sy_in, double sz_in, double mass_in)
{
	std::string tmpname(name_in);
	std::transform(tmpname.begin(), tmpname.end(), tmpname.begin(), ::toupper);
	tmpname = tmpname.substr(0,name_len);
	strncpy(name,tmpname.c_str(), sizeof(name));

	i = i_in;
	ix = ix_in;
	turn = turn_in;
	particleID = particleID_in;
	parentID = parentID_in;
	pdgID = pdgID_in;
	q = q_in;
	weight = weight_in;
	s  = s_in;
	x  = x_in;
	xp = xp_in;
	y  = y_in;
	yp = yp_in;
	z  = z_in;
	dp = dp_in;
	sx = sx_in;
	sy = sy_in;
	sz = sz_in;
	mass  = mass_in;

	BunchDumpTree->Fill();
}

