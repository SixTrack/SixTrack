#include "ApertureDump_root.h"

#include "TTree.h"
#include "TROOT.h"

#include <iostream>
#include <string>
#include <algorithm>
#include <cctype>
#include <cstring>

ApertureDumpRootOutput* RootApertureDumpOutput;

extern "C" void root_ApertureDumpInit()
{
	RootApertureDumpOutput = new ApertureDumpRootOutput();
}

extern "C" void root_DumpAperture(char* apname_in, int apname_len, char* aptype_in, int aptype_len, double s_in, double ap1_in, double ap2_in, double ap3_in, double ap4_in, double ap5_in, double ap6_in, double ap7_in, double ap8_in, double ap9_in, double ap10_in, double ap11_in )
{
    if(RootApertureDumpOutput)
    {
	    RootApertureDumpOutput->DumpAperture(apname_in, apname_len, aptype_in, aptype_len, s_in, ap1_in, ap2_in, ap3_in, ap4_in, ap5_in, ap6_in, ap7_in, ap8_in, ap9_in, ap10_in, ap11_in );
    }
    else
    {
        std::cerr << "RootApertureDumpOutput is not yet created" << std::endl;
        exit(EXIT_FAILURE);
    }
}


//Class functions
ApertureDumpRootOutput::ApertureDumpRootOutput()
{
	//Tree stuff
	ApertureDumpTree = new TTree("ApertureDump","ApertureDumpTree");
	ApertureDumpTree->Branch("name",name,"name[49]/C");
	ApertureDumpTree->Branch("type",type,"type[3]/C");
	ApertureDumpTree->Branch("s",&s,"s/D");
	ApertureDumpTree->Branch("ap1", &ap1, "ap1/D");
	ApertureDumpTree->Branch("ap2", &ap2, "ap2/D");
	ApertureDumpTree->Branch("ap3", &ap3, "ap3/D");
	ApertureDumpTree->Branch("ap4", &ap4, "ap4/D");
	ApertureDumpTree->Branch("ap5", &ap5, "ap5/D");
	ApertureDumpTree->Branch("ap6", &ap6, "ap6/D");
	ApertureDumpTree->Branch("ap7", &ap7, "ap7/D");
	ApertureDumpTree->Branch("ap8", &ap8, "ap8/D");
	ApertureDumpTree->Branch("ap9", &ap9, "ap9/D");
	ApertureDumpTree->Branch("ap10", &ap10, "ap10/D");
	ApertureDumpTree->Branch("ap11", &ap11, "ap11/D");
}

void ApertureDumpRootOutput::DumpAperture(char* apname_in, int apname_len, char* aptype_in, int aptype_len, double s_in, double ap1_in, double ap2_in, double ap3_in, double ap4_in, double ap5_in, double ap6_in, double ap7_in, double ap8_in, double ap9_in, double ap10_in, double ap11_in )
{
	std::string tmpname(apname_in);
	std::transform(tmpname.begin(), tmpname.end(), tmpname.begin(), ::toupper);
    tmpname = tmpname.substr(0,apname_len);
	strncpy(name,tmpname.c_str(), sizeof(name));

	std::string tmptype(aptype_in);
	std::transform(tmptype.begin(), tmptype.end(), tmptype.begin(), ::toupper);
    tmptype = tmptype.substr(0,aptype_len);
	strncpy(type, tmptype.c_str(), sizeof(type));

	s = s_in;
	ap1 = ap1_in;
	ap2 = ap2_in;
	ap3 = ap3_in;
	ap4 = ap4_in;
	ap5 = ap5_in;
	ap6 = ap6_in;
	ap7 = ap7_in;
	ap8 = ap8_in;
	ap9 = ap9_in;
	ap10 = ap10_in;
	ap11 = ap11_in;

	ApertureDumpTree->Fill();
}

