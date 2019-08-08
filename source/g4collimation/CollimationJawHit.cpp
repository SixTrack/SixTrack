#include "CollimationJawHit.h"

//CollimationJawHit::CollimationJawHit() : G4VHit(), TrackID(-1), edep(0.0), StartPosition(G4ThreeVector()), EndPosition(G4ThreeVector())
CollimationJawHit::CollimationJawHit() : G4VHit(),  edep(0.0), sx(0), sy(0), sz(0), ex(0), ey(0), ez(0) 
{

}

CollimationJawHit::~CollimationJawHit()
{

}


void CollimationJawHit::Draw()
{

}

void CollimationJawHit::Print()
{

}
/*
void CollimationJawHit::SetTrackID(int id_in)
{
	TrackID = id_in;
}

void CollimationJawHit::SetPDGid(int id_in)
{
	PDGid = id_in;
}
*/
void CollimationJawHit::SetEdep(double edep_in)
{
	edep = edep_in;
}

void CollimationJawHit::SetStartPosition(G4ThreeVector pos_in)
{
//	StartPosition = pos_in;
	sx = pos_in.x();
	sy = pos_in.y();
	sz = pos_in.z();
}

void CollimationJawHit::SetEndPosition(G4ThreeVector pos_in)
{
//	EndPosition = pos_in;
	ex = pos_in.x();
	ey = pos_in.y();
	ez = pos_in.z();
}
/*
int CollimationJawHit::GetTrackID() const
{
	return TrackID;
}

int CollimationJawHit::GetPDGid() const
{
	return PDGid;
}
*/
double CollimationJawHit::GetEdep() const
{
	return edep;
}
/*
G4ThreeVector CollimationJawHit::GetStartPosition() const
{
	return StartPosition;
}

G4ThreeVector CollimationJawHit::GetEndPosition() const
{
	return EndPosition;
}
*/

