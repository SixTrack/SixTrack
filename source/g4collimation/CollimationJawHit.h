#ifndef _CollimationJawHit_h
#define _CollimationJawHit_h

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"

class CollimationJawHit : public G4VHit
{
public:
	CollimationJawHit();
	~CollimationJawHit();

	void Draw();
	void Print();

//    void SetTrackID(int);
//    void SetPDGid(int);
    void SetEdep(double);
    void SetStartPosition(G4ThreeVector);
    void SetEndPosition(G4ThreeVector);


//    int GetTrackID() const;
//    int GetPDGid() const;
    double GetEdep() const;

//    G4ThreeVector GetStartPosition() const;
//    G4ThreeVector GetEndPosition() const;

	inline float GetStartX() const
	{
		return sx;
	}

	inline float GetStartY() const
	{
		return sy;
	}

	inline float GetStartZ() const
	{
		return sz;
	}

	inline float GetEndX() const
	{
		return ex;
	}

	inline float GetEndY() const
	{
		return ey;
	}

	inline float GetEndZ() const
	{
		return ez;
	}

private:

//	Could be nice if there is more storage
//	int TrackID;
//	int PDGid;

	double edep;

//  Why not use a G4ThreeVector? It stores everything as a double, and for high energy events memory usage rapidly becomes an issue.
//	G4ThreeVector StartPosition;
//	G4ThreeVector EndPosition;

	//Start positions
	float sx;
	float sy;
	float sz;

	//End positions
	float ex;
	float ey;
	float ez;
};

typedef G4THitsCollection<CollimationJawHit> CollimationJawHitsCollection;

#endif

