#ifndef HASHTABLE_H_
#define HASHTABLE_H_

#include "EdbPattern.h"
#include <vector>
#include <EdbDataSet.h>

class HashTable {
	private:
		EdbPVRec* m_pvr;

		std::vector<EdbTrackP*> m_tracks_f; // sorted with first segment plate.
		std::vector<EdbTrackP*> m_tracks_l; // sorted with last segment plate.


	public:
		// constructor
		HashTable(EdbPVRec* pvr);

		// get around tracks.
		std::vector<EdbTrackP*> GetNeighbors(EdbTrackP* track, bool forward=true, bool backward=true);

		void SortTracks();

};

#endif
