#ifndef HASHTABLE_H_
#define HASHTABLE_H_

#include "EdbPattern.h"
#include <vector>
#include <EdbDataSet.h>

class HashTable {
	private:
		EdbPVRec* m_pvr;

		double m_vmin[3]; // minimum key value.
		double m_vmax[3]; // maximum key value.
		double m_div[3]; // 
		int m_n[3];

		std::vector<std::vector<std::vector<std::vector<EdbTrackP*>>>> m_table; // hash table.

	public:
		// constructor
		HashTable(EdbPVRec* pvr);

		std::vector<EdbTrackP*> GetNeighbors(EdbTrackP* track);
};

#endif
