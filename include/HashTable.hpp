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

		// get around tracks.
		std::vector<EdbTrackP*> GetNeighbors(EdbTrackP* track);

		// get index.
		int GetCellIndex(int iv, double v) {
			int idx = (v - m_vmin[iv]) / m_div[iv];
			if (idx >= m_n[iv]) idx = m_n[iv] - 1;

			return idx;
		}

		// find cell.
		std::vector<EdbTrackP*> FindCell(double cell[3]) {
			int i = GetCellIndex(0, cell[0]);
			int j = GetCellIndex(1, cell[1]);
			int k = GetCellIndex(2, cell[2]);

			return m_table[i][j][k];
		}
};

#endif
