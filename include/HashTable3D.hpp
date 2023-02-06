#ifndef HASHTABLE3D_H_
#define HASHTABLE3D_H_

#include "EdbPattern.h"
#include <vector>
#include <EdbDataSet.h>

class HashTable3D {
	private:
		EdbPVRec* m_pvr;

		std::vector<EdbTrackP*> m_neighbors;

		double m_vmin[3];
		double m_vmax[3];
		double m_div[3];
		int m_n[3];

		std::vector<std::vector<std::vector<std::vector<EdbTrackP*>>>> m_table_f;
		std::vector<std::vector<std::vector<std::vector<EdbTrackP*>>>> m_table_b;

	public:
		// constructor
		HashTable3D(EdbPVRec* pvr, int n1, double min1, double max1, int n2, double min2, double max2, int n3, double min3, double max3);

		int getCellIndex(int iv, double v) {
		  int idx = (v - m_vmin[iv])/m_div[iv]; // calculate index
		  if(idx < 0) idx = 0;
		  if(idx >= m_n[iv]) idx = m_n[iv] - 1;
		  return idx;
		}

		std::vector<EdbTrackP*>& findCell(double* cell) {
		  Int_t i = getCellIndex(0, cell[0]);
		  Int_t j = getCellIndex(1, cell[1]);
		  Int_t k = getCellIndex(2, cell[2]);
		  return m_table_f[i][j][k];
		}

		int fillCells(std::vector<EdbTrackP*> vec);
                                   
		// get around tracks.      
		std::vector<EdbTrackP*>& GetNeighbors(double* v, double* r);

		void SortTracks();

};

#endif
