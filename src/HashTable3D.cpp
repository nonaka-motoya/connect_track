#include "HashTable3D.hpp"
#include "EdbPattern.h"

#include <iostream>
#include <vector>

#include <EdbDataSet.h>


HashTable3D::HashTable3D(EdbPVRec* pvr, int n1, double min1, double max1, int n2, double min2, double max2, int n3, double min3, double max3) : m_vmin{min1, min2, min3}, m_vmax{max1, max2, max3}, m_n{n1, n2, n3} {
	printf("Constructor: setting cells\n");

	m_div[0] = (max1 - min1)/n1;
	m_div[1] = (max2 - min2)/n2;
	m_div[2] = (max3 - min3)/n3;

	m_table_f.resize(n1);
	for (Int_t i = 0; i < n1; i++) {
		m_table_f[i].resize(n2);
		for (Int_t j = 0; j < n2; j++) {
			m_table_f[i][j].resize(n3);
		}
	}
}

int HashTable3D::fillCells(std::vector<EdbTrackP*> vec) {
	printf("Filling cells\n");
	Int_t vec_size = vec.size();
	for (Int_t i = 0; i < vec_size; i++) {
		double data[3];
		data[0] = vec[i] -> GetSegmentFirst() -> X();
		data[1] = vec[i] -> GetSegmentFirst() -> Y();
		data[2] = vec[i] -> GetSegmentFirst() -> ScanID().GetPlate();

		std::vector<EdbTrackP*>& cell = findCell(data);
		cell.push_back(vec[i]);
	}
	return 0;
}

std::vector<EdbTrackP*>& HashTable3D::GetNeighbors(double* v, double* r) {
	// minimum and maximum index.
	int ivmin[3], ivmax[3];
	for (int i = 0; i < 3; i++) {
		ivmin[i] = getCellIndex(i, v[i] - r[i]);
		ivmax[i] = getCellIndex(i, v[i] + r[i]);
	}


	// max and min value.
	double vnmin[3] = {v[0] - r[0], v[1] - r[1], v[2] - r[2]};
	double vnmax[3] = {v[0] + r[0], v[1] + r[1], v[2] + r[2]};

	m_neighbors.clear();
	for (Int_t i = ivmin[0]; i <= ivmax[0]; i++) {
		for (Int_t j = ivmin[1]; j <= ivmax[1]; j++) {
			for (Int_t k = ivmin[2]; k <= ivmax[2]; k++) {
				std::vector<EdbTrackP*>& vec = m_table_f[i][j][k];
				Int_t vec_size = vec.size();
				for (Int_t l = 0; l < vec_size; l++) {
					//if (vnmin[0] < vec[l]->getValue(0) and vec[l]->getValue(0) < vnmax[0]) {
					//	if (vnmin[1] < vec[l]->getValue(1) and vec[l]->getValue(1) < vnmax[1]) {
					//	  if (vnmin[2] < vec[l]->getValue(2) and vec[l]->getValue(2) < vnmax[2]) {
					//		m_neighbors.push_back(vec[l]);
					//	  }
					//	}
					//}
				}
			}
		}
	}
  return m_neighbors;
}
