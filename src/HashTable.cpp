#include "HashTable.hpp"
#include <cmath>
#include <vector>

HashTable::HashTable(EdbPVRec* pvr) {
	m_pvr = pvr;

	double xmax = -1e8;
	double ymax = -1e8;
	double pmax = -1e8;
	double xmin = 1e8;
	double ymin = 1e8;
	double pmin = 1e8;

	for (int i=0; i<m_pvr->Ntracks(); i++) {
		EdbTrackP *track = m_pvr -> GetTrack(i);
		for (int j=0; j<track->N(); j++) {
			EdbSegP *seg = track -> GetSegment(j);

			double x = seg -> X();
			double y = seg -> Y();
			int plate = seg -> ScanID().GetPlate();

			if (pmin > plate) pmin = plate;
			if (pmax < plate) pmax = plate;
			if (xmin > x) xmin = x;
			if (xmax < x) xmax = x;
			if (ymin > y) ymin = y;
			if (ymax < y) ymax = y;
		}
	}

	int cellsize = 20;
	int nx = floor((xmax - xmin)/cellsize);
	int ny = floor((ymax - ymin)/cellsize);
	int npl = pmax - pmin + 1;

	m_n[0] = nx;
	m_n[1] = ny;
	m_n[2] = npl;
	m_div[0] = (xmax - xmin)/nx;
	m_div[1] = (ymax - ymin)/ny;
	m_div[2] = (pmax - pmin)/npl;
	m_vmin[0] = xmin;
	m_vmin[1] = ymin;
	m_vmin[2] = pmin - 0.5;
	m_vmax[0] = xmax;
	m_vmax[1] = ymax;
	m_vmax[2] = pmax + 0.5;
}

std::vector<EdbTrackP*> HashTable::GetNeighbors(EdbTrackP *track) {
	std::vector<EdbTrackP*> v_tracks;
	double l_x = track -> GetSegmentLast() -> X();
	double l_y = track -> GetSegmentLast() -> Y();
	double f_x = track -> GetSegmentFirst() -> X();
	double f_y = track -> GetSegmentFirst() -> Y();
	int f_plate = track -> GetSegmentFirst() -> ScanID().GetPlate();
	int l_plate = track -> GetSegmentLast() -> ScanID().GetPlate();

	for (int i=0; i<m_pvr->Ntracks(); i++) {
		EdbTrackP *track_buf = m_pvr -> GetTrack(i);
		EdbSegP *f_seg = track_buf -> GetSegmentFirst();
		EdbSegP *l_seg = track_buf -> GetSegmentLast();

		// search with plate.
		int f_can_plate = f_seg -> ScanID().GetPlate();
		int l_can_plate = l_seg -> ScanID().GetPlate();
		if (abs(l_plate - f_can_plate) <= 6) {

			double f_can_x = f_seg -> X();
			double f_can_y = f_seg -> Y();

			//std::cout << "track ID: " << track -> GetSegmentFirst() -> Track() << "\tcand track ID: " << track_buf -> GetSegmentFirst() -> Track() << "\tdx: " << abs(x-f_x) << std::endl;
			
			
			// search with x, y
			if (abs(l_x - f_can_x) < 100 and abs(l_y - f_can_y) < 100) {
				//std::cout << "last track: " << f_seg -> Track() << std::endl;
				v_tracks.push_back(track_buf);
			}
		} else if (abs(f_plate - l_can_plate) <= 6) {

			double l_can_x = l_seg -> X();
			double l_can_y = l_seg -> Y();


			if (abs(f_x - l_can_x) < 100 and abs(f_y - l_can_y) < 100) {
				//std::cout << "front track: " << f_seg -> Track() << std::endl;
				v_tracks.push_back(track_buf);
			}
		}
	}

	return v_tracks;
}
