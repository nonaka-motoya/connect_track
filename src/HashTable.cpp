#include "HashTable.hpp"
#include <vector>

HashTable::HashTable(EdbPVRec* pvr) {
	m_pvr = pvr;
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
		if (abs(l_plate - f_can_plate) <= 20) {

			double f_can_x = f_seg -> X();
			double f_can_y = f_seg -> Y();

			//std::cout << "track ID: " << track -> GetSegmentFirst() -> Track() << "\tcand track ID: " << track_buf -> GetSegmentFirst() -> Track() << "\tdx: " << abs(x-f_x) << std::endl;
			
			
			// search with x, y
			if (abs(l_x - f_can_x) < 1000 and abs(l_y - f_can_y) < 1000) {
				//std::cout << "last track: " << f_seg -> Track() << std::endl;
				v_tracks.push_back(track_buf);
			}
		} else if (abs(f_plate - l_can_plate) <= 20) {

			double l_can_x = l_seg -> X();
			double l_can_y = l_seg -> Y();


			if (abs(f_x - l_can_x) < 1000 and abs(f_y - l_can_y) < 1000) {
				//std::cout << "front track: " << f_seg -> Track() << std::endl;
				v_tracks.push_back(track_buf);
			}
		}
	}

	return v_tracks;
}
