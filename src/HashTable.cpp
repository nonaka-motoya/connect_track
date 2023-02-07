#include "HashTable.hpp"
#include <algorithm>
#include <cmath>
#include <iterator>
#include <vector>
#include <numeric>

#include <EdbEDAUtil.h>

HashTable::HashTable(EdbPVRec* pvr) {
	m_pvr = pvr;

	for (int i=0; i<m_pvr->Ntracks(); i++) {
		EdbTrackP *track = m_pvr -> GetTrack(i);

		m_tracks_f.push_back(track);
		m_tracks_l.push_back(track);
	}

	std::cout << m_tracks_f.size() << " tracks are read." << std::endl;

	std::cout << "Sort tracks with plate." << std::endl;
	SortTracks();
}

bool CompareFirstPlate(const EdbTrackP* track1, const EdbTrackP* track2) {
	return track1->GetSegmentFirst()->ScanID().GetPlate() < track2->GetSegmentFirst()->ScanID().GetPlate();
}

bool CompareLastPlate(const EdbTrackP* track1, const EdbTrackP* track2) {
	return track1->GetSegmentLast()->ScanID().GetPlate() < track2->GetSegmentLast()->ScanID().GetPlate();
}

void HashTable::SortTracks() {
	std::sort(m_tracks_f.begin(), m_tracks_f.end(), CompareFirstPlate);
	std::sort(m_tracks_l.begin(), m_tracks_l.end(), CompareLastPlate);

	std::cout << "First plate: " << m_tracks_f.front()->GetSegmentFirst()->ScanID().GetPlate() << "\tLast plate: " << m_tracks_l.back()->GetSegmentLast()->ScanID().GetPlate() << std::endl;
}


std::vector<EdbTrackP*> HashTable::GetNeighbors(EdbTrackP* track, bool forward, bool backward) {
	std::vector<EdbTrackP*> tracks;
	double l_x = track -> GetSegmentLast() -> X();
	double l_y = track -> GetSegmentLast() -> Y();
	double f_x = track -> GetSegmentFirst() -> X();
	double f_y = track -> GetSegmentFirst() -> Y();

	int nseg = track -> N();
	EdbSegP* f_seg = track -> GetSegment(1);
	EdbSegP* l_seg = track -> GetSegment(nseg-1);

	//f_x -= 1000 * track -> GetSegmentFirst() -> TX();
	//f_y -= 1000 * track -> GetSegmentFirst() -> TY();
	//l_x += 1000 * track -> GetSegmentLast() -> TX();
	//l_y += 1000 * track -> GetSegmentLast() -> TY();

	int plate = track -> GetSegmentLast() -> ScanID().GetPlate();

	int dpl = 10;

	// forward tracks.
	if (forward) {
		std::vector<EdbTrackP*>::iterator iter_lower = std::lower_bound(m_tracks_f.begin(), m_tracks_f.end(), track,
				[&](const EdbTrackP* track1, EdbTrackP* track2) { return track1->GetSegmentFirst()->ScanID().GetPlate() < track2->GetSegmentLast()->ScanID().GetPlate() - dpl; }
				);

		std::vector<EdbTrackP*>::iterator iter_upper = std::lower_bound(m_tracks_f.begin(), m_tracks_f.end(), track,
				[&](const EdbTrackP* track1, EdbTrackP* track2) { return track1->GetSegmentFirst()->ScanID().GetPlate() < track2->GetSegmentLast()->ScanID().GetPlate() + dpl; }
				);

		int idx_lower = std::distance(m_tracks_f.begin(), iter_lower);
		int idx_upper = std::distance(m_tracks_f.begin(), iter_upper);


		for (int i=idx_lower; i<idx_upper; i++) {
			EdbSegP* cand_f_seg = m_tracks_f[i] -> GetSegment(1);
			double dmin = EdbEDAUtil::CalcDmin(l_seg, cand_f_seg);
			if (dmin > 50) continue;
			if (track -> GetSegmentFirst() -> ID() == m_tracks_f[i] -> GetSegmentFirst() -> ID()) continue;
			tracks.push_back(m_tracks_f[i]);
		}
	}

	
	// backward tracks.
	if (backward) {
		std::vector<EdbTrackP*>::iterator iter_lower = std::lower_bound(m_tracks_l.begin(), m_tracks_l.end(), track,
				[&](const EdbTrackP* track1, EdbTrackP* track2) { return track1->GetSegmentLast()->ScanID().GetPlate() < track2->GetSegmentFirst()->ScanID().GetPlate() - dpl; }
				);

		std::vector<EdbTrackP*>::iterator iter_upper = std::lower_bound(m_tracks_l.begin(), m_tracks_l.end(), track,
				[&](const EdbTrackP* track1, EdbTrackP* track2) { return track1->GetSegmentLast()->ScanID().GetPlate() < track2->GetSegmentFirst()->ScanID().GetPlate() + dpl; }
				);

		int idx_lower = std::distance(m_tracks_l.begin(), iter_lower);
		int idx_upper = std::distance(m_tracks_l.begin(), iter_upper);

		for (int i=idx_lower; i<idx_upper; i++) {
			int cand_nseg = m_tracks_l[i] -> N();
			EdbSegP* cand_l_seg = m_tracks_l[i] -> GetSegment(cand_nseg-1);
			double dmin = EdbEDAUtil::CalcDmin(f_seg, cand_l_seg);
			if (dmin > 50) continue;
			if (track -> GetSegmentFirst() -> ID() == m_tracks_l[i] -> GetSegmentFirst() -> ID()) continue;
			tracks.push_back(m_tracks_l[i]);
		}
	}

	return tracks;
}
