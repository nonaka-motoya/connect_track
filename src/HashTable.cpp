#include "HashTable.hpp"
#include <algorithm>
#include <cmath>
#include <iterator>
#include <vector>
#include <numeric>

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

	std::cout << "First plate: " << m_tracks_f.front()->GetSegmentFirst()->ScanID().GetPlate() << "\tLast plate: " << m_tracks_f.back()->GetSegmentFirst()->ScanID().GetPlate() << std::endl;
}


std::vector<EdbTrackP*> HashTable::GetNeighbors(EdbTrackP* track) {
	int plate = track -> GetSegmentLast() -> ScanID().GetPlate();

	int dpl = 10;

	// forward tracks.
	std::vector<EdbTrackP*>::iterator iter_lower = std::lower_bound(m_tracks_f.begin(), m_tracks_f.end(), track,
			[&](const EdbTrackP* track1, EdbTrackP* track2) { return track1->GetSegmentFirst()->ScanID().GetPlate() < track2->GetSegmentLast()->ScanID().GetPlate() - dpl; }
			);

	std::vector<EdbTrackP*>::iterator iter_upper = std::lower_bound(m_tracks_f.begin(), m_tracks_f.end(), track,
			[&](const EdbTrackP* track1, EdbTrackP* track2) { return track1->GetSegmentFirst()->ScanID().GetPlate() < track2->GetSegmentLast()->ScanID().GetPlate() + dpl; }
			);

	int idx_lower = std::distance(m_tracks_f.begin(), iter_lower);
	int idx_upper = std::distance(m_tracks_f.begin(), iter_upper);

	std::vector<EdbTrackP*> tracks;

	for (int i=idx_lower; i<idx_upper; i++) {
		tracks.push_back(m_tracks_f[i]);
	}

	
	// backward tracks.
	iter_lower = std::lower_bound(m_tracks_l.begin(), m_tracks_l.end(), track,
			[&](const EdbTrackP* track1, EdbTrackP* track2) { return track1->GetSegmentLast()->ScanID().GetPlate() < track2->GetSegmentFirst()->ScanID().GetPlate() - dpl; }
			);

	iter_upper = std::lower_bound(m_tracks_l.begin(), m_tracks_l.end(), track,
			[&](const EdbTrackP* track1, EdbTrackP* track2) { return track1->GetSegmentLast()->ScanID().GetPlate() < track2->GetSegmentFirst()->ScanID().GetPlate() + dpl; }
			);

	idx_lower = std::distance(m_tracks_l.begin(), iter_lower);
	idx_upper = std::distance(m_tracks_l.begin(), iter_upper);

	for (int i=idx_lower; i<idx_upper; i++) {
		tracks.push_back(m_tracks_l[i]);
	}


	return tracks;
}
