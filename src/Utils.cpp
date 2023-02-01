#include "Utils.hpp"

#include "EdbPattern.h"
#include "EdbSegP.h"
#include "TGraph.h"

#include <EdbDataSet.h>
#include <EdbEDAUtil.h>
#include <EdbVertex.h>

#include "HashTable.hpp"


bool Utils::IsSwitch(EdbTrackP* track, int index) {
	EdbSegP* prv_seg = track -> GetSegment(index-1);
	EdbSegP* seg = track -> GetSegment(index);
	EdbSegP* nxt_seg = track -> GetSegment(index+1);

	// if previous MC track ID is different from current one,
	// and if current MC track ID is equal to next one,
	if (prv_seg -> Volume() != seg -> Volume() and seg -> Volume() == nxt_seg -> Volume()) {
		return true;
	}

	return false;
}

int Utils::CheckTrackingQuality(EdbTrackP* track) {

	EdbSegP* first_seg = track -> GetSegmentFirst();

	int num_diff = 0;

	// start from 5th segment.
	for (int i=1; i<track->N()-1; i++) {
		if (Utils::IsSwitch(track, i)) {
			num_diff ++;
		}
	}

	return num_diff;
}

void Utils::PrintTrack(EdbTrackP* track) {
	
	std::cout << "MC event ID: " << track -> GetSegmentFirst() -> MCEvt();
	std::cout << "\tMC track ID: " << track -> GetSegmentFirst() -> Volume();
	std::cout << "\tPDG ID: " << track -> GetSegmentFirst() -> MCTrack();
	std::cout << "\tP: " << track -> GetSegmentFirst() -> P() << "GeV";
	std::cout << "\tplate: " << track -> GetSegmentFirst() -> ScanID().GetPlate();
	std::cout << "\tnpl: " << track -> Npl();
	std::cout << "\tquality: " << Utils::CheckTrackingQuality(track);
	std::cout << std::endl;

	return;
}

int Utils::SearchVertex(EdbTrackP* track, int index, EdbPVRec* pvr) {

	EdbSegP* seg_1ry = track -> GetSegment(index);
	int plate_1ry = seg_1ry -> ScanID().GetPlate();

	EdbVertex* v = new EdbVertex;
	for (int j=0; j<pvr->Ntracks(); j++) {
		EdbTrackP* track_buf = pvr -> GetTrack(j);

		// if track start from next plate of the segment.
		int plate_2ry = track_buf -> GetSegmentFirst() -> ScanID().GetPlate();
		if (plate_2ry <= plate_1ry + 2 and plate_2ry > plate_1ry) {
			
			double dmin = EdbEDAUtil::CalcDmin(seg_1ry, track_buf->GetSegmentFirst());

			if (dmin < 1) {
				EdbVTA* vta = new EdbVTA(track, v);
				v -> AddVTA(vta);
			}
		}
	}
	return v -> Nn();

}


std::pair<double, double> Utils::CalcTrackAngle(EdbTrackP* track, int index) {

	TGraph* grx = new TGraph();
	TGraph* gry = new TGraph();

	for (int i=index-1; i<=index+1; i++) {
		EdbSegP* seg = track -> GetSegment(i);

		grx -> SetPoint(i-index+1, seg->Z(), seg->X());
		gry -> SetPoint(i-index+1, seg->Z(), seg->Y());
	}

	grx -> Fit("pol1", "Q");
	gry -> Fit("pol1", "Q");

	double tx = grx -> GetFunction("pol1") -> GetParameter(1);
	double ty = gry -> GetFunction("pol1") -> GetParameter(1);

	return {tx, ty};
}

double Utils::CalcTrackAngleDiff(EdbTrackP* track, int index) {
	
	std::pair<double, double> prv_theta = Utils::CalcTrackAngle(track, index-2);
	std::pair<double, double> nxt_theta = Utils::CalcTrackAngle(track, index+1);

	double thx1 = prv_theta.first;
	double thy1 = prv_theta.second;
	double thx2 = nxt_theta.first;
	double thy2 = nxt_theta.second;

	double theta = sqrt((thx1-thx2)*(thx1-thx2) + (thy1-thy2)*(thy1-thy2));
	return theta * 1000;
}


// calculate track angle using all segments.
std::pair<double, double> theta(EdbTrackP *track) {

	TGraph* grx = new TGraph();
	TGraph* gry = new TGraph();

	for (int i=0; i<track->N(); i++) {
		EdbSegP* seg = track -> GetSegment(i);

		grx -> SetPoint(i, seg->Z(), seg->X());
		gry -> SetPoint(i, seg->Z(), seg->Y());
	}

	grx -> Fit("pol1", "Q");
	gry -> Fit("pol1", "Q");

	double tx = grx -> GetFunction("pol1") -> GetParameter(1);
	double ty = gry -> GetFunction("pol1") -> GetParameter(1);

	return {tx, ty};
}

// calculate track angle difference between two tracks.
double dtheta(EdbTrackP *track1, EdbTrackP *track2) {

	std::pair<double, double> theta1 = theta(track1);
	std::pair<double, double> theta2 = theta(track2);

	double dtheta = sqrt((theta2.first-theta1.first)*(theta2.first-theta1.first) + (theta2.second-theta1.second)*(theta2.second-theta1.second));

	return dtheta;
}

bool isDuplicate(EdbTrackP* track1, EdbTrackP* track2) {

	for (int i=0; i<track1->N(); i++) {
		int plate1 = track1 -> GetSegment(i) -> ScanID().GetPlate();

		for (int j=0; j<track2->N(); j++) {
			int plate2 = track2 -> GetSegment(j) -> ScanID().GetPlate();

			if (plate1 == plate2) return true;
		}
	}

	return false;
}

EdbSegP *make_virtual_segment(EdbTrackP *track) {

	EdbSegP *seg = new EdbSegP;
	
	TGraph* grx = new TGraph();
	TGraph* gry = new TGraph();

	for (int i=0; i<track->N(); i++) {
		EdbSegP *seg = track -> GetSegment(i);

		grx -> SetPoint(i, seg->Z(), seg->X());
		gry -> SetPoint(i, seg->Z(), seg->Y());
	}

	grx -> Fit("pol1", "Q");
	gry -> Fit("pol1", "Q");

	double x = track -> GetSegment(track->N()/2) -> X();
	double y = track -> GetSegment(track->N()/2) -> Y();
	double z = track -> GetSegment(track->N()/2) -> Z();

	double tx = grx -> GetFunction("pol1") -> GetParameter(1);
	double ty = gry -> GetFunction("pol1") -> GetParameter(1);

	int segid = track -> GetSegmentFirst() -> ID();
	int pid = track -> GetSegmentFirst() -> PID();

	seg -> Set(segid, x, y, tx, ty, 23, 0);
	seg -> SetZ(z);
	seg -> SetPID(pid);

	return seg;
}

// calculate distance between two tracks.
double distance(EdbTrackP *track1, EdbTrackP *track2) {
	//EdbSegP *seg1 = track1 -> GetSegment(track1->N()/2);
	//EdbSegP *seg2 = track2 -> GetSegment(track2->N()/2);

	EdbSegP *seg1 = make_virtual_segment(track1);
	EdbSegP *seg2 = make_virtual_segment(track2);

	double distance = EdbEDAUtil::CalcDmin(seg1, seg2);
	return distance;
}

// there are bugs.
int get_min_chi_index(EdbTrackP *track, std::vector<EdbTrackP*> v_tracks) {

	std::cout << "There are " << v_tracks.size() << " candidates." << std::endl;

	double dist = 1e8;
	int index = 0;

	for (int i=0; i<v_tracks.size(); i++) {

		EdbTrackP *cand_track = v_tracks[i];

		double dist_buf = distance(track, cand_track);


		if (dist_buf < dist) {
			dist = dist_buf;
			index = i;
		}
	}

	return index;
}

void add_tracks(EdbTrackP *track1, EdbTrackP *track2) {
	for (int i=0; i<track2->N(); i++) {
		EdbSegP *seg = track2 -> GetSegment(i);

		track1 -> AddSegment(seg);
	}

	for (int i=0; i<track2->NF(); i++) {
		EdbSegP *segf = track2 -> GetSegmentF(i);

		track1 -> AddSegmentF(segf);
	}

	track1 -> SetNpl();
	track1 -> SetSegmentsTrack();
	track1 -> SetCounters();

	return;
}


EdbPVRec *Utils::ConnectTrack(std::string path) {

	EdbPVRec *connected_pvr = new EdbPVRec;
	

	// read linked_tracks.root.
	EdbDataProc* dproc = new EdbDataProc;
	EdbPVRec* pvr = new EdbPVRec;
	dproc -> ReadTracksTree(*pvr, path.c_str(), "nseg>3");

	// create HashTable
	std::cout << "Fill HashTable." << std::endl;
	HashTable *hashtable = new HashTable(pvr); // hashtable of EdbTrackP.
	
	std::vector<int> removed_track_id; // track ID which should be removed.

	TObjArray *selected = new TObjArray; // array of EdbTrackP which is connected.

	for (int i=0; i<pvr->Ntracks(); i++) {
		EdbTrackP *track = pvr -> GetTrack(i);

		if (abs(track -> GetSegmentFirst() -> MCTrack()) != 13) continue;

		// if track should be removed -> skip.
		int track_id = track-> GetSegmentFirst() -> Track();
		if (std::find(removed_track_id.begin(), removed_track_id.end(), track_id) != removed_track_id.end()) continue;
		std::cout << "target track ID: " << track_id << std::endl;

		while (1) {
			// get neipghbor tracks.
			std::vector<EdbTrackP*> v_tracks = hashtable -> GetNeighbors(track);
			std::cout << "Number of neighbor tracks: " << v_tracks.size() << std::endl;

			// if neighbor track does not exist -> break.
			if (v_tracks.size() == 0) break;

			std::vector<EdbTrackP*> v_tbc_tracks; // tracks to be connected.
			for (int j=0; j<v_tracks.size(); j++) {
				EdbTrackP *cand_track = v_tracks[j];
				int cand_track_id = cand_track -> GetSegmentFirst() -> Track();

				// if candidate track ID is equal to the track -> skip.
				if (track_id == cand_track_id) continue;

				// if candidate track is to be removed -> skip.
				if (std::find(removed_track_id.begin(), removed_track_id.end(), cand_track_id) != removed_track_id.end()) continue;


				// if angle is smaller than threshold angle
				// and if distance between two tracks is smaller than threshold distance
				// and if two tracks have no common segment
				double delta_theta = dtheta(track, cand_track);
				double dist = distance(track, cand_track);
				std::cout << "MC track ID: " << track -> GetSegmentFirst() -> MCTrack() << "\ttrack ID: " << track -> GetSegmentFirst() -> Track() << "\tcandidate track ID: " << cand_track -> GetSegmentFirst() -> Track() << "\tdtheta: " << delta_theta << "\tdistance: " << dist << std::endl;

				if (!isDuplicate(track, cand_track) and delta_theta < 1 and dist < 200) {

					// fill vector of tracks which are to be connected.
					// after connecting this track, need to remove this track from pvr.
					std::cout << "connect!" << std::endl;
					v_tbc_tracks.push_back(cand_track);
				}
			}

			if (v_tbc_tracks.size() == 0) {
				break;
			} else if (v_tbc_tracks.size() == 1) {
				add_tracks(track, v_tbc_tracks.front());
				removed_track_id.push_back(v_tbc_tracks.front() -> GetSegmentFirst() -> Track());
				//track -> PrintNice();
			} else {
				int index = get_min_chi_index(track, v_tbc_tracks);
				add_tracks(track, v_tbc_tracks[index]);
				removed_track_id.push_back(v_tbc_tracks[index] -> GetSegmentFirst() -> Track());
			}
		}

		selected -> Add(track);
		connected_pvr  -> AddTrack(track);
		removed_track_id.push_back(track->GetSegmentFirst()->Track());
	}

	dproc -> MakeTracksTree(connected_pvr, "connect_tracks.root");

	// check pvr
	//pvr -> PrintSummary();
	//for (int i=0; i<pvr->Ntracks(); i++) {
	//	EdbTrackP *track = pvr -> GetTrack(i);
	//	if (track -> GetSegmentFirst() -> MCTrack() == 13) {
	//		track -> PrintNice();
	//	}
	//}

	return connected_pvr;
}
