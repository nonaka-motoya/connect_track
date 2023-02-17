#include <iostream>
#include <vector>
#include <fstream>
#include <ostream>

#include <EdbDataSet.h>

#include "Utils.hpp"

bool isMuonOutgo(std::string path, EdbPVRec* pvr) {

	int event_id = 0;
	sscanf(path.c_str(), "20230203_nuall_enlarged/evt_%d", &event_id);
	event_id += 100000;


	for (int i=0; i<pvr->Ntracks(); i++) {
		EdbTrackP* track = pvr -> GetTrack(i);
		
		if (abs(track -> GetSegmentFirst() -> MCTrack()) == 13 and track -> GetSegmentFirst() -> MCEvt() == event_id ) {
			if (Utils::IsOutgo(track, pvr)) {
				return true;
			}
			break;
		}
	}

	return false;
}


int main() {
	std::ifstream input_file;
	input_file.open("./input_list_enlarged.txt", std::ios::in);
	std::string line_buf;
	int num_file = 0;

	std::ofstream ofs("outgo.txt");

	EdbDataProc* dproc = new EdbDataProc;
	EdbPVRec* pvr = new EdbPVRec;


	// read linked_tracks.root
	while (std::getline(input_file, line_buf)) {
		std::string path = line_buf + "/linked_tracks.root"; // path of linked_tracks.
		dproc -> ReadTracksTree(*pvr, path.c_str(), "t.eP>20");

		bool is_outgo = isMuonOutgo(path, pvr);
		if (is_outgo) {
			ofs << path << std::endl;
		}

		if (pvr -> eTracks) pvr -> eTracks -> Clear();
		num_file ++;
	}
}
