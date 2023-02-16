#include <iostream>
#include <vector>
#include <fstream>
#include <ostream>

#include <EdbDataSet.h>

#include "Utils.hpp"

bool isMuonOutgo(std::string path) {

	int event_id = 0;
	sscanf(path.c_str(), "20230107_nuall/evt_%d", &event_id);
	event_id += 100000;

	EdbDataProc* dproc = new EdbDataProc;
	EdbPVRec* pvr = new EdbPVRec;

	dproc -> ReadTracksTree(*pvr, path.c_str(), "t.eP>20");

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
	input_file.open("input_list_2.txt", std::ios::in);
	std::string line_buf;
	int num_file = 0;

	std::ofstream ofs("outgo.txt");


	// read linked_tracks.root
	while (std::getline(input_file, line_buf)) {
		std::string path = line_buf + "/linked_tracks.root"; // path of linked_tracks.
		bool is_outgo = isMuonOutgo(path);
		if (is_outgo) {
			ofs << path << std::endl;
		}
		num_file ++;
	}

	isMuonOutgo("20230107_nuall/evt_2950_pl1_300/linked_tracks.root");
	isMuonOutgo("20230107_nuall/evt_1015_pl1_300/linked_tracks.root");

}
