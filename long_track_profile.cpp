#include <iostream>
#include <fstream>
#include <sstream>

#include "TH1.h"
#include "TRint.h"

#include "EdbDataSet.h"

#include "Utils.hpp"

int bin_max = 300;
TCut c = "npl>=200&&t.eP>0"; // cut values in reading linked_tracks

std::map<int, std::vector<int>> uniqueID;
int ntrack;

void initTruth(std::string filename) {

	std::cout << "Read " << filename << std::endl;

	TFile* file = new TFile(filename.c_str(), "READ");
	TTree* tree = (TTree*) file -> Get("m_NuMCTruth_tree");

	int event_id, pdg_id;
	std::vector<int> *trackid_out_particle = 0;

	tree -> SetBranchAddress("m_event_id_MC", &event_id);
	tree -> SetBranchAddress("m_pdg_id", &pdg_id);
	tree -> SetBranchAddress("m_trackid_out_particle", &trackid_out_particle);

	for (int i=0; i<tree->GetEntries(); i++) {
		tree -> GetEntry(i);

		if (abs(pdg_id) != 14) continue;

		std::vector<int> trackid_buf;
		for (int j=0; j<trackid_out_particle->size(); j++) {
			trackid_buf.push_back(trackid_out_particle->at(j));
		}

		if (uniqueID.find(event_id+100000) == uniqueID.end()) {
			uniqueID[event_id+100000] = trackid_buf;
		} else {
			uniqueID[event_id+10000].insert(uniqueID[event_id+100000].end(), trackid_buf.begin(), trackid_buf.end());
		}
	}
}

bool isTrack(EdbTrackP* track) {
	int event_id = track -> GetSegmentFirst() -> MCEvt();
	int track_id = track -> GetSegmentFirst() -> Volume();

	//std::cout << "event id: " << event_id << "\ttrack id: " << track_id << std::endl;

	auto iter = uniqueID.find(event_id);

	if (iter == uniqueID.end()) {
		return false;
	} else {
		std::vector<int> v_trackid = iter -> second;
		if (std::find(v_trackid.begin(), v_trackid.end(), track_id) == v_trackid.end()) return false;
	}

	std::cout << "matching event ID, track ID: " << event_id << ", " << track_id << std::endl;
	return true;
}


TH1D* long_track_profile(std::string path) {

	const char delim = '_';
	std::stringstream ss;
	ss << path;
	std::string tmp;
	for (int i=0; i<=2; i++) {
		std::getline(ss, tmp, delim);
		std::cout << tmp << std::endl;
	}

	std::string tmp2 = "s.eMCEvt==";
	int event_id = atoi(tmp.c_str()) + 100000;
	tmp2 += std::to_string(event_id);

	TCut cut;
	cut += c;
	cut += tmp2.c_str();

	TH1D* hist = new TH1D("hist", path.c_str(), 2*bin_max, -bin_max, bin_max);
	
	EdbDataProc* dproc = new EdbDataProc;
	EdbPVRec* pvr = new EdbPVRec;
	dproc -> ReadTracksTree(*pvr, path.c_str(), cut); // read linked_tracks with cut npl >= 200.
	
	for (int i=0; i<pvr->Ntracks(); i++) {
		EdbTrackP* track = pvr -> GetTrack(i);

		if (isTrack(track)) {
			hist -> Fill(track->GetSegmentFirst() -> MCTrack());
			ntrack ++;
			Utils::PrintTrack(track);
		}
	}

	return hist;
}

int main() {

	for (int i=2; i<10; i++) {
		std::string filename = Form("20221202_nuall_200003-00000-dev/FaserMC-MC22_Genie_all_600fbInv_v1-200003-0000%d-dev-EM.root", i);
		initTruth(filename);
	}

	TRint app("app", 0, 0);

	std::string title = c.GetTitle();
	title += ";mrad;counts";

	TH1D* hist = new TH1D("hist", title.c_str(),2*bin_max, -bin_max, bin_max);
	
	// buffer of histograms
	TH1D* hist_buf;
	ntrack = 0;

	std::ifstream input_file;
	input_file.open("input_list_2.txt", std::ios::in);
	std::string line_buf;

	// read linked_tracks.root
	while (std::getline(input_file, line_buf)) {
		std::string path = line_buf + "/linked_tracks.root"; // path of linked_tracks.
		hist_buf = long_track_profile(path);
		std::cout << "Read: " << path << std::endl << std::endl << std::endl;
		hist -> Add(hist_buf);
	}

	std::cout << "number of long track: " << ntrack << std::endl;
	std::cout << "numeber of uniqueID: " << uniqueID.size() << std::endl;

	hist -> Draw();
	app.Run();

	return 0;
}
