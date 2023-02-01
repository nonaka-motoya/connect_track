#include <algorithm>
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <string>
#include <set>

#include "EdbEDAUtil.h"
#include "EdbPVRec.h"
#include "Rtypes.h"
#include "TH1.h"
#include "TObjArray.h"
#include "TRint.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TString.h"
#include "TSystem.h"

#include "EdbDataSet.h"

#include "Utils.hpp"
#include "HashTable.hpp"

int bin_max = 10;
TCut cut = "nseg>3";

std::map<int, std::vector<int>> uniqueID; // <evend ID, array(track_id)>.


// fill unique ID using truth information.
void initTruth(std::string filename) {

	std::cout << "Read " << filename << std::endl;

	TFile* file = new TFile(filename.c_str(), "READ");
	TTree* tree = (TTree*) file -> Get("m_NuMCTruth_tree");

	int event_id, pdg_id;
	std::vector<int> *trackid_out_particle = 0;
	std::vector<int> *pdg_out_particle= 0;

	tree -> SetBranchAddress("m_event_id_MC", &event_id);
	tree -> SetBranchAddress("m_pdg_id", &pdg_id);
	tree -> SetBranchAddress("m_trackid_out_particle", &trackid_out_particle);
	tree -> SetBranchAddress("m_pdg_out_particle", &pdg_out_particle);

	for (int i=0; i<tree->GetEntries(); i++) {
		tree -> GetEntry(i);

		// is nu mu?
		if (abs(pdg_id) != 14) continue;

		std::vector<int> trackid_buf;
		for (int j=0; j<trackid_out_particle->size(); j++) {
			trackid_buf.push_back(trackid_out_particle->at(j));
			//if (abs(pdg_out_particle->at(j)) == 13) num_nu_mu_cc ++;
		}

		if (uniqueID.find(event_id+100000) == uniqueID.end()) {
			// if event iD is not included uniqueID.
			uniqueID[event_id+100000] = trackid_buf;
		} else {
			// if event ID is included uniqueID -> add track ID.
			uniqueID[event_id+10000].insert(uniqueID[event_id+100000].end(), trackid_buf.begin(), trackid_buf.end());
		}
	}
}

// check is track included uniqueID.
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

	//std::cout << "matching event ID, track ID: " << event_id << ", " << track_id << std::endl;
	return true;
}


int main() {

	for (int i=0; i<3; i++) {
		std::string filename = Form("20221202_nuall_200003-00000-dev/FaserMC-MC22_Genie_all_600fbInv_v1-200003-0000%d-dev-EM.root", i);
		initTruth(filename);
	}

	TRint app("app", 0, 0);


	std::ifstream input_file;
	input_file.open("input_list_2.txt", std::ios::in);
	std::string line_buf;
	int num_file = 0;

	// read linked_tracks.root
	while (std::getline(input_file, line_buf)) {
		std::string path = line_buf + "/linked_tracks.root"; // path of linked_tracks.
		num_file ++;
	}

	// development
	//std::string path = "20230107_nuall/evt_2950_pl1_300/linked_tracks.root"; // path of linked_tracks.
	std::string path = "20230107_nuall/evt_2503_pl1_300/linked_tracks.root"; // path of linked_tracks.
	Utils::ConnectTrack(path);

	app.Run();

	return 0;
}
