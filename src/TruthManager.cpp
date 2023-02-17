#include "TruthManager.hpp"
#include <cstdlib>
#include <fstream>
#include <string>

TruthManager::TruthManager(std::string input_path) : truth_path_(input_path) {

	std::cout << "Read truth files." << std::endl;

	std::ifstream truth_input_file;
	truth_input_file.open(input_path);

	if (!truth_input_file) {
		std::cerr << "Invalid file!" << std::endl;
		exit(1);
	}

	std::string line_buf;
	while (std::getline(truth_input_file, line_buf)) {
		ReadTruth(line_buf);
	}

}

TruthManager::~TruthManager() {}

void TruthManager::ReadTruth(std::string path) {
	
	std::cout << "Read " << path << std::endl;

	TFile* file = new TFile(path.c_str(), "READ");
	TTree* tree = (TTree*) file -> Get("m_NuMCTruth_tree");

	int event_id, pdg_id;
	float vz_decay;
	std::vector<int> *trackid_out_particle = 0;
	std::vector<int> *pdg_out_particle= 0;

	tree -> SetBranchAddress("m_event_id_MC", &event_id);
	tree -> SetBranchAddress("m_pdg_id", &pdg_id);
	tree -> SetBranchAddress("m_trackid_out_particle", &trackid_out_particle);
	tree -> SetBranchAddress("m_pdg_out_particle", &pdg_out_particle);
	tree -> SetBranchAddress("m_vz_decay", &vz_decay);

	for (int i=0; i<tree->GetEntries(); i++) {
		tree -> GetEntry(i);

		if (abs(pdg_id) != 14) continue;

		std::vector<int> trackid_buf;
		for (int j=0; j<trackid_out_particle->size(); j++) {
			trackid_buf.push_back(trackid_out_particle->at(j));
		}

		if (unique_id_.find(event_id+100000) == unique_id_.end()) {
			unique_id_[event_id+100000] = trackid_buf;
		} else {
			unique_id_[event_id+100000].insert(unique_id_[event_id+100000].end(), trackid_buf.begin(), trackid_buf.end());
		}
	}
}

bool TruthManager::IsTrack(EdbTrackP* track) {

	int event_id = track -> GetSegmentFirst() -> MCEvt();
	int track_id = track -> GetSegmentFirst() -> Volume();
	double z = track -> GetSegmentFirst() -> Z();

	//std::cout << "event id: " << event_id << "\ttrack id: " << track_id << std::endl;

	auto iter = unique_id_.find(event_id);

	if (iter == unique_id_.end()) {
		return false;
	} else {
		std::vector<int> v_trackid = iter -> second;
		if (std::find(v_trackid.begin(), v_trackid.end(), track_id) == v_trackid.end()) return false;
	}

	//std::cout << "matching event ID, track ID: " << event_id << ", " << track_id << std::endl;
	return true;
}
