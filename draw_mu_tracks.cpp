#include <iostream>
#include <fstream>
#include <set>

#include "TCanvas.h"
#include "TH1.h"
#include "TRint.h"
#include "TCut.h"

#include "Utils.hpp"
#include "EdbDataSet.h"

TCut cut = "nseg>3";

std::map<int, std::vector<int>> uniqueID; // <evend ID, array(track_id)>.
std::set<std::string> s_u200nodup; // event include muon whose npl is smaller than 200 but which is not duplicated.
std::set<std::string> s_dupl; // event include muon whose npl is smaller than 200 but which is not duplicated.

TH1I* dup_hist = new TH1I("num duplicated", ";#duplicate;counts", 10+1, -0.5, 10+0.5);
TH1D* mom_hist = new TH1D("mom", ";mom (GeV);counts", 100, 0, 10000);
TH2D* dup_mom_hist = new TH2D("mom", ";#duplicate;mom (GeV)",10, 0, 10, 100, 0, 1000);

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


// count duplicated track.
std::pair<TH1I*, TH2I*> make_hist(std::string path) {
	TH1I* hist1 = new TH1I("h1", "h1", 10+1, -0.5, 10+0.5);
	TH2I* hist2 = new TH2I("h2", "h2", 10+1, -0.5, 10+0.5, 300+1, -0.5, 300+0.5);

	EdbDataProc* dproc = new EdbDataProc;
	EdbPVRec* pvr = new EdbPVRec;
	dproc -> ReadTracksTree(*pvr, path.c_str(), cut); // read linked_tracks with cut npl >= 200.
	//pvr = Utils::ConnectTrack(path);

	// specify event ID from file name.
	int event_id = 0;
	sscanf(path.c_str(), "20230107_nuall/evt_%d", &event_id);
	//sscanf(path.c_str(), "20230203_nuall_enlarged/evt_%d", &event_id);
	event_id += 100000;
	
	std::cout << "Event : " << event_id << std::endl;
	std::cout << "Read: " << path << std::endl;

	// map for count duplicated track <MC track  ID, number of tracks>
	std::map<int, int> mu_num_duplicate_tracks;

	bool over200 = false;

	for (int i=0; i<pvr->Ntracks(); i++) {
		EdbTrackP* track = pvr -> GetTrack(i);

		// is track outgoing from nu mu.
		if (!isTrack(track)) continue;

		// is track in specific event.
		if (event_id != track -> GetSegmentFirst() -> MCEvt()) continue;

		int f_pdg_id = track -> GetSegmentFirst() -> MCTrack();
		int c_pdg_id = track -> GetSegment(track->N()/2) -> MCTrack();
		int track_id = track -> GetSegmentFirst() -> Volume();
		int npl = track -> Npl();
		double mom = track -> GetSegmentFirst() -> P();

		if (npl>200) over200 = true;

		if (abs(f_pdg_id) == 13 and abs(c_pdg_id)==13) {
			auto iter = mu_num_duplicate_tracks.find(track_id);

			// if track is not found.
			if (iter == mu_num_duplicate_tracks.end()) {
				mu_num_duplicate_tracks[track_id] = 1;
				mom_hist -> Fill(mom);
			} else {
				mu_num_duplicate_tracks[track_id] ++;
			}
		}
	}

	for (auto iter = mu_num_duplicate_tracks.begin(); iter != mu_num_duplicate_tracks.end(); ++iter) {
		dup_hist -> Fill(iter -> second);
		if (iter -> second == 1 and !over200) {
			s_u200nodup.insert(path);
		} else if (iter -> second == 2) {
			s_dupl.insert(path);
		}
	}

	return {hist1, hist2};
}

int main() {

	for (int i=0; i<10; i++) {
		std::string filename = Form("20221202_nuall_200003-00000-dev/FaserMC-MC22_Genie_all_600fbInv_v1-200003-0000%d-dev-EM.root", i);
		initTruth(filename);
	}

	TRint app("app", 0, 0);

	std::string title = std::string(cut.GetTitle());
	title += ";#duplicated;counts";


	
	// buffer of histograms
	std::pair<TH1I*, TH2I*> hist_buf;

	std::ifstream input_file;
	input_file.open("input_list_2.txt", std::ios::in);
	//input_file.open("input_list_enlarged.txt", std::ios::in);
	std::string line_buf;
	int num_file = 0;

	// read linked_tracks.root
	while (std::getline(input_file, line_buf)) {
		std::string path = line_buf + "/linked_tracks.root"; // path of linked_tracks.
		hist_buf = make_hist(path);
		num_file ++;
	}

	//std::ofstream ofs("broken_muon.csv");

	std::cout << "========================================" << std::endl;
	std::cout << "Number of event include muon whose npl < 200 but which is not duplicated: " << s_u200nodup.size() << std::endl;
	for (auto p: s_u200nodup) {
		std::cout << p << std::endl;
	}

	std::cout << "========================================" << std::endl;
	std::cout << "Number of event include duplicated muon track: " << s_dupl.size() << std::endl;
	for (auto p: s_dupl) {
		std::cout << p << std::endl;
		//ofs << p << std::endl;
	}
	std::cout << "========================================" << std::endl;

	TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
	TCanvas* c2 = new TCanvas("c2", "c2", 600, 600);

	c1 -> cd();
	dup_hist -> Draw();

	c2 -> cd();
	mom_hist -> Draw();

	app.Run();

	return 0;
}
