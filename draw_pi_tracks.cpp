#include <iostream>
#include <fstream>
#include <set>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TRint.h"
#include "TCut.h"
#include "TObjArray.h"
#include "TStyle.h"

#include "EdbDataSet.h"

#include "HashTable.hpp"
#include "Utils.hpp"

TCut cut = "nseg>3";

std::map<int, std::vector<int>> uniqueID; // <evend ID, array(track_id)>
std::map<int, int> p_uniqueID;
std::set<std::string> s_u200nodup; // event include muon whose npl is smaller than 200 but which is not duplicated.
std::set<std::string> s_dupl; // event include muon whose npl is smaller than 200 but which is not duplicated.

TH1D* pos_hist = new TH1D("pos", ";micron;counts", 25, 0, 100);
TH1D* ang_hist = new TH1D("ang", ";mrad;counts", 100, 0, 100);
TH1I* pdg_hist = new TH1I("pdg", ";pdg;count", 600, -300, 300);
TH2D* pos_ang_hist = new TH2D("pos ang", ";micron;mrad", 50, 0, 50, 50, 0, 5);

void initTruth(std::string filename, int tree_number) {

	std::cout << "Read " << filename << std::endl;

	TFile* file = new TFile(filename.c_str(), "READ");
	TTree* tree = (TTree*) file -> Get("m_NuMCTruth_tree");

	int event_id, pdg_id, track_id;
	int num_in;
	std::vector<int> *trackid_out_particle = nullptr;
	std::vector<int> *pdg_out_particle= nullptr;

	tree -> SetBranchAddress("m_event_id_MC", &event_id);
	tree -> SetBranchAddress("m_track_id", &track_id);
	tree -> SetBranchAddress("m_pdg_id", &pdg_id);
	tree -> SetBranchAddress("m_trackid_out_particle", &trackid_out_particle);
	tree -> SetBranchAddress("m_pdg_out_particle", &pdg_out_particle);
	tree -> SetBranchAddress("m_num_in_particle", &num_in);

	for (int i=0; i<tree->GetEntries(); i++) {
		tree -> GetEntry(i);

		if (num_in != 0) continue;
		p_uniqueID[tree_number*1e5+event_id] = track_id;

		std::vector<int> trackid_buf;
		for (int j=0; j<trackid_out_particle->size(); j++) {
			trackid_buf.push_back(trackid_out_particle->at(j));
			//if (abs(pdg_out_particle->at(j)) == 13) num_nu_mu_cc ++;
		}

		if (uniqueID.find(tree_number*1e5+event_id) == uniqueID.end()) {
			// if event iD is not included uniqueID.
			uniqueID[tree_number*1e5+event_id] = trackid_buf;
		} else {
			// if event ID is included uniqueID -> add track ID.
			uniqueID[tree_number*1e5+event_id].insert(uniqueID[event_id+tree_number*1e5].end(), trackid_buf.begin(), trackid_buf.end());
		}
	}
	std::cout << p_uniqueID.size() << std::endl;
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
	return true;
}

bool isPrimaryTrack(EdbTrackP* track) {
	int event_id = track -> GetSegmentFirst() -> MCEvt();
	int track_id = track -> GetSegmentFirst() -> Volume();

	auto iter = p_uniqueID.find(event_id);
	if (iter == p_uniqueID.end()) {
		return false;
	} else if (iter -> second != track_id) {
		return false;
	}
	return true;

}


void make_hist(std::string path) {

	EdbDataProc* dproc = new EdbDataProc;
	EdbPVRec* pvr = new EdbPVRec;
	dproc -> ReadTracksTree(*pvr, path.c_str(), cut); // read linked_tracks with cut npl >= 200.

	std::cout << "Read: " << path << std::endl;

	std::map<int, EdbTrackP*> primary_tracks;
	std::map<int, std::vector<EdbTrackP*>> secondary_tracks;

	int num_primary = 0;

	for (int i=0; i<pvr->Ntracks(); i++) {
		EdbTrackP* track = pvr -> GetTrack(i);
		int event_id = track -> GetSegmentFirst() -> MCEvt();

		if (track->GetSegmentFirst()->Volume() == 10001 and track->GetSegmentFirst()->MCEvt()<1e5) {
			num_primary++;
		} 
		if (isPrimaryTrack(track)) {
			primary_tracks[event_id] = track;
			pdg_hist -> Fill(track->GetSegmentFirst()->MCTrack());
			//std::cout << "Mom: " << track -> GetSegmentFirst() -> P() << "\tTrack ID: " << track -> GetSegmentFirst() -> Volume() << "\tevent ID: " << track -> GetSegmentFirst() -> MCEvt() << std::endl;
		} else if (isTrack(track)) {
			secondary_tracks[event_id].push_back(track);
		} else {
			continue;
		}
	}

	for (auto iter: primary_tracks) {
		int event_id = iter.first;
		EdbTrackP* track = iter.second;
		
		auto secondary_iter = secondary_tracks.find(event_id);
		if (secondary_iter != secondary_tracks.end()) {
			std::vector<EdbTrackP*> cand_tracks = secondary_iter -> second;
			double min_dist = 1e9;
			double min_dtheta = 1e9;
			for (int i=0; i<cand_tracks.size(); i++) {
				EdbTrackP* cand_track = cand_tracks[i];
				
				double dist = Utils::Distance(track, cand_track);
				double dth = Utils::Dtheta(track, cand_track)*1000;

				if (dth < min_dtheta) {
					min_dist = dist;
					min_dtheta = dth;
				}

			}
			pos_hist -> Fill(min_dist);
			ang_hist -> Fill(min_dtheta);
			pos_ang_hist -> Fill(min_dist, min_dtheta);
		}
	}

	TObjArray* selected = new TObjArray();
	selected -> Add(primary_tracks[2]);

	for (int i=0; i<secondary_tracks[2].size(); i++) {
		selected -> Add(secondary_tracks[2][i]);
	}

	dproc -> MakeTracksTree(*selected, 0, 0, "hadrons.root");


	std::cout << "Number of Primary: " << num_primary << std::endl;
	return;
}

int main() {

	for (int i=0; i<1; i++) {
	std::string filename = Form("100032/ntp/s0008/FaserMC-MC22_PG_pion_in_fasernu_100GeV-100032-0000%d-s0008-NTUP.root", i);
	//std::string filename = Form("100033/ntp/s0008/FaserMC-MC22_PG_pion_in_fasernu_300GeV-100033-0000%d-s0008-NTUP.root", i);
		initTruth(filename, i);
	}

	TRint app("app", 0, 0);

	

	std::ifstream input_file;
	input_file.open("input_list_2.txt", std::ios::in);
	std::string line_buf;
	int num_file = 0;

	// read linked_tracks.root
	//std::string path = "100032/pseudo_reco/s0008/pseudo_linked_tracks.root"; // 100 GeV
	std::string path = "/data/FASER/MC/sim/mc22/particle_gun/100033/pseudo_reco/s0008/pseudo_linked_tracks.root";
	make_hist(path);


	TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
	TCanvas* c2 = new TCanvas("c2", "c2", 600, 600);
	TCanvas* c3 = new TCanvas("c3", "c3", 600, 600);
	TCanvas* c4 = new TCanvas("c4", "c4", 600, 600);

	c1 -> cd();
	pos_hist -> Draw();

	c2 -> cd();
	ang_hist -> Draw();

	c3 -> cd();
	pdg_hist -> Draw();

	c4 -> cd();
	pos_ang_hist -> SetMarkerStyle(8);
	pos_ang_hist -> Draw();

	gStyle -> SetOptStat(0);

	app.Run();

	return 0;
}
