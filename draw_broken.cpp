#include <iostream>
#include <fstream>
#include <set>

#include "TCanvas.h"
#include "TH1.h"
#include "TRint.h"
#include "TCut.h"
#include "TStyle.h"

#include "Utils.hpp"
#include "EdbDataSet.h"

TCut cut = "nseg>3&&t.eP>20";

std::map<int, std::vector<int>> uniqueID; // <evend ID, array(track_id)>.
std::set<std::string> s_u200nodup; // event include muon whose npl is smaller than 200 but which is not duplicated.
std::set<std::string> s_dupl; // event include muon whose npl is smaller than 200 but which is not duplicated.
							  //
EdbDataProc* dproc = new EdbDataProc;
EdbPVRec* pvr = new EdbPVRec;

TH1D* pos_hist = new TH1D("pos", ";micron;counts", 25, 0, 100);
TH1D* ang_hist = new TH1D("ang", ";mrad;counts", 25, 0, 10);
TH2D* mom_ang_hist = new TH2D("mom ang", ";mrad;Mom (GeV)", 50, 0, 10, 100, 0, 2000);
TH2D* pos_ang_hist = new TH2D("pos ang", ";micron;mrad", 50, 0, 50, 50, 0, 5);
TH1D* broken_hist = new TH1D("broken place", ";plate;counts", 300, 0, 300);

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
void make_hist(std::string path) {
	dproc -> ReadTracksTree(*pvr, path.c_str(), cut); // read linked_tracks with cut npl >= 200.
	//pvr = Utils::ConnectTrack(path);

	// specify event ID from file name.
	int event_id = 0;
	sscanf(path.c_str(), "20230107_nuall/evt_%d", &event_id);
	event_id += 100000;
	
	std::cout << "Event : " << event_id << std::endl;
	std::cout << "Read: " << path << std::endl;

	// map for count duplicated track <MC track  ID, number of tracks>


	std::vector<EdbTrackP*> v_muons;
	for (int i=0; i<pvr->Ntracks(); i++) {
		EdbTrackP* track = pvr -> GetTrack(i);

		// is track outgoing from nu mu.
		if (!isTrack(track)) continue;

		// is track in specific event.
		if (event_id != track -> GetSegmentFirst() -> MCEvt()) continue;

		int pdg_id = track -> GetSegmentFirst() -> MCTrack();
		if (abs(pdg_id) == 13) v_muons.push_back(track);
	}

	std::sort(v_muons.begin(), v_muons.end(),
			[] (const EdbTrackP* lhs, const EdbTrackP* rhs) {return lhs->GetSegmentFirst()->ScanID().GetPlate() < rhs->GetSegmentFirst()->ScanID().GetPlate();}
			);

	std::cout << "number of muon: " << v_muons.size() << std::endl;
	if (v_muons.size() == 0) return;
	for (int i=0; i<v_muons.size()-1; i++) {
		EdbTrackP* prv_track = v_muons[i];
		EdbTrackP* nxt_track = v_muons[i+1];

		double dist = Utils::Distance(prv_track, nxt_track);
		double ang = Utils::Dtheta(prv_track, nxt_track)*1000;
		double mom = prv_track -> GetSegmentFirst() -> P();

		std::cout << "Broken at " << prv_track -> GetSegmentLast() -> ScanID().GetPlate() << std::endl;
		std::cout << "distance: " << dist << " micron." << "\tangle: " << ang << " mrad" << "\tMom: " << mom << " GeV" << std::endl;

		pos_hist -> Fill(dist);
		ang_hist -> Fill(ang);
		mom_ang_hist -> Fill(ang, mom);
		pos_ang_hist -> Fill(dist, ang);
		broken_hist -> Fill(prv_track -> GetSegmentLast() -> ScanID().GetPlate());
	}

	return;
}

int main() {

	for (int i=0; i<10; i++) {
		std::string filename = Form("20221202_nuall_200003-00000-dev/FaserMC-MC22_Genie_all_600fbInv_v1-200003-0000%d-dev-EM.root", i);
		initTruth(filename);
	}

	TRint app("app", 0, 0);

	std::ifstream input_file("broken_muon.csv");
	std::string line_buf;
	int num_file = 0;

	// read linked_tracks.root
	while (std::getline(input_file, line_buf)) {
		std::string path = line_buf; // path of linked_tracks.
		make_hist(path);
		num_file ++;
	}

	TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
	TCanvas* c2 = new TCanvas("c2", "c2", 600, 600);
	TCanvas* c3 = new TCanvas("c3", "c3", 600, 600);
	TCanvas* c4 = new TCanvas("c4", "c4", 600, 600);
	TCanvas* c5 = new TCanvas("c5", "c5", 600, 600);

	c1 -> cd();
	pos_hist -> Draw();

	c2 -> cd();
	ang_hist -> Draw();

	c3 -> cd();
	mom_ang_hist -> SetMarkerStyle(8);
	mom_ang_hist -> Draw();

	c4 -> cd();
	pos_ang_hist -> SetMarkerStyle(8);
	pos_ang_hist -> Draw();

	c5 -> cd();
	broken_hist -> Draw();

	gStyle -> SetOptStat(0);

	app.Run();

	return 0;
}
