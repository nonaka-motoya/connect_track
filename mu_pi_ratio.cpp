#include <iostream>
#include <fstream>
#include <limits>
#include <ostream>
#include <sstream>

#include "EdbPVRec.h"
#include "Rtypes.h"
#include "TH1.h"
#include "TRint.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"

#include "EdbDataSet.h"

#include "Utils.hpp"

int bin_max = 310;
TCut cut = "nseg>3";
int n_long_mu = 0;
int n_long_pi = 0;
std::vector<std::string> v_path;

std::map<int, std::vector<int>> uniqueID;

	
void initTruth(std::string filename) {

	std::cout << "Read " << filename << std::endl;

	TFile* file = new TFile(filename.c_str(), "READ");
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

		if (uniqueID.find(event_id+100000) == uniqueID.end()) {
			uniqueID[event_id+100000] = trackid_buf;
		} else {
			uniqueID[event_id+100000].insert(uniqueID[event_id+100000].end(), trackid_buf.begin(), trackid_buf.end());
		}
	}
}

bool isTrack(EdbTrackP* track) {
	int event_id = track -> GetSegmentFirst() -> MCEvt();
	int track_id = track -> GetSegmentFirst() -> Volume();
	double z = track -> GetSegmentFirst() -> Z();

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

void track_length_hist(EdbPVRec* pvr, int event_id, TH1D* mu_hist, TH1D* pi_hist) {
	std::cout << "Event : " << event_id << std::endl;

	std::map<int, std::vector<EdbTrackP*>> primary_muons; // <trackID, EdbTrack>.
	std::map<int, std::vector<EdbTrackP*>> primary_pions; // <trackID, EdbTrack>.

	for (int i=0; i<pvr->Ntracks(); i++) {
		EdbTrackP* track = pvr -> GetTrack(i);

		int pdg_id = track -> GetSegmentFirst() -> MCTrack();
		int track_id = track -> GetSegmentFirst() -> Volume();

		if (!isTrack(track)) continue;

		if (event_id != track -> GetSegmentFirst() -> MCEvt()) continue;

		if (abs(track -> GetSegmentFirst() -> MCTrack()) == 13 and Utils::IsOutgo(track, pvr)) return;

		if (abs(pdg_id)==13) {
			primary_muons[track_id].push_back(track);
			//if (track -> Npl() < 3) Utils::PrintTrack(track);
			if (track -> Npl() > 200) {
				n_long_mu ++;
			} else {
				//v_path.push_back(path);
			}
		} else if (abs(pdg_id) == 211) {
			primary_pions[track_id].push_back(track);
			if (track -> Npl() > 200) n_long_pi ++;
		}
	}

	for (auto iter: primary_muons) {
		std::vector<EdbTrackP*> tracks = iter.second;
		std::sort(tracks.begin(), tracks.end(),
				[](EdbTrackP* lhs, EdbTrackP* rhs) { return lhs -> GetSegmentFirst() -> PID() < rhs -> GetSegmentFirst() -> PID(); });
		mu_hist -> Fill(tracks.front() -> Npl());
	}

	for (auto iter: primary_pions) {
		std::vector<EdbTrackP*> tracks = iter.second;
		std::sort(tracks.begin(), tracks.end(),
				[](EdbTrackP* lhs, EdbTrackP* rhs) { return lhs -> GetSegmentFirst() -> PID() < rhs -> GetSegmentFirst() -> PID(); });
		pi_hist -> Fill(tracks.front() -> Npl());
	}


	return ;
}

void make_mu_hist(EdbPVRec* pvr, int event_id, TH1D* mu_hist) {

	std::map<int, std::vector<EdbTrackP*>> primary_muons; // <trackID, EdbTrack>.

	for (int i=0; i<pvr->Ntracks(); i++) {
		EdbTrackP* track = pvr -> GetTrack(i);

		int pdg_id = track -> GetSegmentFirst() -> MCTrack();
		int track_id = track -> GetSegmentFirst() -> Volume();

		if (!isTrack(track)) continue;

		if (event_id != track -> GetSegmentFirst() -> MCEvt()) continue;

		if (abs(track -> GetSegmentFirst() -> MCTrack()) == 13 and Utils::IsOutgo(track, pvr)) return;

		if (abs(pdg_id)==13) {
			primary_muons[track_id].push_back(track);
			//if (track -> Npl() < 3) Utils::PrintTrack(track);
			if (track -> Npl() > 200) {
				n_long_mu ++;
			} else {
				//v_path.push_back(path);
			}
		}
	}

	for (auto iter: primary_muons) {
		std::vector<EdbTrackP*> tracks = iter.second;
		std::sort(tracks.begin(), tracks.end(),
				[](EdbTrackP* lhs, EdbTrackP* rhs) { return lhs -> GetSegmentFirst() -> PID() < rhs -> GetSegmentFirst() -> PID(); });
		mu_hist -> Fill(tracks.front() -> Npl());
	}

	return ;
}

void make_pi_hist(EdbPVRec* pvr, int event_id, TH1D* pi_hist) {


	std::map<int, std::vector<EdbTrackP*>> primary_pions; // <trackID, EdbTrack>.

	for (int i=0; i<pvr->Ntracks(); i++) {
		EdbTrackP* track = pvr -> GetTrack(i);

		int pdg_id = track -> GetSegmentFirst() -> MCTrack();
		int track_id = track -> GetSegmentFirst() -> Volume();

		if (!isTrack(track)) continue;

		if (event_id != track -> GetSegmentFirst() -> MCEvt()) continue;

		if (abs(pdg_id) == 211) {
			primary_pions[track_id].push_back(track);
			if (track -> Npl() > 200) n_long_pi ++;
		}
	}

	for (auto iter: primary_pions) {
		std::vector<EdbTrackP*> tracks = iter.second;
		std::sort(tracks.begin(), tracks.end(),
				[](EdbTrackP* lhs, EdbTrackP* rhs) { return lhs -> GetSegmentFirst() -> PID() < rhs -> GetSegmentFirst() -> PID(); });
		pi_hist -> Fill(tracks.front() -> Npl());
	}


	return ;
}


int main() {

	for (int i=0; i<10; i++) {
		std::string filename = Form("20221202_nuall_200003-00000-dev/FaserMC-MC22_Genie_all_600fbInv_v1-200003-0000%d-dev-EM.root", i);
		initTruth(filename);
	}

	TRint app("app", 0, 0);

	std::string title = std::string(cut.GetTitle());
	title += ";#plate;#tracks";
	title = ";#plate;#tracks";

	TH1D *mu_hist = new TH1D("mu hist", title.c_str(), bin_max/2, 0, bin_max);
	TH1D *pi_hist = new TH1D("pi hist", title.c_str(), bin_max/2, 0, bin_max);
	
	// buffer of histograms
	std::vector<TH1D*> hist_buf;

	std::ifstream input_file;
	input_file.open("input_list_2.txt", std::ios::in);
	std::string line_buf;
	int num_file = 0;

	// read linked_tracks.root
	while (std::getline(input_file, line_buf)) {
		std::string path = line_buf + "/linked_tracks.root"; // path of linked_tracks.

		EdbDataProc* dproc = new EdbDataProc;
		EdbPVRec* pvr = new EdbPVRec;

		// specify event ID.
		int event_id = 0;
		sscanf(path.c_str(), "20230107_nuall/evt_%d", &event_id);
		event_id += 100000;
		//dproc -> ReadTracksTree(*pvr, path.c_str(), cut); // read linked_tracks with cut npl >= 200.
		
		pvr = Utils::ConnectTrack(path);
	

		//track_length_hist(pvr, event_id, mu_hist, pi_hist);
		make_mu_hist(pvr, event_id, mu_hist);
		make_pi_hist(pvr, event_id, pi_hist);
		num_file ++;
	}

	std::cout << "numeber of uniqueID: " << uniqueID.size() << std::endl;
	std::cout << num_file << " files are read." << std::endl;

	TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);

	for (auto p: v_path) {
		std::cout << p << std::endl;
	}

	gStyle -> SetOptStat(0);
	// cumulative histograms.
	c1 -> cd();
	int num_mu = mu_hist -> GetEntries();
	int num_pi = pi_hist -> GetEntries();

	TH1 *mu_cum_hist;
	TH1 *pi_cum_hist;

	TLegend* cum_legend = new TLegend(0.7,  0.7, 0.9, 0.9);
	mu_cum_hist = mu_hist -> GetCumulative(kFALSE);
	pi_cum_hist = pi_hist -> GetCumulative(kFALSE);
	mu_cum_hist -> Scale(1./num_mu);
	pi_cum_hist -> Scale(1./num_pi);
	mu_cum_hist -> SetTitle(";#plate;#tracks (a.u.)");
	mu_cum_hist -> Draw();
	pi_cum_hist -> SetLineColor(kRed);
	pi_cum_hist -> Draw("SAME");
	cum_legend -> Draw();

	cum_legend -> AddEntry(mu_cum_hist, "Muon");
	cum_legend -> AddEntry(pi_cum_hist, "Pion");
	cum_legend -> Draw();


	std::cout << "Number of muon: " << num_mu << std::endl;
	std::cout << "Number of muon npl > 200: " << n_long_mu << std::endl;
	std::cout << "Number of pion npl > 200: " << n_long_pi << std::endl;


	app.Run();

	return 0;
}
