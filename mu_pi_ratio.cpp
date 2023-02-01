#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>

#include "EdbPVRec.h"
#include "Rtypes.h"
#include "TH1.h"
#include "TRint.h"
#include "TLegend.h"
#include "TStyle.h"

#include "EdbDataSet.h"

#include "Utils.hpp"

int bin_max = 300;
TCut cut = "nseg>3";
int num_nu_mu_cc = 0;
int n_long_mu = 0;
std::vector<std::string> v_path;

std::map<int, std::vector<int>> uniqueID;

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

		if (abs(pdg_id) != 14) continue;

		std::vector<int> trackid_buf;
		for (int j=0; j<trackid_out_particle->size(); j++) {
			trackid_buf.push_back(trackid_out_particle->at(j));
			//if (abs(pdg_out_particle->at(j)) == 13) num_nu_mu_cc ++;
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

	//std::cout << "matching event ID, track ID: " << event_id << ", " << track_id << std::endl;
	return true;
}

std::pair<TH1D*, TH1D*> track_length_hist(std::string path) {

	// specify event ID.
	int event_id = 0;
	sscanf(path.c_str(), "20230107_nuall/evt_%d", &event_id);
	event_id += 100000;
	
	EdbDataProc* dproc = new EdbDataProc;
	EdbPVRec* pvr = new EdbPVRec;
	dproc -> ReadTracksTree(*pvr, path.c_str(), cut); // read linked_tracks with cut npl >= 200.
	
	//pvr = Utils::ConnectTrack(path);
	

	TH1D* mu_hist = new TH1D("mu_hist", "mu_hist", bin_max+1, -0.5, bin_max+0.5);
	TH1D* pi_hist = new TH1D("pi_hist", "pi_hist", bin_max+1, -0.5, bin_max+0.5);

	std::cout << "Event : " << event_id << std::endl;

	std::cout << "Read: " << path << std::endl;

	for (int i=0; i<pvr->Ntracks(); i++) {
		EdbTrackP* track = pvr -> GetTrack(i);

		int pdg_id = track -> GetSegmentFirst() -> MCTrack();

		if (!isTrack(track)) continue;

		if (event_id != track -> GetSegmentFirst() -> MCEvt()) continue;

		if (track -> GetSegmentFirst() -> ScanID().GetPlate() > 100) continue;

		if (abs(pdg_id)==13) {
			mu_hist -> Fill(track->Npl());
			if (track -> Npl() < 3) Utils::PrintTrack(track);
			if (track -> Npl() > 200) {
				n_long_mu ++;
			} else {
				v_path.push_back(path);
			}
		} else if (abs(pdg_id) == 211) {
			pi_hist -> Fill(track->Npl());
		}
	}

	return {mu_hist, pi_hist};
}


int main() {

	for (int i=0; i<10; i++) {
		std::string filename = Form("20221202_nuall_200003-00000-dev/FaserMC-MC22_Genie_all_600fbInv_v1-200003-0000%d-dev-EM.root", i);
		initTruth(filename);
	}

	TRint app("app", 0, 0);

	std::string title = std::string(cut.GetTitle());
	title += ";npl;counts";
	title = ";npl;counts";

	TH1D* mu_hist = new TH1D("mu hist", title.c_str(), bin_max+1, -0.5, bin_max+0.5);
	TH1D* pi_hist = new TH1D("pi hist", title.c_str(), bin_max+1, -0.5, bin_max+0.5);
	
	// buffer of histograms
	std::pair<TH1D*, TH1D*> hist_buf;

	std::ifstream input_file;
	input_file.open("input_list_2.txt", std::ios::in);
	std::string line_buf;
	int num_file = 0;

	// read linked_tracks.root
	while (std::getline(input_file, line_buf)) {
		std::string path = line_buf + "/linked_tracks.root"; // path of linked_tracks.
		hist_buf = track_length_hist(path);
		mu_hist -> Add(hist_buf.first);
		pi_hist -> Add(hist_buf.second);
		num_file ++;
	}

	std::cout << "numeber of uniqueID: " << uniqueID.size() << std::endl;
	std::cout << num_file << " files are read." << std::endl;
	std::cout << "Number of NuMuCC: " << num_nu_mu_cc << std::endl;

	mu_hist -> GetYaxis() -> SetRangeUser(0, 450);

	TLegend* legend = new TLegend(0.8,  0.8, 0.9, 0.9);
	legend -> AddEntry(mu_hist, "muon");
	legend -> AddEntry(pi_hist, "pion");

	mu_hist -> GetCumulative(kFALSE) -> Draw();
	pi_hist -> SetLineColor(kRed);
	pi_hist -> GetCumulative(kFALSE) -> Draw("SAME");

	for (auto p: v_path) {
		std::cout << p << std::endl;
	}

	std::cout << "Number of muon penetrate 200 plate: " << mu_hist -> GetCumulative(kFALSE) -> GetBinContent(200) << std::endl;
	std::cout << "Number of pion penetrate 200 plate: " << pi_hist -> GetCumulative(kFALSE) -> GetBinContent(200) << std::endl;
	std::cout << "Number of muon npl > 200: " << n_long_mu << std::endl;

	gStyle -> SetOptStat(0);

	legend -> Draw();
	app.Run();

	return 0;
}
