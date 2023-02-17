#include <iostream>
#include <fstream>
#include <limits>
#include <ostream>
#include <sstream>
#include <string>

#include "EdbPVRec.h"
#include "Rtypes.h"
#include "TH1.h"
#include "TRint.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"

#include "EdbDataSet.h"

#include "Utils.hpp"
#include "TruthManager.hpp"

int bin_max = 310;
TCut cut = "nseg>3&&t.eP>20";
//TCut cut = "t.eP>20";
int n_long_mu = 0;
int n_long_pi = 0;
std::vector<std::string> v_path;

TruthManager* tm;
	

void make_mu_hist(EdbPVRec* pvr, int event_id, TH1D* mu_hist) {

	std::map<int, std::vector<EdbTrackP*>> primary_muons; // <trackID, EdbTrack>.

	for (int i=0; i<pvr->Ntracks(); i++) {
		EdbTrackP* track = pvr -> GetTrack(i);

		int pdg_id = track -> GetSegmentFirst() -> MCTrack();
		int track_id = track -> GetSegmentFirst() -> Volume();

		if (!tm -> IsTrack(track)) continue;

		if (event_id != track -> GetSegmentFirst() -> MCEvt()) continue;

		if (abs(track -> GetSegmentFirst() -> MCTrack()) == 13 and Utils::IsOutgo(track, pvr)) return;


		if (abs(pdg_id)==13) {
			primary_muons[track_id].push_back(track);
			//if (track -> Npl() < 3) Utils::PrintTrack(track);
			if (track -> Npl() > 200) {
				n_long_mu ++;
			} else {
				//v_path.push_back(std::to_string(event_id));
			}
		}
	}

	int plate_buf;
	for (auto iter: primary_muons) {
		std::vector<EdbTrackP*> tracks = iter.second;
		std::sort(tracks.begin(), tracks.end(),
				[](EdbTrackP* lhs, EdbTrackP* rhs) { return lhs -> GetSegmentFirst() -> PID() < rhs -> GetSegmentFirst() -> PID(); });

		int plate = Utils::HasKink(tracks.front());
		if (plate != -1) {
			v_path.push_back(std::to_string(event_id));
			std::cout << "Kink at " << plate << std::endl;
			mu_hist -> Fill(plate);
			continue;
		} else {
			mu_hist -> Fill(tracks.front() -> Npl());
		}


		//mu_hist -> Fill(tracks.front() -> Npl());
	}
	return ;
}

void make_pi_hist(EdbPVRec* pvr, int event_id, TH1D* pi_hist) {


	std::map<int, std::vector<EdbTrackP*>> primary_pions; // <trackID, EdbTrack>.

	for (int i=0; i<pvr->Ntracks(); i++) {
		EdbTrackP* track = pvr -> GetTrack(i);

		int pdg_id = track -> GetSegmentFirst() -> MCTrack();
		int track_id = track -> GetSegmentFirst() -> Volume();

		if (!tm -> IsTrack(track)) continue;

		if (event_id != track -> GetSegmentFirst() -> MCEvt()) continue;

		//if (abs(track -> GetSegmentFirst() -> MCTrack()) == 211 and Utils::IsOutgo(track, pvr)) continue;


		if (abs(pdg_id) == 211) {
			primary_pions[track_id].push_back(track);
			if (track -> Npl() > 200) n_long_pi ++;
		}
	}

	int plate_buf;
	for (auto iter: primary_pions) {
		std::vector<EdbTrackP*> tracks = iter.second;
		std::sort(tracks.begin(), tracks.end(),
				[](EdbTrackP* lhs, EdbTrackP* rhs) { return lhs -> GetSegmentFirst() -> PID() < rhs -> GetSegmentFirst() -> PID(); });
		int plate = Utils::HasKink(tracks.front());
		if (plate != -1) {
			pi_hist -> Fill(plate);
			std::cout << "Kink at " << plate << std::endl;
			continue;
		} else {
			pi_hist -> Fill(tracks.front() -> Npl());
		}


		//pi_hist -> Fill(tracks.front() -> Npl());
	}


	return ;
}


int main() {

	tm = new TruthManager("./truth_input.txt");

	TRint app("app", 0, 0);

	TH1D *mu_hist = new TH1D("mu hist", "", bin_max/2, 0, bin_max);
	TH1D *pi_hist = new TH1D("pi hist", "", bin_max/2, 0, bin_max);
	
	// buffer of histograms
	std::vector<TH1D*> hist_buf;

	std::ifstream input_file;
	input_file.open("./input_list_enlarged.txt", std::ios::in);
	std::string line_buf;
	int num_file = 0;

	EdbDataProc* dproc = new EdbDataProc;
	EdbPVRec* pvr = new EdbPVRec;

	// read linked_tracks.root
	while (std::getline(input_file, line_buf)) {
		std::string path = line_buf + "/linked_tracks.root"; // path of linked_tracks.

		// specify event ID.
		int event_id = 0;
		sscanf(path.c_str(), "20230203_nuall_enlarged/evt_%d", &event_id);
		event_id += 100000;
		//dproc -> ReadTracksTree(*pvr, path.c_str(), cut); // read linked_tracks with cut npl >= 200.
		pvr = Utils::ConnectTrack(path, 70, 0.007, "nseg>3&&t.eP>20");
	

		//track_length_hist(pvr, event_id, mu_hist, pi_hist);
		make_mu_hist(pvr, event_id, mu_hist);
		make_pi_hist(pvr, event_id, pi_hist);
		num_file ++;
		if (pvr->eTracks) pvr -> eTracks -> Clear();
	}

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
	mu_cum_hist -> SetTitle("");
	mu_cum_hist -> Draw();
	pi_cum_hist -> SetLineColor(kRed);
	pi_cum_hist -> Draw("SAME");

	mu_cum_hist -> GetXaxis() -> SetRangeUser(0, 200);
	mu_cum_hist -> GetYaxis() -> SetRangeUser(0, 1.4);

	cum_legend -> AddEntry(mu_cum_hist, "Muon");
	cum_legend -> AddEntry(pi_cum_hist, "Pion");
	cum_legend -> Draw();


	std::cout << "Number of muon: " << num_mu << std::endl;
	std::cout << "Number of pion: " << num_pi << std::endl;
	std::cout << "Number of muon npl > 200: " << mu_cum_hist -> GetBinContent(101) << std::endl;
	std::cout << "Number of pion npl > 200: " << pi_cum_hist -> GetBinContent(101) << std::endl;


	app.Run();

	return 0;
}
