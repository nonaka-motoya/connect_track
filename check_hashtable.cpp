#include <iostream>

#include <TObjArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TRint.h>

#include <EdbDataSet.h>

#include "HashTable.hpp"
#include "Utils.hpp"

int main() {

	TRint app("app", 0, 0);

	EdbDataProc* dproc = new EdbDataProc;
	EdbPVRec* pvr = new EdbPVRec;

	dproc -> ReadTracksTree(*pvr, "20230107_nuall/evt_2503_pl1_300/linked_tracks.root");

	int iplmax = -1e9;
	int iplmin = 1e9;
	double xmin = 1e9;
	double xmax = -1e9;
	double ymin = 1e9;
	double ymax = -1e9;

	std::vector<EdbTrackP*> all_tracks;
	
	for (int i = 0; i<pvr->Ntracks(); i++) {
		EdbTrackP* track = pvr -> GetTrack(i);
		all_tracks.push_back(track);

		for (int j=0; j<track->N(); j++) {
			EdbSegP* seg = track -> GetSegment(j);
			int sPlate = seg ->Plate();
			float sX = seg ->X();
			float sY = seg ->Y();

			if (iplmin > sPlate) iplmin = sPlate;
			if (iplmax < sPlate) iplmax = sPlate;
			if (xmin > sX) xmin = sX;
			if (xmax < sX) xmax = sX;
			if (ymin > sY) ymin = sY;
			if (ymax < sY) ymax = sY;

		}

	}
	xmin -= 1; xmax += 1; ymin -= 1; ymax += 1;
	printf("ipl %d - %d\n", iplmin, iplmax);


	HashTable* hashtable = new HashTable(pvr);

	TObjArray* selected = new TObjArray;

	TH1D* pos_hist = new TH1D("pos", "", 100, 0, 100);
	TH1D* ang_hist = new TH1D("ang", "", 100, 0, 100);
	TH2D* pos_ang_hist = new TH2D("pos ang", "", 100, 0, 100, 100, 0, 100);

	EdbTrackP* track;

	track = pvr -> FindTrack(10650);
	selected -> Add(track);


	std::vector<EdbTrackP*> tracks;


	double r[3] = {6.1, 20, 20};

	tracks = hashtable -> GetNeighbors(track);
	for (int i=0; i<tracks.size(); i++) {
		EdbTrackP* cand_track = tracks[i];

		double dist = Utils::Distance(track, cand_track);
		double dtheta = Utils::Dtheta(track, cand_track) * 1000;

		pos_hist -> Fill(dist);
		ang_hist -> Fill(dtheta);
		pos_ang_hist -> Fill(dist, dtheta);

		selected -> Add(tracks[i]);
	}

	dproc -> MakeTracksTree(*selected, 0, 0, "neighbors.root");

	TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
	TCanvas* c2 = new TCanvas("c2", "c2", 600, 600);
	TCanvas* c3 = new TCanvas("c3", "c3", 600, 600);

	c1 -> cd();
	pos_hist -> Draw();

	c2 -> cd();
	ang_hist -> Draw();

	c3 -> cd();
	pos_ang_hist -> Draw("COLZ");

	app.Run();

	return 0;
}
