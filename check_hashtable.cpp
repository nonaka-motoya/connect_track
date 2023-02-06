#include <iostream>

#include <TObjArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TRint.h>

#include <EdbDataSet.h>

#include "HashTable3D.hpp"
#include "Utils.hpp"

int main() {

	TRint app("app", 0, 0);

	EdbDataProc* dproc = new EdbDataProc;
	EdbPVRec* pvr = new EdbPVRec;

	dproc -> ReadTracksTree(*pvr, "./20221024_dset_TS_zone3_p171-210_finealign_p171-210_pm11000_from_pc36/linked_tracks.root");

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

	double cellsize = 20; // 50 micron

	HashTable3D* hashtable = new HashTable3D(
			pvr,
			iplmax - iplmin + 1, iplmin - 0.5, iplmax + 0.5,
			floor((xmax - xmin)/cellsize) + 1, xmin, xmax,
			floor((ymax - ymin)/cellsize) + 1, ymin, ymax
			);

	hashtable -> fillCells(all_tracks);


	TObjArray* selected = new TObjArray;

	TH1D* pos_hist = new TH1D("pos", "", 100, 0, 100);
	TH1D* ang_hist = new TH1D("ang", "", 100, 0, 100);
	TH2D* pos_ang_hist = new TH2D("pos ang", "", 100, 0, 100, 100, 0, 100);

	EdbTrackP* track;

	track = pvr -> GetTrack(0);
	selected -> Add(track);


	std::vector<EdbTrackP*> tracks;

	double v[3];
	v[0] = track -> GetSegmentLast() -> ScanID().GetPlate();
	v[1] = track -> GetSegmentLast() -> X();
	v[2] = track -> GetSegmentLast() -> Y();

	double r[3] = {6.1, 20, 20};

	tracks = hashtable -> GetNeighbors(v, r);
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
