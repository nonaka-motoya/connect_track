#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <ostream>

#include <TCut.h>

#include <EdbDataSet.h>

#include "TruthManager.hpp"
#include "Utils.hpp"

class ClassifyNumu : private TruthManager {
	private:
		std::string input_file_path_;
		std::vector<std::string> input_files_;
		std::string event_id_prefix_;

		std::string output_file_path_;
		std::ofstream ofs_;

		EdbDataProc* dproc_;
		EdbPVRec* pvr_;

		TCut cut_;

	public:
		// Constructor.
		ClassifyNumu(std::string input_file_path, std::string truth_file_path);

		// Destructor.
		~ClassifyNumu();

		// Setter & Getter
		void SetOutputFilePath(std::string path) { output_file_path_ = path; }
		void SetEventIDPrefix(std::string prefix) { event_id_prefix_ = prefix; }
		void SetCut(TCut cut) { cut_ = cut; }

		// Run ClassifyNumu.
		void Run();

		// Classify specific event.
		void Classify(std::string, int event_id);

		void ClearPVR();

		bool DoesMuonExist(int event_id);

		double MuonEnergy(int event_id);
		double NumberOfMuon(int event_id);
		bool DoesMuonOutgo(int event_id);
		int MuonLastPlate(int event_id);

		
};

// --------------------

ClassifyNumu::ClassifyNumu(std::string input_file_path, std::string truth_file_path) : TruthManager(truth_file_path), input_file_path_(input_file_path), cut_("") {

	std::cout << "Constructor of ClassifyNumu called." << std::endl;

	std::ifstream ifs;
	ifs.open(input_file_path_);

	if (!ifs) {
		std::cerr << "Invalid file! Check the input file." << std::endl;
		exit(1);
	}

	std::string line_buf;
	while (std::getline(ifs, line_buf)) {
		input_files_.push_back(line_buf);
	}

	std::cout << "Read " << input_files_.size() << " files." << std::endl;

	dproc_ = new EdbDataProc;
	pvr_ = new EdbPVRec;
}

// --------------------

void ClassifyNumu::Run() {

	if (output_file_path_.empty()) {
		std::cerr << "Error! Set OutputFilePath." << std::endl;
		exit(1);
	}

	if (event_id_prefix_.empty()) {
		std::cerr << "Error! Event ID prefix does not exist." << std::endl;
		exit(1);
	}


	ofs_.open(output_file_path_);
	ofs_ << "Event ID, IsNuMuCC, Energy, Number of muon, Is outgo, Last Plate" << std::endl;

	for (auto iter : input_files_) {
		int event_id = 0;
		sscanf(iter.c_str(), event_id_prefix_.c_str(), &event_id);
		event_id += 100000;

		std::string path = iter + "/linked_tracks.root";

		Classify(path, event_id);
	}

}


// --------------------

void ClassifyNumu::Classify(std::string path, int event_id) {

	bool is_muon_exist;
	double energy;
	int num_mu_duplicate;
	bool is_outgo;
	int last_plate;

	dproc_ -> ReadTracksTree(*pvr_, path.c_str(), cut_);

	is_muon_exist = DoesMuonExist(event_id);
	energy = MuonEnergy(event_id);
	num_mu_duplicate = NumberOfMuon(event_id);
	is_outgo = DoesMuonOutgo(event_id);
	last_plate = MuonLastPlate(event_id);
	

	ofs_ << event_id;
	ofs_ << ", " << (is_muon_exist ? 1 : 0);
	ofs_ << ", " <<  energy;
	ofs_ << ", " << num_mu_duplicate;
	ofs_ << ", " << (is_outgo ? 1 : 0);
	ofs_ << ", " << last_plate;

	ofs_ << std::endl;

	ClearPVR();
}

// --------------------

void ClassifyNumu::ClearPVR() {
	if (pvr_ -> eTracks) pvr_ -> eTracks -> Clear();
}

// --------------------

bool ClassifyNumu::DoesMuonExist(int event_id) {

	for (int i=0; i<pvr_->Ntracks(); i++) {
		EdbTrackP* track = pvr_ -> GetTrack(i);

		if (!this -> IsTrack(track)) continue;
		if (track -> GetSegmentFirst() -> MCEvt() == event_id and abs(track -> GetSegmentFirst() -> MCTrack()) == 13) return true;
	}

	return false;
}

// --------------------

double ClassifyNumu::MuonEnergy(int event_id) {
	
	for (int i=0; i<pvr_->Ntracks(); i++) {
		EdbTrackP* track = pvr_ -> GetTrack(i);

		if (!this -> IsTrack(track)) continue;

		if (track -> GetSegmentFirst() -> MCEvt() == event_id and abs(track -> GetSegmentFirst() -> MCTrack()) == 13)
			return track -> GetSegmentFirst() -> P();
	}

	return -999;
}

// --------------------

double ClassifyNumu::NumberOfMuon(int event_id) {

	int num_mu_duplicate = 0;

	for (int i=0; i<pvr_->Ntracks(); i++) {
		EdbTrackP* track = pvr_ -> GetTrack(i);

		if (!this -> IsTrack(track) or track -> GetSegmentFirst() -> MCEvt() != event_id) continue;

		int pdg_middle_seg = track -> GetSegment(track->N()/2) -> MCTrack();
		if (abs(pdg_middle_seg) == 13) num_mu_duplicate++;
	}

	return num_mu_duplicate;
}

// --------------------

bool ClassifyNumu::DoesMuonOutgo(int event_id) {
	for (int i=0; i<pvr_->Ntracks(); i++) {
		EdbTrackP* track = pvr_ -> GetTrack(i);

		if (!this -> IsTrack(track) or track -> GetSegmentFirst() -> MCEvt() != event_id) continue;

		int pdg_middle_seg = track -> GetSegment(track->N()/2) -> MCTrack();
		if (abs(pdg_middle_seg) == 13) return Utils::IsOutgo(track, pvr_);
	}

	return false;
}

// --------------------

int ClassifyNumu::MuonLastPlate(int event_id) {
	std::vector<EdbTrackP*> tracks;

	for (int i=0; i<pvr_->Ntracks(); i++) {
		EdbTrackP* track = pvr_ -> GetTrack(i);

		if (!this -> IsTrack(track) or track -> GetSegmentFirst() -> MCEvt() != event_id) continue;

		int pdg_middle_seg = track -> GetSegment(track->N()/2) -> MCTrack();
		if (abs(pdg_middle_seg) == 13) tracks.push_back(track);
	}

	if (tracks.size() == 0) return -999;

	if (tracks.size() > 1) {
		std::sort(tracks.begin(), tracks.end(), 
				[] (EdbTrackP* lhs, EdbTrackP* rhs) { return lhs -> GetSegmentFirst() -> PID() < rhs -> GetSegmentFirst() -> PID(); }
			);
	}
	int last_plate = tracks.front() -> GetSegmentLast() -> ScanID().GetPlate();

	return last_plate;
}

// --------------------

int main() {

	ClassifyNumu* cnm = new ClassifyNumu("./input_list_2.txt", "./truth_input.txt");

	cnm -> SetOutputFilePath("./output_classify.csv");
	cnm -> SetEventIDPrefix("20230107_nuall/evt_%d");
	cnm -> SetCut("nseg>3");
	cnm -> Run();

}
