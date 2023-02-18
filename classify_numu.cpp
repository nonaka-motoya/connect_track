#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <ostream>

#include <EdbDataSet.h>

#include "EdbDataSet.h"
#include "TruthManager.hpp"

class ClassifyNumu : private TruthManager {
	private:
		std::string input_file_path_;
		std::vector<std::string> input_files_;
		std::string event_id_prefix_;

		std::string output_file_path_;
		std::ofstream ofs_;

		EdbDataProc* dproc_;
		EdbPVRec* pvr_;

	public:
		// Constructor.
		ClassifyNumu(std::string input_file_path, std::string truth_file_path);

		// Destructor.
		~ClassifyNumu();

		// Setter & Getter
		void SetOutputFilePath(std::string path) { output_file_path_ = path; }
		void SetEventIDPrefix(std::string prefix) { event_id_prefix_ = prefix; }

		// Run ClassifyNumu.
		void Run();

		// Classify specific event.
		void Classify(std::string, int event_id);

		void ClearPVR();

		
};

// --------------------

ClassifyNumu::ClassifyNumu(std::string input_file_path, std::string truth_file_path) : TruthManager(truth_file_path), input_file_path_(input_file_path) {

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

	for (auto iter : input_files_) {
		int event_id = 0;
		sscanf(iter.c_str(), event_id_prefix_.c_str(), &event_id);

		std::string path = iter + "/linked_tracks.root";

		Classify(path, event_id);
	}

}


// --------------------

void ClassifyNumu::Classify(std::string path, int event_id) {

	dproc_ -> ReadTracksTree(*pvr_, path.c_str(), "");

	ofs_ << event_id << ", " << pvr_ -> Ntracks() << std::endl;

	ClearPVR();
}

// --------------------

void ClassifyNumu::ClearPVR() {
	if (pvr_ -> eTracks) pvr_ -> eTracks -> Clear();
}

// --------------------

int main() {

	ClassifyNumu* cmn = new ClassifyNumu("./input_list_2.txt", "./truth_input.txt");

	cmn -> SetOutputFilePath("./output_classify.csv");
	cmn -> SetEventIDPrefix("20230107_nuall/evt_%d");
	cmn -> Run();

}
