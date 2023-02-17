#ifndef TRUTHMANAGER_H
#define TRUTHMANAGER_H

#include <iostream>
#include <vector>
#include <map>

#include <EdbDataSet.h>

class TruthManager {
	private:
		std::string truth_path_;
		std::vector<std::string> truth_file_path_;
		std::map<int, std::vector<int>> unique_id_;

	public:
		TruthManager(std::string input_path);
		~TruthManager();

		void ReadTruth(std::string path);

		bool IsTrack(EdbTrackP* track);

};


#endif
