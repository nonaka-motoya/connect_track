#include <algorithm>
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <string>
#include <set>

#include "EdbEDAUtil.h"
#include "EdbPVRec.h"
#include "Rtypes.h"
#include "TH1.h"
#include "TObjArray.h"
#include "TRint.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TString.h"
#include "TSystem.h"

#include "EdbDataSet.h"

#include "Utils.hpp"
#include "HashTable.hpp"

int main() {

	TRint app("app", 0, 0);


	std::ifstream input_file;
	input_file.open("input_list_2.txt", std::ios::in);
	std::string line_buf;
	int num_file = 0;

	// development
	std::string path = "20230107_nuall/evt_2950_pl1_300/linked_tracks.root"; // path of linked_tracks.
	//std::string path = "20230107_nuall/evt_4812_pl1_300/linked_tracks.root"; // path of linked_tracks.
	Utils::ConnectTrack(path);
	std::cout << "Read " << path << std::endl;

	app.Run();

	return 0;
}
