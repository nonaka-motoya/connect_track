#ifndef UTILS_H
#define UTILS_H

#include <EdbDataSet.h>

class Utils {

	public:
		static int CheckTrackingQuality(EdbTrackP* track);
		static void PrintTrack(EdbTrackP* track);
		static int SearchVertex(EdbTrackP* track, int index, EdbPVRec* pvr);
		static bool IsSwitch(EdbTrackP* track, int index);

		// calculate track angle.
		// use segments in [i-1, i+1].
		// return pair of tangent.
		static std::pair<double, double> CalcTrackAngle(EdbTrackP* track, int index);

		// calculate track angle difference (mrad).
		// fit x-z (y-z) space using segments in [i-3, i-1] and [i, i+2].
		static double CalcTrackAngleDiff(EdbTrackP* track, int index);

		static EdbPVRec *ConnectTrack(std::string path, double distance, double angle);
		static std::pair<double, double> Theta(EdbTrackP *track);
		static double Dtheta(EdbTrackP *track1, EdbTrackP *track2);
		static EdbSegP *MakeVirtualSegment(EdbTrackP *track);
		static double Distance(EdbTrackP *track1, EdbTrackP *track2);
		static bool IsOutgo(EdbTrackP* track, EdbPVRec* pvr);
		static bool HasKink(EdbTrackP* track);
};


#endif
