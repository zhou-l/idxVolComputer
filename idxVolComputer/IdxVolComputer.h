#pragma once
#include <QString>
#include <vector>
class VolumeData;
class PointDataHandler;
class XPCPdataAnalyzer;

class IdxVolComputer
{
public:
	IdxVolComputer();
	~IdxVolComputer();

	void processSpatialDomainContinuousIdxPt(const std::vector<VolumeData*>& volData);
	static bool loadVolumeFile(const QString& fileName, VolumeData** vol, float& minVal, float& maxVal, bool doNormalize = true);
private:

///////////////////////////////////////////////////
// Data processor: 
//////////////////////////////////////////////////
 //The processor should be either a VolumeSampler or a simple PointDataHandler based on input!
	PointDataHandler* _dataHandler;
	bool                _isVolHandler; // a flag indicating if it's a volume handler
	////////////////////////////////////////////////
	XPCPdataAnalyzer* _xpcpAnalyzer;// Process raw data

};

