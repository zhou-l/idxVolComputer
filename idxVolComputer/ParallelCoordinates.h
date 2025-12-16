#pragma once
#include "geometry.h"
#include <string>
#include <map>
#include "MyVectors.h"

class PCPpatch
{
public:
	float atr1_ext; // extent of attribute 1 
	float atr2_ext; // extent of attribute 2
};

class ParallelCoordinates
{
public:

	ParallelCoordinates(const std::vector<std::string>& attributes);
	virtual ~ParallelCoordinates();
	void changeAxisOrder();

	// get normalized parallel coordinates line from one normalized 2D point
	// PCP -> Scatterplot duality
	// 1. Point->Lines
	

	void cartesian2PCP(const Point& pt, const std::string& atr1, const std::string& atr2, Point& vp1, Point& vp2); // point to line (Use attribute strings)
	void cartesian2PCP(const std::vector<float>& cart_val, const std::vector<std::string>& atrOrder, std::vector<Point>& pcp_pts); // multivariate point to lines
	
	// 2. Lines -> Point duality
	// The output is in the [0,2*(_axes.size()-1)*d]!
	int cartesian2PCP(const Point& pt1, const Point& pt2, const std::string& atr1, const std::string& atr2, Point& vp, bool flipY = false); // line to point: return 0, no intersection, 1, unique intersectoin, 2, overlaps, 3, overlap outside the segments; flipY is a flag for whether to flip the Y components in pt1 & pt2
	///////////////////////////////////////////////////

	// Compute statistics (TODO: add more later)
	void calcAveragedAccumDir(int intervalId, double theta); // compute the averaged direction for given interval id
	const std::vector<double>& getAvgDirAccum() { return _avgDirList; }

	void flipPtY(Point& pt);
	const std::vector<Point>&   getAxes() const { return _axes; }
protected:
	std::map<std::string, int>  _attribs; // all attributes 
	std::vector<Point>          _axes;    // location of vertical axes (normalized)
	int                         _ndims;   // number of dimensions
	std::vector<double>         _avgDirList; // a list of the average direction between axes.
	UINT64                      _sampleCnt;  // keep track of how many samples have been processed.

	Point                       _boundTopLeft; 
	Point                       _boundBotRight; 

	float                       _axesYrange; // full range of Y axis in [0,1]
	float                       _axesYbase;  // base y position of the Y axis in [0,1]

	std::vector<PCPpatch>       _patch_list;
	float                       _d; // spacing between two axes
};

