#include "ParallelCoordinates.h"
#include "opencv.hpp"
using namespace std;

ParallelCoordinates::ParallelCoordinates(const vector<string>& attributes):
_d(1.0f)// default spacing is 1
{
	_boundTopLeft = Point(0.f, 0.f);
	_boundBotRight = Point(1.f, 1.f);
	_ndims = int(attributes.size());
	_sampleCnt = 0;
	float axisInterval = (_boundBotRight.x - _boundTopLeft.x) / float(_ndims - 1);
	for(int i = 0; i < _ndims; i++)
	{
		Point axis;
		axis.x = _boundTopLeft.x + axisInterval * float(i);
		axis.y = _boundBotRight.y;
		_axes.push_back( axis );
		_attribs[attributes[i]] = i;
	}
	_avgDirList.resize(_ndims - 1, 0.0); // Set direction to 0 by default
	_axesYrange = (_boundBotRight.y - _boundTopLeft.y);
	_axesYbase = _boundTopLeft.y;
}


ParallelCoordinates::~ParallelCoordinates()
{
}

void ParallelCoordinates::cartesian2PCP(const Point& pt, const string& atr1, const string& atr2, Point& vp1, Point& vp2)
{
	// end points in the PCP space
	vp1.x = _axes[_attribs[atr1]].x;
	vp1.y = pt.x;

	vp2.x = _axes[_attribs[atr2]].x;
	vp2.y = pt.y;
}

int ParallelCoordinates::cartesian2PCP(const Point& pt1, const Point& pt2, const string& atr1, const string& atr2, Point& vp, bool flipY)
{
	// point to line
	Point vp1, vp2, vp3, vp4;
	if (flipY)
	{
		Point flipYpt1 = pt1;
		Point flipYpt2 = pt2;
		flipPtY(flipYpt1);
		flipPtY(flipYpt2);
		cartesian2PCP(flipYpt1, atr1, atr2, vp1, vp2);
		cartesian2PCP(flipYpt2, atr1, atr2, vp3, vp4);
	}
	else
	{
		cartesian2PCP(pt1, atr1, atr2, vp1, vp2);
		cartesian2PCP(pt2, atr1, atr2, vp3, vp4);

	}


	Segment S0, S1;
	S0.P0 = vp1;
	S0.P1 = vp2;
	
	S1.P0 = vp3;
	S1.P1 = vp4;

	Point ep;
	// compute the intersection of two segments
	int intersect_res = intersect2D_2Segments(S0, S1, vp, ep);

	// If the two segments do not intersect but they are also not parallel, we still can return the intersection point (outside the axes).

	//vp.x = MAX(MIN(vp.x, _axes[_attribs[atr2]].x), _axes[_attribs[atr1]].x);
	//vp.y = MAX(MIN(vp.y, _boundBotRight.y), _boundTopLeft.y);
	return intersect_res;
}



void ParallelCoordinates::cartesian2PCP(const vector<float>& cart_val, const vector<string>& atrOrder, vector<Point>& pcp_pts)
{
	pcp_pts.clear();
	pcp_pts.resize(cart_val.size());
	for (size_t i = 0; i < cart_val.size(); i++)
	{
		int axisid = _attribs[atrOrder[i]]; // get the correct axis
		pcp_pts[axisid].x = _axes[axisid].x;
		pcp_pts[axisid].y = _axesYrange * cart_val[i] + _axesYbase;
	}
}


void ParallelCoordinates::flipPtY(Point& pt)
{
	pt.y = (1.0f - pt.y) * _axesYrange + _axesYbase;
}

void ParallelCoordinates::changeAxisOrder()
{}

void ParallelCoordinates::calcAveragedAccumDir(int intervalId, double theta)
{
	double denom = 1.0 / double(_sampleCnt + 1);
	_avgDirList[intervalId] = denom * (double(_sampleCnt) * _avgDirList[intervalId] + theta);
	_sampleCnt++;
}