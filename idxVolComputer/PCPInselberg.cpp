#include "PCPInselberg.h"

using namespace std;

PCPInselberg::PCPInselberg(const std::vector<std::string>& attributes, bool repeat):
ParallelCoordinates(attributes),
_repeat(repeat)
{
	if (_repeat)
		_d_max = float(attributes.size() * 2 - 1); // the distance between the last axis and the first axis
	else
		_d_max = float(attributes.size() - 1);
	_num_attribs = int(attributes.size());
}


PCPInselberg::~PCPInselberg()
{
}

void PCPInselberg::cartesian2PCP(const Point& pt, float di, Point& vp1, Point& vp2)
{
	vp1.x = di;
	vp1.y = pt.x;

	vp2.x = di + _d;
	vp2.y = pt.y;
}

void PCPInselberg::cartesian2PCP(const vector<float>& cart_val, vector<Point>& pcp_pts)
{
	pcp_pts.resize(cart_val.size());
	for (size_t i = 0; i < cart_val.size(); i++)
	{
		pcp_pts[i].x = float(i);
		pcp_pts[i].y = _axesYrange * cart_val[i] + _axesYbase;
	}
}

void PCPInselberg::cartesian2PCP_repeatAxes(const vector<float>& cart_val, vector<Point>& pcp_repeat_pts)
{
	pcp_repeat_pts.resize(cart_val.size() * 2);
	for (size_t i = 0; i < cart_val.size(); i++)
	{
		pcp_repeat_pts[i].x = float(i);
		pcp_repeat_pts[i].y = _axesYrange * cart_val[i] + _axesYbase;

		pcp_repeat_pts[cart_val.size() + i].x = float(cart_val.size() + i); 
		pcp_repeat_pts[cart_val.size() + i].y = _axesYrange * cart_val[i] + _axesYbase;
	}
}

void PCPInselberg::cartesian2PCP_subspace(const vector<float>& full_cart_val, int startSubspaceId, int endSubspaceId, vector<Point>& subspace_pcp_pts)
{
	int subspaceDim = endSubspaceId - startSubspaceId;
	subspace_pcp_pts.resize(subspaceDim);
	for (size_t i = startSubspaceId; i < endSubspaceId; i++)
	{
		subspace_pcp_pts[i - startSubspaceId].x = float(i);
		subspace_pcp_pts[i - startSubspaceId].y = _axesYrange * full_cart_val[i] + _axesYbase;
	}
}

// Convert multivariate point to lines (in order) with repeating axes for subspace values. 
void PCPInselberg::cartesian2PCP_subspace_repeatAxes(const vector<float>& full_cart_val, int startSubspaceId, int endSubspaceId, vector<Point>& subspace_pcp_repeat_pts)
{
	int subspaceDim = endSubspaceId - startSubspaceId;
	subspace_pcp_repeat_pts.resize(subspaceDim * 2);
	for (size_t i = startSubspaceId; i < endSubspaceId; i++)
	{
		subspace_pcp_repeat_pts[i - startSubspaceId].x = float(i);
		subspace_pcp_repeat_pts[i - startSubspaceId].y = _axesYrange * full_cart_val[i] + _axesYbase;

		subspace_pcp_repeat_pts[subspaceDim + i - startSubspaceId].x = float(full_cart_val.size() + i);
		subspace_pcp_repeat_pts[subspaceDim + i - startSubspaceId].y = _axesYrange * full_cart_val[i] + _axesYbase;
	}
}


int PCPInselberg::cartesian2PCP(const Point& pt1, const Point& pt2, float di, Point& vp)
{
	// point to line
	Point vp1, vp2, vp3, vp4;
	cartesian2PCP(pt1, di, vp1, vp2);
	cartesian2PCP(pt2, di, vp3, vp4);

	Segment S0, S1;
	S0.P0 = vp1;
	S0.P1 = vp2;

	S1.P0 = vp3;
	S1.P1 = vp4;

	Point ep;
	// compute the intersection of two segments
	return  intersect2D_2Segments(S0, S1, vp, ep);
}



bool PCPInselberg::calcFlatOneOrderHigher(const vector<Point>& pcp_pm1flat_1, const vector<Point>& pcp_pm1flat_2, std::vector<Point>& pcp_flat_p)
{
	pcp_flat_p.resize(pcp_pm1flat_1.size() - 1);
	int k = 1;
	for (size_t i = 0; i < pcp_pm1flat_1.size() - 1; i++)
	{
		// compute the intersections formed by line segements connecting neighboring points of p-1 flat 1 and p-1 flat 2

		// First p-1 flat 
		Point p11 = pcp_pm1flat_1[i];
		Point p12 = pcp_pm1flat_1[i + 1];
		// Second p-1 flat
		Point p21 = pcp_pm1flat_2[i];
		Point p22 = pcp_pm1flat_2[i + 1];

		Segment s1;
		s1.P0 = p11;
		s1.P1 = p12;

		Segment s2;
		s2.P0 = p21;
		s2.P1 = p22;

		// index point for the new p-flat
		Point pcp_flat_p_pt;
		Point dummy;
		int inter = intersect2D_2Segments(s1, s2, pcp_flat_p_pt, dummy);
		k = k * inter;
		if (inter == 0)
		{
			pcp_flat_p[i] = Point(INFINITY, INFINITY);
		}
		else
		{
			pcp_flat_p[i] = pcp_flat_p_pt;
		}

	}
	if (k == 0)
		return false;
	else
		return true; // A-ok!
}


void PCPInselberg::calcOneFlatGeneralForm(const vector<float>& cart_pm1flat_1, const vector<float>& cart_pm1flat_2, std::vector<Point>& pcp_oneflat)
{
	pcp_oneflat.resize(cart_pm1flat_1.size() - 1);
	for (size_t i = 0; i < cart_pm1flat_1.size() - 1; i++)
	{
		// compute the intersections formed by line segements connecting neighboring points of p-1 flat 1 and p-1 flat 2

		// First 2D cartesian point 
		Point p1(cart_pm1flat_1[i], cart_pm1flat_1[i+1]);
	
		// Second p-1 flat
		Point p2(cart_pm1flat_2[i], cart_pm1flat_2[i + 1]);

		// general form of the line
		float c1 = p2.y - p1.y;
		float c2 = -(p2.x - p1.x);
		float c3 = p1.y * (p2.x - p1.x) - p1.x * (p2.y - p1.y);
		
		pcp_oneflat[i] = Point(c1, c2, c3);
	}

}

void PCPInselberg::calcOneFlatParametricForm(const vector<float>& cart_pm1flat_1, const vector<float>& cart_pm1flat_2, std::vector<Point>& pcp_oneflat)
{
	pcp_oneflat.resize(cart_pm1flat_1.size() - 1);
	for (size_t i = 0; i < cart_pm1flat_1.size() - 1; i++)
	{
		// compute the intersections formed by line segements connecting neighboring points of p-1 flat 1 and p-1 flat 2

		// First 2D cartesian point 
		Point p1(cart_pm1flat_1[i], cart_pm1flat_1[i + 1]);

		// Second p-1 flat
		Point p2(cart_pm1flat_2[i], cart_pm1flat_2[i + 1]);

		// parametric form with normalization
		Vector v = p1 - p2;
		v = v / (v.Length());

		pcp_oneflat[i] = Point(v.x, v.y);
	}

}

void PCPInselberg::calc2FlatFrom3DPlane(const PlaneCoeffs& plane, vector<FLOATVECTOR2>& indexPts, bool oneIndexPt)
{
	indexPts.clear();
	if (oneIndexPt)
		indexPts.resize(1); // have 1 index point
	else
		indexPts.resize(4); // have 4 index points

	float oneOverDenom = 1.0f / (plane.c1 + plane.c2 + plane.c3);
	
	FLOATVECTOR2 pi;
	pi.x = (plane.c2 + 2.0f * plane.c3) * oneOverDenom;
	pi.y = plane.c0 * oneOverDenom;
	indexPts[0] = pi;
	if (!oneIndexPt)
	{ // additional 3 index points
		pi.x += 3.0f * plane.c1 * oneOverDenom;
		indexPts[1] = pi;

		pi.x += 3.0f * plane.c2 * oneOverDenom;
		indexPts[2] = pi;

		pi.x += 3.0f * plane.c3 * oneOverDenom;
		indexPts[3] = pi;
	}

}

void PCPInselberg::calc2FlatFrom3DPlane_GeneralForm(const PlaneCoeffs& plane, vector<FLOATVECTOR3>& indexPts, bool oneIndexPt)
{
	indexPts.clear();
	if (oneIndexPt)
		indexPts.resize(1); // have 1 index point
	else
		indexPts.resize(4); // have 4 index points

	float oneOverDenom = 1.0f / (plane.c1 + plane.c2 + plane.c3);

	FLOATVECTOR3 pi;
	pi.x = (plane.c2 + 2.0f * plane.c3) * oneOverDenom;
	pi.y = plane.c0 * oneOverDenom;
	pi.z = 1;
	indexPts[0] = pi;
	if (!oneIndexPt)
	{ // additional 3 index points
		pi.x += 3.0f * plane.c1 * oneOverDenom;
		indexPts[1] = pi;

		pi.x += 3.0f * plane.c2 * oneOverDenom;
		indexPts[2] = pi;

		pi.x += 3.0f * plane.c3 * oneOverDenom;
		indexPts[3] = pi;
	}

}