#pragma once
#include "ParallelCoordinates.h"

// Extension to basic PCP with 1-flat and 2-flat
struct PCPbtAxesData // data item between PCP axes 
{

	std::vector< float >  x; // N-D value from the cartesian coordinates
	float i; // the axis id
	float c1, c2, c3, c0; // coefficients for 2-flat
};

// PCP using Inselberg's original cooridnates:
// |     |     |     |     |     |
// |     |     |     |     |     |
// |     |     |     |     |     |
// X1(Y) X2    X3    X1'   X2'   X3'
// d1=0  d2=1  d3=2  d1'=3 d2'=4 d3'=5
// 0                             _d_max
// The spacing between neighboring axes is 1


class PCPInselberg : public ParallelCoordinates
{
public:
	PCPInselberg(const std::vector<std::string>& attributes, bool repeat = true);
	virtual ~PCPInselberg();
	
	// Accessors
	bool isRepeat() const { return _repeat;  }
	bool toggleRepeat() { _repeat = !_repeat; return _repeat; } // toggle repeat axes and return current state
	
	int numAttribs() const { return _num_attribs; }
	// maximum x distance to the first axis
	float max_d() const { return _d_max; }
	// Dualities: ALL outputs are in [0, _d_max]
	// 1. Point->lines
	// Given cartesian point and the distance to Y axis for vp1's x coordinate (According to Inselberg's book) 
	void cartesian2PCP(const Point& pt, float di, Point& vp1, Point& vp2);
	void cartesian2PCP(const std::vector<float>& cart_val, std::vector<Point>& pcp_pts); // multivariate point to lines (no repeat)
	// multivariate point to lines (in order) with repeating axes. 
	void cartesian2PCP_repeatAxes(const std::vector<float>& cart_val, std::vector<Point>& pcp_repeat_pts);
	
	// Convert multivariate point to lines (in order) for subspace values. 
	void cartesian2PCP_subspace(const vector<float>& full_cart_val, int startSubspaceId, int endSubspaceId, vector<Point>& subspace_pcp_pts); // without repeating axes 
	void cartesian2PCP_subspace_repeatAxes(const std::vector<float>& full_cart_val, int startSubspaceId, int endSubspaceId, std::vector<Point>& subspace_pcp_repeat_pts); // with repeating axes
	
	
	// 2. Lines ->point
	int cartesian2PCP(const Point& pt1, const Point& pt2, float di, Point& vp); // Given cartesian point and the distance to Y axis for vp1's x coordinate (According to Inselberg's book) 
	

	///////////////////////////////////////
	// Construction of p-1 flat from two p-2 flats
	///////////////////////////////////////

	// NOTE: this function uses the 2D point coordinates in the xy-image plane
	bool calcFlatOneOrderHigher(const std::vector<Point>& pcp_pm1flat_1, 
									 const std::vector<Point>& pcp_pm1flat_2,
									 std::vector<Point>& pcp_pflat); // Return false if any of the subspace lines do not intersect!
	// NOTE: this function uses the 3D homogeneous coordinates!
	void calcOneFlatGeneralForm(const std::vector<float>& cart_pm1flat_1,
		                        const std::vector<float>& cart_pm1flat_2,
								std::vector<Point>& pcp_oneflat); // calculate the 1-flat indexed point from two Cartesian points(!) with the general line form; returns the 3D line coordinates of the indexed point (include all cases)!!!
	///////////////////////////////////////////////////////////
	void calcOneFlatParametricForm(const std::vector<float>& cart_pm1flat_1,
		const std::vector<float>& cart_pm1flat_2,
		std::vector<Point>& pcp_oneflat); // calculate the 1-flat indexed point from two Cartesian points(!) with the parametric line form; returns the direction of the line. 
	// 2-flats: use normal vectors to compute the indexed points
	void calc2FlatFrom3DPlane(const PlaneCoeffs& plane, std::vector<FLOATVECTOR2>& indexPts, bool oneIndexPt = true);
	// use 3D homogeneous coordinates 
	// NOTE: this is WRONG at the moment!
	void calc2FlatFrom3DPlane_GeneralForm(const PlaneCoeffs& plane, std::vector<FLOATVECTOR3>& indexPts, bool oneIndexPt = true);
private:
	float   _d_max; // the maximum spacing for the first vertical axis (Y)
	bool    _repeat; // repeating axes? 
	int     _num_attribs; // number of attributes
};

