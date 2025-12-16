#ifndef GLOBAL_H_
#define GLOBAL_H_
#include <cmath>
#include <algorithm>
using std::min;
using std::max;
using std::swap;
using std::sort;
#include <iostream>
#include <vector>
#include <string>
using std::ostream;
#include "MersenneTwister.h"
#include <assert.h>
#include "MyVectors.h"
#include "KDtree.h"
#include "XPCPSample.h"
#include "Eigen/Core"
#include "MyAdvColor.h"
#include "QGradientStops"
#include "myProgDef.h"
#include "Params.h"
//#include "VolRenParams.h"
#include "RenderingParameters.h"
//#include "Eigen/SVD"
#ifdef M_PI
#undef M_PI
#endif
#define M_PI    3.14159265358979323846f
#define INV_PI  0.31830988618379067154f
#define INV_TWOPI  0.15915494309189533577f
#ifndef INFINITY
#define INFINITY FLT_MAX
#endif
#define PBRT_VERSION 1.0
#define RAY_EPSILON 1e-3f
#define BUF_OVERLAP_ID 100
#define MAX_CLUSTER_NUM 8
#define MAX_VOLS 10
typedef unsigned int u_int;
// Global forward declarations
class Vector;
class Point;
class Normal;
class Ray;
class RayDifferential;
class BBox;
class Transform;
class ANNkd_tree;
struct Sample;
class pFlatPt;
class PCPInselberg;

inline float Lerp(float t, float v1, float v2) {
	return (1.f - t) * v1 + t * v2;
}
inline float Clamp(float val, float low, float high) {
	if (val < low) return low;
	else if (val > high) return high;
	else return val;
}
inline int Clamp(int val, int low, int high) {
	if (val < low) return low;
	else if (val > high) return high;
	else return val;
}
inline int Mod(int a, int b) {
    int n = int(a/b);
    a -= n*b;
    if (a < 0)
        a += b;
    return a;
}
inline float Radians(float deg) {
	return ((float)M_PI/180.f) * deg;
}
inline float Degrees(float rad) {
	return (180.f/(float)M_PI) * rad;
}
inline float Log2(float x) {
	static float invLog2 = 1.f / logf(2.f);
	return logf(x) * invLog2;
}
inline int Log2Int(float v) {
	return ((*(int *) &v) >> 23) - 127;
}
inline bool IsPowerOf2(int v) {
	return (v & (v - 1)) == 0;
}
inline unsigned int RoundUpPow2(unsigned int v) {
	v--;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	return v+1;
}


// Test if a 2D point is inside a polygon
//http://alienryderflex.com/polygon/
//  Globals which should be set before calling this function:
//
//  int    polyCorners  =  how many corners the polygon has (no repeats)
//  float  polyX[]      =  horizontal coordinates of corners
//  float  polyY[]      =  vertical coordinates of corners
//  float  x, y         =  point to be tested
//
//  (Globals are used in this example for purposes of speed.  Change as
//  desired.)
//
//  The function will return YES if the point x,y is inside the polygon, or
//  NO if it is not.  If the point is exactly on the edge of the polygon,
//  then the function may return YES or NO.
//
//  Note that division by zero is avoided because the division is protected
//  by the "if" clause which surrounds it.

inline bool pointInPolygon(float x, float y, const std::vector<FLOATVECTOR2>& poly)
{
	int polyCorners = poly.size();
	int   i, j = polyCorners - 1;
	bool  oddNodes = false;

	for (i = 0; i<polyCorners; i++) {
		if ((poly[i].y < y && poly[j].y >= y
			|| poly[j].y < y && poly[i].y >= y)
			&& (poly[i].x <= x || poly[j].x <= x)) {
			if (poly[i].x + (y - poly[i].y) / (poly[j].y - poly[i].y)*(poly[j].x - poly[i].x)<x) {
				oddNodes = !oddNodes;
			}
		}
		j = i;
	}

	return oddNodes;
}

// Global variables
extern MTRand g_randGen;
static bool   g_isGlobalInit = false;
// init global variables
extern void  global_init();
// cleanup global variables
extern void  global_cleanup();
extern float RandomFloat();
extern u_int RandomUInt() ;
extern void  RandomNormDistr(std::vector<float>& y, int num, float m, float s); // random number follows Gaussian distribution with  

extern SCVecFieldParam2D g_scVec2DParam;
extern float g_samplePortion;
extern int   g_sampleNeighborSize; // number of neighbors [in the Spatial Domain] for regression computation
extern FLOATVECTOR3 g_spatialSampleNeighborSize; // radii of the axis aligned ellipse of the spatial neihborhood; 
//By default, the sampling region is a sphere with the radius of g_sampleNeighborSize

extern int   g_numNN_valDomain; // number of nearest neighbors [in the Value Domain] for regression computation
extern int   g_pcpDispBufWidth;  // width of the display buffer of the PCP plot
extern int   g_pcpDispBufHeight; // height of the display buffer of the PCP plot
extern int   g_splomDispBufWidth; // width of the display buffer of the SPLOM
extern int   g_splomDispBufHeight; // height of the display buffer of the SPLOM

// A global variable that stores the whole eXtended PCP data (RAW values + 1-flat + 2-flat)!!
extern std::vector<XPCPSample>         g_xpcpData; 
// Stores the major eigen vector of the neighborhood fitting
extern std::vector<Eigen::VectorXf>    g_majEigData; 
// Have a global variable of Inselberg's PCP is handy for coordinates conversion!
extern PCPInselberg*                   g_pcp;
// Separate vector storage for xtended pcp data
extern std::vector<std::vector<float>> g_sampled_raw;

// unique raw data
extern std::vector<std::vector<float>> g_unique_sampl_raw; // the unique sampled raw data! 
extern std::vector<FLOATVECTOR2>       g_attrMinMax; // min/max of every attribute

// It would be easier to have the KDtrees as global variables
extern std::vector<zeroFlatPt*> g_0flat_list;
extern KD<zeroFlatPt*>*       g_kdTree_rawData;
// pFlatPt representation for p-flats
extern std::vector<pFlatPt*>  g_1flat_list;
extern std::vector<pFlatPt*>  g_2flat_list;
extern KD<pFlatPt*>*          g_kdTree_1flat;
extern KD<pFlatPt*>*          g_kdTree_2flat;
// Another set of indexed points for comparison
extern std::vector<XPCPSample>         g_xpcpDataS2;
extern std::vector<pFlatPt*>  g_1flat_listS2;
extern std::vector<pFlatPt*>  g_2flat_listS2;
// Query result
extern std::vector< std::vector<UINT64>>    g_xpcp_query_result;	// the polyline layer query result if we don't have KD-tree
extern std::vector< std::vector<zeroFlatPt*> > g_0flat_query_result; // an alternative of polyline layer if we have KD-tree
extern std::vector< std::vector<pFlatPt*> > g_1flat_query_result;
extern std::vector< std::vector<pFlatPt*> > g_2flat_query_result;
// Color associated with query brushes
extern std::vector< FLOATVECTOR3 > g_xpcp_brushColors;
extern std::vector< FLOATVECTOR3 > g_1flat_brushColors;
extern std::vector< FLOATVECTOR3 > g_2flat_brushColors;
// colors for rendering
extern INTVECTOR4 g_pcp_basic_layer_color;
extern INTVECTOR4 g_pcp_1flat_layer_color;
extern INTVECTOR4 g_pcp_2flat_layer_color;
extern INTVECTOR4 g_sc_basic_layer_color;
extern INTVECTOR4 g_sc_1flat_layer_color;
extern std::vector<INTVECTOR4> g_pcp_1flat_colormap;
extern std::vector<INTVECTOR4> g_pcp_2flat_colormap;
extern std::vector<INTVECTOR4> g_clusterColors;

extern UINT64     g_min_rec_num_thres_for_sample;
extern int        g_sc_pen_width;
extern int        g_pcp_pen_width;
extern float      g_pcp_pen_alpha;
extern float      g_sc_pen_alpha;
extern double     g_pcp_render_scale; // the parameter controlling accumulation to color mapping for pcp rendering 
extern int        g_num_sample_per_update;
extern double     g_splom_render_scale; 
extern double     g_pFlatIdPt_thres; // Threshold for the p-flat indexed point strength
extern int        g_trendData_stride;
extern FLOATVECTOR3 g_pcp_layer_weights; // weights for blending layers of the PCP
extern bool       g_use_repeat_pcp; // Are we using repeating pcp (by default, no)
// Display area in Inselberg's coordinates for parallel coordinates
extern FLOATVECTOR2 g_minPCInselbergPt;
extern FLOATVECTOR2 g_maxPCInselbergPt;
// Attribute names
extern std::vector<std::string> g_attrib_names; // names of all attributes
// 
extern std::vector<FLOATVECTOR3> g_brushColormap; // a predefined colormap for brushes

extern bool       g_antialiasLines;  // draw anti-aliased lines???

extern int        g_superSampleNum; // sample number for supersampling

extern bool       g_processFullSpace; // Are we using the full space to find NN & eigenvectors ???
extern int        g_uncertain_samples_per_distr; // number of samplings to do to get different neighbors in a neighborhood.

extern DAT_VOL_INFO  g_datVolInfo; // the variable that records the dat vector file
extern bool          g_isVecData; // if it's a vector data
//======================================================
// Tranformation related
extern bool          g_isScaleXform; // use scaling function in transformation? !!! DO NOT turn off other than to show the unscaled effects.

extern NN_QUERY_METHOD g_nnQueryMethod; // the method for neighborhood query
// PCP render modes
extern int            g_pcp_num_samples_skip_DS_LINE_MODE; // how many samples we skip when drawing the simplified pcp
extern int            g_pcp_num_samples_skip_LINE_IDXPT_MODE; // how many samples to skip when drawing line+idx pt mode?
extern double         g_pcp_downsample_max_log_val;
// colormap variables
extern std::vector<MyAdvColor>     g_cmColormap; // the current colormap as an global parameter
extern std::vector<MyAdvColorNode> g_cmCtrlPoints; // control points for the colormap

extern std::vector<QGradientStops>  g_cm_pcp_tfs; // base layer colormaps  
extern std::vector<QGradientStops>  g_alpha_pcp_tfs; // base layer alpha maps
extern int                         g_cmColormapLen;
extern COLOR_SCHEME g_visColorScheme;
extern int          g_num_segCurve_halfdist;
// hdr variable
extern bool            g_isHDR;
extern double        g_boundScale; // the scaling factor at bounds
extern Params        g_params;
//extern VolRenParams  g_ext_rendering_parameters;
extern RenderingParameters  g_ext_rendering_parameters;

// the opengl colormap 1d texture
extern unsigned int g_texCm1Did;
extern unsigned int g_texCM1DXform;
extern unsigned int g_texSFCorder;
extern bool g_isOpenglInit; // if opengl has been initialized?
// helper functions
const int MAX_OCTREE_LEVELS = 8;
const int MAX_NUM_VOLS = 6;
const int MAX_NUM_BRUSHES = 20;
// evaluate spline functions

inline double poly_func(const vector<double>& p, double x)
{
	double y = 0.0;
	int deg = p.size() - 1;
	for (size_t i = 0; i < p.size(); i++)
	{
		y += p[i] * pow(x, deg);
		deg--;
	}
	return y;

}
inline double spline_interp(double x, const std::vector<std::vector<double>>& coefs, const std::vector<double>& breaks)
{
	if (x < breaks[0] )
	{
		return poly_func(coefs[0], x - breaks[0]);
	}
	else if (x >= breaks.back())
	{
		return poly_func(coefs.back(), x - breaks.back());
	}
	
	// find the corresponding ctrl points for x
	vector<double>::const_iterator it = lower_bound(breaks.begin(), breaks.end(), x);
	int cid = it - 1 - breaks.begin();
	double x1 = breaks[cid];
	double y = poly_func(coefs[cid], x - x1);
	return y;
}



#endif