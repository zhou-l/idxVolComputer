#include "global.h"

MTRand g_randGen;
SCVecFieldParam2D g_scVec2DParam;
float g_samplePortion;
int   g_sampleNeighborSize;
int   g_numNN_valDomain;
int   g_pcpDispBufWidth; 
int   g_pcpDispBufHeight;
int   g_splomDispBufWidth;
int   g_splomDispBufHeight;
UINT64 g_min_rec_num_thres_for_sample;

FLOATVECTOR3 g_spatialSampleNeighborSize;

KD<zeroFlatPt*>*  g_kdTree_rawData;
KD<pFlatPt*>* g_kdTree_1flat;
KD<pFlatPt*>* g_kdTree_2flat;
PCPInselberg* g_pcp;

std::vector<std::vector<float>> g_sampled_raw;
std::vector<std::vector<float>> g_unique_sampl_raw;
std::vector<FLOATVECTOR2>       g_attrMinMax; 
std::vector<XPCPSample>         g_xpcpData;
std::vector<Eigen::VectorXf>    g_majEigData;
std::vector<Eigen::VectorXf>    g_trendLines;

std::vector<zeroFlatPt*> g_0flat_list;
std::vector<pFlatPt*> g_1flat_list;
std::vector<pFlatPt*> g_2flat_list;

std::vector<pFlatPt*> g_contIdxPt_list;
// indexed points data set 2 for comparison 
std::vector<XPCPSample>         g_xpcpDataS2;
std::vector<zeroFlatPt*> g_0flat_listS2;
std::vector<pFlatPt*>	 g_1flat_listS2;
std::vector<pFlatPt*>    g_2flat_listS2;

// Query results
std::vector< std::vector<UINT64>>  g_xpcp_query_result;
std::vector<std::vector<zeroFlatPt*>> g_0flat_query_result;
std::vector<std::vector<pFlatPt*>> g_1flat_query_result;
std::vector<std::vector<pFlatPt*>> g_2flat_query_result;
std::vector<FLOATVECTOR3> g_xpcp_brushColors;
std::vector<FLOATVECTOR3> g_1flat_brushColors;
std::vector<FLOATVECTOR3> g_2flat_brushColors;

INTVECTOR4 g_pcp_basic_layer_color;
INTVECTOR4 g_pcp_1flat_layer_color;
INTVECTOR4 g_pcp_2flat_layer_color;
std::vector<INTVECTOR4> g_pcp_1flat_colormap;
std::vector<INTVECTOR4> g_pcp_2flat_colormap;
std::vector<INTVECTOR4> g_clusterColors;
FLOATVECTOR3 g_pcp_layer_weights;

vector<string> g_attrib_names;
bool       g_use_repeat_pcp;

FLOATVECTOR2 g_minPCInselbergPt;
FLOATVECTOR2 g_maxPCInselbergPt;

int        g_pcp_pen_width;
float      g_pcp_pen_alpha;
float      g_sc_pen_alpha;
double     g_pcp_render_scale;
double     g_splom_render_scale;
double     g_pFlatIdPt_thres; // Threshold for the p-flat indexed point strength

int        g_num_sample_per_update;

INTVECTOR4 g_sc_basic_layer_color;
INTVECTOR4 g_sc_1flat_layer_color;
int        g_sc_pen_width;
int        g_trendData_stride;

std::vector<FLOATVECTOR3> g_brushColormap;

bool                g_antialiasLines;
int                 g_superSampleNum; 

bool                g_processFullSpace; 
DAT_VOL_INFO        g_datVolInfo;
bool                g_isVecData;
bool                g_isScaleXform;
NN_QUERY_METHOD     g_nnQueryMethod;
PC_RENDER_MODE      g_pcp_RenderMode;
int                 g_uncertain_samples_per_distr; 
int                 g_pcp_num_samples_skip_LINE_IDXPT_MODE;
int                 g_pcp_num_samples_skip_DS_LINE_MODE;
double              g_pcp_downsample_max_log_val;

//std::vector<MyAdvColor>     g_cmColormap; // the current colormap as an global parameter
//std::vector<MyAdvColorNode> g_cmCtrlPoints; // control points for the colormap
//std::vector<QGradientStops> g_cm_pcp_tfs;// base layer COLOR 1D transfer functions 
//std::vector<QGradientStops> g_alpha_pcp_tfs; // base layer ALPHA 1D transfer functions
bool                 g_isHDR;
int                  g_cmColormapLen;
int                  g_num_segCurve_halfdist; // number of segments of half the distance between two neighoring axes
COLOR_SCHEME         g_visColorScheme;
double               g_boundScale; // the scaling factor at bounds
Params               g_params;
//VolRenParams         g_ext_rendering_parameters;
//RenderingParameters  g_ext_rendering_parameters;

unsigned int g_texCm1Did;
unsigned int g_texCM1DXform;
unsigned int g_texSFCorder;
bool g_isOpenglInit; // if opengl has been initialized?

void global_init()
{
	g_randGen.seed();
	g_isGlobalInit = true;
	// Setup param for the vector field on 2D scatterplot
	g_scVec2DParam.num_seeds = 1000;
	g_scVec2DParam.max_pts_per_curve = 100;
	g_scVec2DParam.rej_sampl_M = 50.0; 
	g_scVec2DParam.max_zero_velo_pts = 2;
	g_scVec2DParam.streamline_stepsize = 1.0f;
	g_samplePortion = 0.002f;//0.002f;//0.005f;
	g_sampleNeighborSize = 5;
	g_spatialSampleNeighborSize = FLOATVECTOR3(g_sampleNeighborSize, g_sampleNeighborSize, g_sampleNeighborSize);
	g_numNN_valDomain = 1000;// 1000;
	// KD trees
	g_kdTree_rawData = NULL;
	g_kdTree_1flat = NULL;
	g_kdTree_2flat = NULL;
	// pcp 
	g_pcp = NULL;
	// size of plots
	//g_pcpDispBufWidth = 4096;
	//g_pcpDispBufHeight = 4096;

	g_params.setPCPDispBufSize(4096, 4096);

	g_splomDispBufWidth = 2048;
	g_splomDispBufHeight = 2048;
	// rendering colors: use colorbrewer 3-class qualitative scheme
	g_pcp_basic_layer_color = INTVECTOR4(0x8d, 0xa0, 0xcb, 60);
	g_pcp_1flat_layer_color = INTVECTOR4(0x1f, 0x78, 0xb4, 60);
	g_pcp_2flat_layer_color = INTVECTOR4(0x66, 0xc2, 0xa5, 60);
	g_sc_basic_layer_color = INTVECTOR4(0, 0, 255, 100);
	g_sc_1flat_layer_color = INTVECTOR4(0xe4, 0x1a, 0xdf, 60);
	// width of the scatterplot drawing pen
	int num_pts_short_edge = 500;
	g_sc_pen_width = int(1.0 / double(num_pts_short_edge) * double(MIN(g_splomDispBufWidth, g_splomDispBufHeight)));
	g_pcp_pen_width = /*2 * */int(1.0 / double(num_pts_short_edge) * double(g_params.getPCPDispBufSize().minVal())); // Should config point size in the config menu.
	// alpha of pcp pen
	g_pcp_pen_alpha = 1.0f;
	g_sc_pen_alpha = 1.0f;
	// pcp render mapping param
	g_pcp_render_scale = 1.0;
	g_splom_render_scale = 1.0;
	// Default to not selecting any layer!
	g_pcp_layer_weights = FLOATVECTOR3(1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 3.0f);

	// number of samples to process per update
	g_num_sample_per_update = 10000;
	// Controls draw one trend for how many samples
	g_trendData_stride = 20;
	// Set pcp repeat mode
	g_use_repeat_pcp = false;
	// Set use full space
	g_processFullSpace = true; 

	// Threshold to determine at least how many records we have in data so that we need to sample only part of the data
	g_min_rec_num_thres_for_sample = 40000;

	// HACK: set pre-defined colormap for brushes
	//g_brushColormap.push_back(FLOATVECTOR3(0, 0, 1));
	//g_brushColormap.push_back(FLOATVECTOR3(1, 0, 0));
	//g_brushColormap.push_back(FLOATVECTOR3(0, 1, 0));
	//g_brushColormap.push_back(FLOATVECTOR3(0, 1, 1));
	//g_brushColormap.push_back(FLOATVECTOR3(1, 0, 1));
	//g_brushColormap.push_back(FLOATVECTOR3(1, 1, 0));

	// set threshold to 0 (show all p-flats)
	g_pFlatIdPt_thres = 0.0;

	// Set display area for Inselberg's coordinates
	//g_minPCInselbergPt = FLOATVECTOR2(-0.1f, -0.5f);
	//g_maxPCInselbergPt = FLOATVECTOR2(1.1f, 1.5f);
	// Use the config below for 2 attributes
	//g_minPCInselbergPt = FLOATVECTOR2(-0.15f, -1.25f);
	//g_maxPCInselbergPt = FLOATVECTOR2(1.15f, 1.25f);

	// No vector data by default
	g_isVecData = false;
	// use transform for idx points and pcp polylines
	//g_params.XformMethod() = XF_NONE;//XF_IDX_ONLY;

	// use scaling in transformation? ALWAYS true unless you want to see the unscaled effects!
	g_isScaleXform = true; 

	//g_params.XformMethod() = XF_IDX_ONLY;
	// set tehe nearest neighbor query methos
	g_nnQueryMethod = NNQ_KNN;
	//g_nnQueryMethod = NNQ_RBALL;
	// How many neighborhood configurations needed for the uncertian point data? 
	g_uncertain_samples_per_distr = 200;
	// Draw anti-aliasing lines? By default: no
	g_antialiasLines = false;
	// To generate figures, use supersampling

	g_superSampleNum = 9; // number of samples per sampling point

	// number of samples to skip
	g_pcp_num_samples_skip_DS_LINE_MODE = 30;
	g_pcp_num_samples_skip_LINE_IDXPT_MODE = 1;
	//g_pcp_num_samples_skip_LINE_IDXPT_MODE = g_pcp_num_samples_skip_DS_LINE_MODE;// 10;

	g_pcp_downsample_max_log_val = 1.0;
	// Set size of density plots
	g_isHDR = false;
	g_cmColormapLen = (g_isHDR) ? 1024 : 256;
	// Set color scheme of PCP visualization
	g_visColorScheme = CS_LIGHT;
	// number of segments of a curve within half the distance of 2 neighboring axes
	g_num_segCurve_halfdist = 40;
	// scaling factors at bounds
	g_boundScale = 1.0;

	g_isOpenglInit = false;
	g_texCm1Did = 0;
	g_texCM1DXform = 0;
	g_texSFCorder = 0;
}

void global_cleanup()
{
	if (g_kdTree_rawData)
	{
		g_kdTree_rawData->clear();
		SAFE_DELETE(g_kdTree_rawData);
	//	g_kdTree_rawData = NULL;
	}
	
	if (g_kdTree_1flat){
		g_kdTree_1flat->clear();
		SAFE_DELETE(g_kdTree_1flat);
		//g_kdTree_1flat = NULL;
	}

	if (g_kdTree_2flat){
		g_kdTree_2flat->clear();
		SAFE_DELETE(g_kdTree_2flat);
		//g_kdTree_2flat = NULL;
	}

	// Clear p-flat lists
	
	if (!g_0flat_list.empty())
	{
		for (vector<zeroFlatPt*>::iterator IT = g_0flat_list.begin(); IT != g_0flat_list.end(); ++IT)
			SAFE_DELETE(*IT);
		g_0flat_list.clear();
	}

	if (!g_1flat_list.empty())
	{
		for (vector<pFlatPt*>::iterator IT = g_1flat_list.begin(); IT != g_1flat_list.end(); ++IT)
			SAFE_DELETE(*IT);
		g_1flat_list.clear();
	}

	if (!g_2flat_list.empty())
	{
		for (vector<pFlatPt*>::iterator IT = g_2flat_list.begin(); IT != g_2flat_list.end(); ++IT)
			SAFE_DELETE(*IT);
		g_2flat_list.clear();
	}

	if (!g_1flat_listS2.empty())
	{
		for (vector<pFlatPt*>::iterator IT = g_1flat_listS2.begin(); IT != g_1flat_listS2.end(); ++IT)
			SAFE_DELETE(*IT);
		g_1flat_listS2.clear();
	}

	if (!g_2flat_listS2.empty())
	{
		for (vector<pFlatPt*>::iterator IT = g_2flat_listS2.begin(); IT != g_2flat_listS2.end(); ++IT)
			SAFE_DELETE(*IT);
		g_2flat_listS2.clear();
	}

	if (!g_0flat_query_result.empty())
	{
		// Remove all elements but no deletion since all nodes have been deleted in previous operations
		for (size_t i = 0; i < g_0flat_query_result.size(); i++)
		{
			g_0flat_query_result[i].clear();
		}
		g_0flat_query_result.clear();
	}

	if (!g_1flat_query_result.empty())
	{
		// Remove all elements but no deletion since all nodes have been deleted in previous operations
		for (size_t i = 0; i < g_1flat_query_result.size(); i++)
		{
			g_1flat_query_result[i].clear();
		}
		g_1flat_query_result.clear();
	}

	if (!g_2flat_query_result.empty())
	{
		// Remove all elements but no deletion since all nodes have been deleted in previous operations
		for (size_t i = 0; i < g_2flat_query_result.size(); i++)
		{
			g_2flat_query_result[i].clear();
		}
		g_2flat_query_result.clear();
	}

	// pcp
	if (g_pcp)
	{
		SAFE_DELETE(g_pcp);
	}
}

float RandomFloat(){
	if( !g_isGlobalInit )
	{
		global_init();
	}
	return g_randGen.rand();
}
u_int RandomUInt() {
	if( !g_isGlobalInit )
	{
		global_init();
	}
	return g_randGen.randInt();
}

float box_muller(float m, float s)	/* normal random variate generator */
{				        /* mean m, standard deviation s */
	float x1, x2, w, y1;
	static float y2;
	static int use_last = 0;

	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0f * RandomFloat() - 1.0f;
			x2 = 2.0f * RandomFloat() - 1.0f;
			w = x1 * x1 + x2 * x2;
		} while (w >= 1.0);

		w = sqrt((-2.0 * log(w)) / w);
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	return(m + y1 * s);
}

void RandomNormDistr(std::vector<float>& y, int num, float m, float s)
{
	if (!g_isGlobalInit)
	{
		global_init();
	}
	// Do Box-muller transformation
	y.resize(num);
	for (size_t i = 0; i < y.size(); i++)
	{
		y[i] = box_muller(m, s);
	}
}

