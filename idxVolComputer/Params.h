#ifndef PARAMS_H
#define PARAMS_H

#include "MyVectors.h"
#include "myProgDef.h"
//#include "opencv.hpp"
#include "VolumeData.h"
//#include "Vertex.h"

enum VOL_RENDER_SHADING_TECH
{
	VRS_DIROCC_VFS = 0, // directional occlusion shading with vertex+fragment shaders
	VRS_DIROCC_CS, // directional occlusion shading with compute shaders
	VRS_PHONG, // traditional phong shading. WARNING: the gradient is not defined for ND TF case!
	VRS_DIROCC_CS_1D // 1D TF version of dirOcc shading (For testing only)!
};


class Params
{
public:
	Params()
	{
		_pcpDispBufSize = INTVECTOR2(2048, 1024);
		_minPCInselbergPt = FLOATVECTOR2(-0.15f, -1.25f);
		_maxPCInselbergPt = FLOATVECTOR2(1.15f, 1.25f);
		_isTimeVaryingData = false;

		// Setup param for the vector field on 2D scatterplot
	/*	_scVec2DParam.num_seeds =

	;
		_scVec2DParam.max_pts_per_curve = 100;
		_scVec2DParam.rej_sampl_M = 50.0;
		_scVec2DParam.max_zero_velo_pts = 2;
		_scVec2DParam.streamline_stepsize = 1.0f;*/
		//_samplePortion = 0.002f;//0.005f;
		//_sampleNeighborSize = 5;
		//_numNN_valDomain = 1000;// 1000;
		//// KD trees
		//_kdTree_rawData = NULL;
		//_kdTree_1flat = NULL;
		//_kdTree_2flat = NULL;
		//// pcp 
		//_pcp = NULL;
		//// size of plots
		////_pcpDispBufWidth = 4096;
		////_pcpDispBufHeight = 4096;

		//_params.setPCPDispBufSize(4096, 4096);

		//_splomDispBufWidth = 2048;
		//_splomDispBufHeight = 2048;
		// rendering colors: use colorbrewer 3-class qualitative scheme
		Pcp_basic_layer_color(INTVECTOR4(0x8d, 0xa0, 0xcb, 60));
		Pcp_1flat_layer_color(INTVECTOR4(0x1f, 0x78, 0xb4, 60));
		_pcp_2flat_layer_color = INTVECTOR4(0x66, 0xc2, 0xa5, 60);
		Sc_basic_layer_color(INTVECTOR4(0, 0, 255, 100));
		Sc_1flat_layer_color(INTVECTOR4(0xe4, 0x1a, 0xdf, 60));
		// width of the scatterplot drawing pen
		int num_pts_short_edge = 500;
		//_sc_pen_width = int(1.0 / double(num_pts_short_edge) * double(MIN(_splomDispBufWidth, _splomDispBufHeight)));
		//_pcp_pen_width = /*2 * */int(1.0 / double(num_pts_short_edge) * double(_params.getPCPDispBufSize().minVal())); // Should config point size in the config menu.
		// alpha of pcp pen
		Pcp_pen_alpha(1.0f);
		_sc_pen_alpha = 1.0f;
		// pcp render mapping param
		Pcp_render_scale(1.0);
		Splom_render_scale(1.0);
		// Default to none for subspace selection
		FirstAxis_subspace(-1);
		// Default to not selecting any layer!
		Pcp_selectedLayer(0);
		Pcp_layer_weights(FLOATVECTOR3(1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 3.0f));
		// Max log val 
		Max_log_val(1.0);
		// number of samples to process per update
		Num_sample_per_update(10000);
		// Controls draw one trend for how many samples
		TrendData_stride(20);
		// Set pcp repeat mode
		Use_repeat_pcp(false);
		// Set use full space
		ProcessFullSpace(true);

		// Threshold to determine at least how many records we have in data so that we need to sample only part of the data
		Min_rec_num_thres_for_sample(40000);

		// HACK: set pre-defined colormap for brushes
		//_brushColormap.push_back(FLOATVECTOR3(0, 0, 1));
		//_brushColormap.push_back(FLOATVECTOR3(1, 0, 0));
		//_brushColormap.push_back(FLOATVECTOR3(0, 1, 0));
		//_brushColormap.push_back(FLOATVECTOR3(0, 1, 1));
		//_brushColormap.push_back(FLOATVECTOR3(1, 0, 1));
		//_brushColormap.push_back(FLOATVECTOR3(1, 1, 0));

		// set threshold to 0 (show all p-flats)
		PFlatIdPt_thres(0.0);

		// Set display area for Inselberg's coordinates
		//_minPCInselbergPt = FLOATVECTOR2(-0.1f, -0.5f);
		//_maxPCInselbergPt = FLOATVECTOR2(1.1f, 1.5f);
		// Use the config below for 2 attributes
		//_minPCInselbergPt = FLOATVECTOR2(-0.15f, -1.25f);
		//_maxPCInselbergPt = FLOATVECTOR2(1.15f, 1.25f);

		// No vector data by default
		IsVecData(false);
		// Not trajectory data by default
		IsTrajData(false);
		// use transform for idx points and pcp polylines
		//_xformMethod = XF_NONE;//XF_IDX_ONLY;
		XformMethod(XF_IDX_LINE);
		// use scaling in transformation? ALWAYS true unless you want to see the unscaled effects!
		IsScaleXform(true);
		// Show full line for each data?
		IsFullPolyline(false);
		//_xformMethod = XF_IDX_ONLY;
		// set tehe nearest neighbor query methos
		NnQueryMethod(NNQ_KNN);
		//_nnQueryMethod = NNQ_RBALL;
		// How many neighborhood configurations needed for the uncertian point data? 
		Uncertain_samples_per_distr(20);
		// Draw anti-aliasing lines? By default: no
		AntialiasLines(false);
		// To generate figures, use supersampling
		SuperSampleLines(false); // 
		SuperSampleNum(9); // number of samples per sampling point

		// Set PC render mode
		//_pcp_RenderMode = PR_LINE_IDX_PT;
		Pcp_RenderMode(PR_LINE_IDX_PT);
		// set compare mode
		Pcp_compare_mode(0);

		// number of samples to skip
		Pcp_num_samples_skip_DS_LINE_MODE(30);
		Pcp_num_samples_skip_LINE_IDXPT_MODE(1);
		//_pcp_num_samples_skip_LINE_IDXPT_MODE = _pcp_num_samples_skip_DS_LINE_MODE;// 10;

		Pcp_downsample_max_log_val(1.0);
		// Set size of density plots
		Pcp_densPlotBufSize(INTVECTOR2(1024, 2048));
		IsHDR(false);
		CmColormapLen((IsHDR()) ? 1024 : 256);
		// Set color scheme of PCP visualization
		VisColorScheme(CS_LIGHT);
		// number of segments of a curve within half the distance of 2 neighboring axes
		Num_segCurve_halfdist(40);
		// scaling factors at bounds
		BoundScale(1.0);

		DensMapCornerThres(0.01f);
		DensMapGamma(2.2f);

		// use downsampling
		UseMonteCarlo(false);
		DownSampleDim(INTVECTOR3(50, 50, 10));
		DownSamplePortion(FLOATVECTOR3(0.1f, 0.1f, 0.1f));

		// run timing
		Timing(false);
		// normalize data
		IsNormalizeData(true);

		// clustering
		doClustering(false);
		numClusters(5);

		// Continuous indexed points
		IsComputeContIdxPts(false);
		IsIdxPtDensityMap(true);
		density1flatImgsPath("."); // using the current path as density images path
		Is2flatIdxPtDensityMap(false);
		density2flatImgsPath("."); // using the current path as density images path
		SpatialNeighborNum(100);
		gridSampleRate(0.2f);

		// The margin of the PC canvas
		setPcCanvasMarginMinPortion(0.12f);

		// Set colors for indexed points difference mode
		colorBothIdxPts(FLOATVECTOR3(5, 113, 176));
		colorVdIdxPts(FLOATVECTOR3(202, 0, 32));

		// volume rendering parameters
		volRenSampleRate(1.0f);
		volRenShadeTech(VRS_DIROCC_CS);
		voxelRescaleFactors(FLOATVECTOR3(1.0f, 1.0f, 1.0f));
		volRenIsUseExtOptimization(false);
		// Query controls
		isQueryPC(true);
		isQuerySPLOM(true);

		// Coordinates drawing controls
		isDrawPCAxisHorizontal(true);

		// Is compare value domain?
		IsCompareValueDomain(true);
	}

	INTVECTOR2 getPCPDispBufSize() const
	{
		return _pcpDispBufSize;
	}

	void setPCPDispBufSize(int width, int height)
	{
		_pcpDispBufSize = INTVECTOR2(width, height);
	}


	float pcCanvasMarginMinPortion() const { return _minPCcanvasPortion; }
	void setPcCanvasMarginMinPortion(float minPortion) { _minPCcanvasPortion = minPortion; }

	FLOATVECTOR2 minPCInselbergPt() const { return _minPCInselbergPt; }
	void setMinPCInselbergPt(const FLOATVECTOR2& minPt) { _minPCInselbergPt = minPt; }

	FLOATVECTOR2 maxPCInselbergPt() const { return _maxPCInselbergPt; }
	void setMaxPCInselbergPt(const FLOATVECTOR2& maxPt) { _maxPCInselbergPt = maxPt; }

	bool isTimeVarying() const { return _isTimeVaryingData; }
	void setIsTimeVarying(bool isTimeVarying) { _isTimeVaryingData = isTimeVarying; }

	std::vector<zeroFlatPt*> zeroflat_list() const { return _0flat_list; }
	void zeroflat_list(std::vector<zeroFlatPt*> val) { _0flat_list = val; }


	KD<zeroFlatPt*>* KdTree_rawData() const { return _kdTree_rawData; }
	void KdTree_rawData(KD<zeroFlatPt*>* val) { _kdTree_rawData = val; }

	std::vector<pFlatPt*> oneflat_list() const { return _1flat_list; }
	void oneflat_list(std::vector<pFlatPt*> val) { _1flat_list = val; }


	std::vector<pFlatPt*> twoflat_list() const { return _2flat_list; }
	void twoflat_list(std::vector<pFlatPt*> val) { _2flat_list = val; }

	KD<pFlatPt*>* KdTree_1flat() const { return _kdTree_1flat; }
	void KdTree_1flat(KD<pFlatPt*>* val) { _kdTree_1flat = val; }

	KD<pFlatPt*>* KdTree_2flat() const { return _kdTree_2flat; }
	void KdTree_2flat(KD<pFlatPt*>* val) { _kdTree_2flat = val; }


	std::vector<std::vector<UINT64>> Xpcp_query_result() const { return _xpcp_query_result; }
	void Xpcp_query_result(std::vector<std::vector<UINT64>> val) { _xpcp_query_result = val; }

	std::vector<std::vector<zeroFlatPt*>> zeroflat_query_result() const { return _0flat_query_result; }
	void zeroflat_query_result(std::vector<std::vector<zeroFlatPt*>> val) { _0flat_query_result = val; }

	std::vector<std::vector<pFlatPt*>> oneflat_query_result() const { return _1flat_query_result; }
	void oneflat_query_result(std::vector<std::vector<pFlatPt*>> val) { _1flat_query_result = val; }


	std::vector<FLOATVECTOR3> twoflat_brushColors() const { return _2flat_brushColors; }
	void twoflat_brushColors(std::vector<FLOATVECTOR3> val) { _2flat_brushColors = val; }


	INTVECTOR4 Pcp_basic_layer_color() const { return _pcp_basic_layer_color; }
	void Pcp_basic_layer_color(INTVECTOR4 val) { _pcp_basic_layer_color = val; }

	INTVECTOR4 Pcp_1flat_layer_color() const { return _pcp_1flat_layer_color; }
	void Pcp_1flat_layer_color(INTVECTOR4 val) { _pcp_1flat_layer_color = val; }

	INTVECTOR4 Sc_basic_layer_color() const { return _sc_basic_layer_color; }
	void Sc_basic_layer_color(INTVECTOR4 val) { _sc_basic_layer_color = val; }

	INTVECTOR4 Sc_1flat_layer_color() const { return _sc_1flat_layer_color; }
	void Sc_1flat_layer_color(INTVECTOR4 val) { _sc_1flat_layer_color = val; }


	std::vector<INTVECTOR4> Pcp_1flat_colormap() const { return _pcp_1flat_colormap; }
	void Pcp_1flat_colormap(std::vector<INTVECTOR4> val) { _pcp_1flat_colormap = val; }


	std::vector<INTVECTOR4> Pcp_2flat_colormap() const { return _pcp_2flat_colormap; }
	void Pcp_2flat_colormap(std::vector<INTVECTOR4> val) { _pcp_2flat_colormap = val; }

	UINT64 Min_rec_num_thres_for_sample() const { return _min_rec_num_thres_for_sample; }
	void Min_rec_num_thres_for_sample(UINT64 val) { _min_rec_num_thres_for_sample = val; }


	int FirstAxis_subspace() const { return _firstAxis_subspace; }
	void FirstAxis_subspace(int val) { _firstAxis_subspace = val; }

	bool isDrawPCAxisHorizontal() const { return _drawCoordHor; }
	void isDrawPCAxisHorizontal(bool val) { _drawCoordHor = val; }

	int Pcp_selectedLayer() const { return _pcp_selectedLayer; }
	void Pcp_selectedLayer(int val) { _pcp_selectedLayer = val; }

	int Pcp_compare_mode() const { return _pcp_compareMode; }
	void Pcp_compare_mode(int mode) { _pcp_compareMode = mode; }
	void Pcp_compare_mode_cycle() {
		if (_pcp_compareMode == 2)
			_pcp_compareMode = 0;
		else
			_pcp_compareMode++;
	}

	int Sc_pen_width() const { return _sc_pen_width; }
	void Sc_pen_width(int val) { _sc_pen_width = val; }


	int Pcp_pen_width() const { return _pcp_pen_width; }
	void Pcp_pen_width(int val) { _pcp_pen_width = val; }
	float Pcp_pen_alpha() const { return _pcp_pen_alpha; }
	void Pcp_pen_alpha(float val) { _pcp_pen_alpha = val; }
	double Pcp_render_scale() const { return _pcp_render_scale; }
	void Pcp_render_scale(double val) { _pcp_render_scale = val; }
	int Num_sample_per_update() const { return _num_sample_per_update; }
	void Num_sample_per_update(int val) { _num_sample_per_update = val; }
	double Splom_render_scale() const { return _splom_render_scale; }
	void Splom_render_scale(double val) { _splom_render_scale = val; }
	double PFlatIdPt_thres() const { return _pFlatIdPt_thres; }
	void PFlatIdPt_thres(double val) { _pFlatIdPt_thres = val; }

	int TrendData_stride() const { return _trendData_stride; }
	void TrendData_stride(int val) { _trendData_stride = val; }

	FLOATVECTOR3 Pcp_layer_weights() const { return _pcp_layer_weights; }
	void Pcp_layer_weights(FLOATVECTOR3 val) { _pcp_layer_weights = val; }

	double Max_log_val() const { return _max_log_val; }
	void Max_log_val(double val) { _max_log_val = val; }

	bool Use_repeat_pcp() const { return _use_repeat_pcp; }
	void Use_repeat_pcp(bool val) { _use_repeat_pcp = val; }

	std::vector<std::string> Attrib_names() const { return _attrib_names; }
	void Attrib_names(std::vector<std::string> val) { _attrib_names = val; }

	std::vector<FLOATVECTOR3> BrushColormap() const { return _brushColormap; }
	void BrushColormap(std::vector<FLOATVECTOR3> val) { _brushColormap = val; }

	bool AntialiasLines() const { return _antialiasLines; }
	void AntialiasLines(bool val) { _antialiasLines = val; }

	bool SuperSampleLines() const { return _superSampleLines; }
	void SuperSampleLines(bool val) { _superSampleLines = val; }

	int SuperSampleNum() const { return _superSampleNum; }
	void SuperSampleNum(int val) { _superSampleNum = val; }

	bool ProcessFullSpace() const { return _processFullSpace; }
	void ProcessFullSpace(bool val) { _processFullSpace = val; }

	int Uncertain_samples_per_distr() const { return _uncertain_samples_per_distr; }
	void Uncertain_samples_per_distr(int val) { _uncertain_samples_per_distr = val; }

	DAT_VOL_INFO DatVolInfo() const { return _datVolInfo; }
	void DatVolInfo(DAT_VOL_INFO val) { _datVolInfo = val; }

	bool IsVecData() const { return _isVecData; }
	void IsVecData(bool val) { _isVecData = val; }

	bool IsTrajData() const { return _isTrajData; }
	void IsTrajData(bool val) { _isTrajData = val; }

	XFORM_OPTION XformMethod() const { return _xformMethod; }
	void XformMethod(XFORM_OPTION val) { _xformMethod = val; }

	bool IsFullPolyline() const { return _isFullPolyline; }
	void IsFullPolyline(bool val) { _isFullPolyline = val; }

	bool IsScaleXform() const { return _isScaleXform; }
	void IsScaleXform(bool val) { _isScaleXform = val; }

	NN_QUERY_METHOD NnQueryMethod() const { return _nnQueryMethod; }
	void NnQueryMethod(NN_QUERY_METHOD val) { _nnQueryMethod = val; }

	PC_RENDER_MODE Pcp_RenderMode() const { return _pcp_RenderMode; }
	void Pcp_RenderMode(PC_RENDER_MODE val) { _pcp_RenderMode = val; }

	INTVECTOR2 Pcp_densPlotBufSize() const { return _pcp_densPlotBufSize; }
	void Pcp_densPlotBufSize(INTVECTOR2 val) { _pcp_densPlotBufSize = val; }

	int Pcp_num_samples_skip_DS_LINE_MODE() const { return _pcp_num_samples_skip_DS_LINE_MODE; }
	void Pcp_num_samples_skip_DS_LINE_MODE(int val) { _pcp_num_samples_skip_DS_LINE_MODE = val; }

	int Pcp_num_samples_skip_LINE_IDXPT_MODE() const { return _pcp_num_samples_skip_LINE_IDXPT_MODE; }
	void Pcp_num_samples_skip_LINE_IDXPT_MODE(int val) { _pcp_num_samples_skip_LINE_IDXPT_MODE = val; }

	double Pcp_downsample_max_log_val() const { return _pcp_downsample_max_log_val; }
	void Pcp_downsample_max_log_val(double val) { _pcp_downsample_max_log_val = val; }

	//std::vector<MyAdvColor> CmColormap() const { return _cmColormap; }
	//void CmColormap(std::vector<MyAdvColor> val) { _cmColormap = val; }

	//std::vector<MyAdvColorNode> CmCtrlPoints() const { return _cmCtrlPoints; }
	//void CmCtrlPoints(std::vector<MyAdvColorNode> val) { _cmCtrlPoints = val; }

	//std::vector<QGradientStops> Cm_pcp_tfs() const { return _cm_pcp_tfs; }
	//void Cm_pcp_tfs(std::vector<QGradientStops> val) { _cm_pcp_tfs = val; }

	int CmColormapLen() const { return _cmColormapLen; }
	void CmColormapLen(int val) { _cmColormapLen = val; }

	COLOR_SCHEME VisColorScheme() const { return _visColorScheme; }
	void VisColorScheme(COLOR_SCHEME val) { _visColorScheme = val; }

	int Num_segCurve_halfdist() const { return _num_segCurve_halfdist; }
	void Num_segCurve_halfdist(int val) { _num_segCurve_halfdist = val; }

	bool IsHDR() const { return _isHDR; }
	void IsHDR(bool val) { _isHDR = val; }

	double BoundScale() const { return _boundScale; }
	void BoundScale(double val) { _boundScale = val; }

	bool IsDrawMaxLocs() const { return _isDrawMaxLocs; }
	void IsDrawMaxLocs(bool val) { _isDrawMaxLocs = val; }

	float DensMapGamma() const { return _densMapGamma; }
	void DensMapGamma(float val) { _densMapGamma = val; }

	float DensMapCornerThres() const { return _densMapCornerThres; }
	void DensMapCornerThres(float val) { _densMapCornerThres = val; }

	std::vector<INTVECTOR4> Density_map_colormap() const { return _density_map_colormap; }
	void Density_map_colormap(std::vector<INTVECTOR4> val) { _density_map_colormap = val; }


	//const std::vector<QImage>& DensMapRenderMask() const { return _densMapRenderMask; }
	//void DensMapRenderMask(const std::vector<QImage>& val) { _densMapRenderMask = val; }

	INTVECTOR2 Pcp_densPlotAccumBufSize() const { return _pcp_densPlotAccumBufSize; }
	void Pcp_densPlotAccumBufSize(INTVECTOR2 val) { _pcp_densPlotAccumBufSize = val; }

	INTVECTOR2 Pcp_densPlotDispBufSize() const { return _pcp_densPlotDispBufSize; }
	void Pcp_densPlotDispBufSize(INTVECTOR2 val) { _pcp_densPlotDispBufSize = val; }

	//const std::vector<cv::Mat>& DensMapROILocs() const { return _densMapROILocs; }
	//void DensMapROILocs(const std::vector<cv::Mat>& val) { _densMapROILocs = val; }

	bool UseMonteCarlo() const { return _useMonteCarlo; }
	void UseMonteCarlo(bool useMonteCarlo) { _useMonteCarlo = useMonteCarlo; }

	// Parameters of the density map of continuous indexed points
	bool IsIdxPtDensityMap() const { return _isIdxPtDensityMap; }
	void IsIdxPtDensityMap(bool val) { _isIdxPtDensityMap = val; }

	bool Is2flatIdxPtDensityMap() const { return _is2flatIdxPtDensityMap; }
	void Is2flatIdxPtDensityMap(bool val) { _is2flatIdxPtDensityMap = val; }

	bool IsComputeContIdxPts() const { return _isCompContIdxPts; }
	void IsComputeContIdxPts(bool val) { _isCompContIdxPts = val; }

	bool IsCompareValueDomain() const { return _isCompareValDomain; }
	void IsCompareValueDomain(bool val) { _isCompareValDomain = val; }

	FLOATVECTOR3 SpatialNeighborSize() const { return _spatialNeighborSize; }
	void SpatialNeighborSize(const FLOATVECTOR3& snSize) { _spatialNeighborSize = snSize; }
	void SpatialNeighborSizeX(float xs) { _spatialNeighborSize.x = xs; }
	void SpatialNeighborSizeY(float ys) { _spatialNeighborSize.y = ys; }
	void SpatialNeighborSizeZ(float zs) { _spatialNeighborSize.z = zs; }

	int SpatialNeighborNum() const { return _spatialNeighborNum; }
	void SpatialNeighborNum(int num) { _spatialNeighborNum = num; }

	const INTVECTOR3& DownSampleDim() const { return _downSampleDim; }
	void DownSampleDim(const INTVECTOR3& dim) { _downSampleDim = dim; }

	const FLOATVECTOR3& DownSamplePortion() const { return _downSamplePortion; }
	void DownSamplePortion(const FLOATVECTOR3& portion) { _downSamplePortion = portion; }

	const std::vector<UINT64>& OutlierIds() const { return _outlierIds; }
	void OutlierIds(const std::vector<UINT64>& val) { _outlierIds = val; }

	const bool Timing() const { return _runTiming; }
	void Timing(const bool& val) { _runTiming = val; }

	const bool IsNormalizeData() const { return _isNormalizeData; }
	void IsNormalizeData(const bool& val) { _isNormalizeData = val; }

	const bool doClustering() const { return _doClustering; }
	void doClustering(const bool& val) { _doClustering = val; }

	const int numClusters() const { return _numClusters; }
	void numClusters(const int& numClusters) { _numClusters = numClusters; }

	const std::vector<int>& clusterLabelData() const { return _labelData; }
	void clusterLabelData(const std::vector<int>& data) { _labelData = data; }

	const std::string& fileDir() const { return _fileDir; }
	void fileDir(const std::string& dir) { _fileDir = dir; }

	const std::string& fileName() const { return _fileName; }
	void fileName(const std::string& fn) { _fileName = fn; }

	// SFC brushes
	const std::vector<FLOATVECTOR2>& brushedSFCranges() const { return  _brushedSFCranges; }// rectangular range selected in the SFC
	void setBrushedSFCranges(std::vector<FLOATVECTOR2>& range) { _brushedSFCranges = range; }

	const std::vector<FLOATVECTOR3>& brushColors() const { return _brushCols; }
	void setBrushColors(std::vector<FLOATVECTOR3>& brushCols) { _brushCols = brushCols; }

	// Brushed volumes
	std::vector<VolumeData*>& brushVols_unsafe() { return _brushVols; }
	const std::vector<VolumeData*>& brushVols() const { return _brushVols; }
	void setBrushVols(const std::vector<VolumeData*>& val) { _brushVols = val; }

	// Indexed points volumes: 1-flats
	std::vector<VolumeData*>& idxPtVols_unsafe() { return _idxPtVols; }
	const std::vector<VolumeData*>& idxPtVols() const { return _idxPtVols; }
	void setIdxPtVols(const std::vector<VolumeData*>& val) { _idxPtVols = val; }
	void setIdxPtVol(VolumeData* val, int i) {
		if (i >= 0 && i < _idxPtVols.size())
			_idxPtVols[i] = val;
	}

	// Indexed points volumes Set 2: 1-flats
	std::vector<VolumeData*>& idxPtVolsS2_unsafe() { return _idxPtVolsS2; }
	const std::vector<VolumeData*>& idxPtVolsS2() const { return _idxPtVolsS2; }
	void setIdxPtVolsS2(const std::vector<VolumeData*>& val) { _idxPtVolsS2 = val; }
	void setIdxPtVolS2(VolumeData* val, int i) {
		if (i >= 0 && i < _idxPtVolsS2.size())
			_idxPtVolsS2[i] = val;
	}

	// 2-flats indexed points volumes
	std::vector<VolumeData*>& idxPt2flatsVols_unsafe() { return _idxPtVols2Flats; }
	const std::vector<VolumeData*>& idxPt2flatsVols() const { return _idxPtVols2Flats; }
	void setIdxPt2flatsVols(const std::vector<VolumeData*>& val) { _idxPtVols2Flats = val; }
	void setIdxPt2flatVol(VolumeData* val, int i) {
		if (i >= 0 && i < _idxPtVols2Flats.size())
			_idxPtVols2Flats[i] = val;
	}

	// Multi-attribute volumes
	std::vector<VolumeData*>& dataVols_unsafe() { return _dataVols; }
	const std::vector<VolumeData*>& dataVols() const { return _dataVols; }
	void setDataVols(const std::vector<VolumeData*>& val) { _dataVols = val; }
	void setDataVol(VolumeData* val, int i) {
		if (i >= 0 && i < _dataVols.size())
			_dataVols[i] = val;
	}

	const float gridSampleRate() const { return _gridSampleRate; }
	void gridSampleRate(const float& sampleRate) { _gridSampleRate = sampleRate; }

	const std::vector<FLOATVECTOR3>& getBrushColormap() const { return _brushColormap; }
	void setBrushColormap(const std::vector<FLOATVECTOR3>& colormap) { _brushColormap = colormap; }

	const int VolBlockSize() const { return _volBlockSize; }
	void VolBlockSize(int volBkSize) { _volBlockSize = MAX(volBkSize, 2); } // min block size is 2*2*2!

	void setOctreeLevel(int level) { _currOctreeLevel = level; }
	int octreeLevel() const { return _currOctreeLevel; }

	//const std::vector<Vertex>& vertices() const {
	//	return _vertices;
	//}
	//void setVertices(const std::vector<Vertex>& vertices) {
	//	_vertices = vertices;
	//}

	const std::vector<std::vector<float> >& pointData() const { return _pointData; }
	void setPointData(const std::vector<std::vector<float> >& pointData) { _pointData = pointData; }

	const std::vector<std::vector<float> >& ensembleData() const { return _ensembleData; }
	void setEnsembleData(const std::vector<std::vector<float> >& val) { _ensembleData = val; }

	void setNumRuns(int numRuns) { _numRuns = numRuns; }
	int numRuns() const { return _numRuns; }

	const FLOATVECTOR3 colorVdIdxPts() const { return _colorVdIdxPtOnly; }
	void colorVdIdxPts(const FLOATVECTOR3& val) { _colorVdIdxPtOnly = val; }

	const FLOATVECTOR3 colorBothIdxPts() const { return _colorBothIdxPts; }
	void colorBothIdxPts(const FLOATVECTOR3& val) { _colorBothIdxPts = val; }

	const float volRenSampleRate() const { return _volRenSampleRate; }
	void volRenSampleRate(const float& val) { _volRenSampleRate = val; }

	const bool volRenIsUseExtOptimization() const { return _volRenUseExtOptimization; }
	void volRenIsUseExtOptimization(const bool& val) { _volRenUseExtOptimization = val; }

	const FLOATVECTOR3 voxelRescaleFactors() const { return _voxelRescaleFactors; }
	void voxelRescaleFactors(const FLOATVECTOR3& val) { _voxelRescaleFactors = val; }

	const VOL_RENDER_SHADING_TECH volRenShadeTech() const { return _vrShadeTech; }
	void volRenShadeTech(const VOL_RENDER_SHADING_TECH& val) { _vrShadeTech = val; }

	const bool isQueryPC() const { return _isQueryPC; }
	void isQueryPC(const bool& val) { _isQueryPC = val; }

	const bool isQuerySPLOM() const { return _isQuerySPLOM; }
	void isQuerySPLOM(const bool& val) { _isQuerySPLOM = val; }

	// Accessors of eigenvolumes
	void setMajEigVols(const std::vector<VolumeData*>& vols) { _majEigVols = vols; }
	const std::vector<VolumeData*>& majEigVols() const { return _majEigVols; }

	void setSecEigVols(const std::vector<VolumeData*>& vols) { _secEigVols = vols; }
	const std::vector<VolumeData*>& secEigVols() const { return _secEigVols; }
	// density maps of indexed points
	void density1flatImgsPath(const std::string& val) { _density1flatImgsPath = val; }
	const std::string& density1flatImgsPath() const { return _density1flatImgsPath; }

	void density2flatImgsPath(const std::string& val) { _density2flatImgsPath = val; }
	const std::string& density2flatImgsPath() const { return _density2flatImgsPath; }
	//std::vector<std::vector<float>>& global_sampled_raw_unsafe() { return g_sampled_raw; }
	//const std::vector<std::vector<float>>& global_sampled_raw() const { return g_sampled_raw; }

	const int Ndim_orgData() const { return _Ndim; }
	void Ndim_orgData(const int& Ndim) { _Ndim = Ndim; }
private:

	INTVECTOR2          _pcpDispBufSize;
	// Display area in Inselberg's coordinates for parallel coordinates
	FLOATVECTOR2 _minPCInselbergPt;
	FLOATVECTOR2 _maxPCInselbergPt;
	// Margin of the PC image plane in terms of portion
	float        _minPCcanvasPortion;
	bool                _isTimeVaryingData;

	// It would be easier to have the KDtrees as global variables
	std::vector<zeroFlatPt*> _0flat_list;
	KD<zeroFlatPt*>*       _kdTree_rawData;
	// pFlatPt representation for p-flats
	std::vector<pFlatPt*>  _1flat_list;
	std::vector<pFlatPt*>  _2flat_list;
	KD<pFlatPt*>*          _kdTree_1flat;
	KD<pFlatPt*>*          _kdTree_2flat;
	// Query result
	std::vector< std::vector<UINT64>>    _xpcp_query_result;	// the polyline layer query result if we don't have KD-tree
	std::vector< std::vector<zeroFlatPt*> > _0flat_query_result; // an alternative of polyline layer if we have KD-tree
	std::vector< std::vector<pFlatPt*> > _1flat_query_result;
	std::vector< std::vector<pFlatPt*> > _2flat_query_result;
	// Color associated with query brushes
	std::vector< FLOATVECTOR3 > _xpcp_brushColors;
	std::vector< FLOATVECTOR3 > _1flat_brushColors;
	std::vector< FLOATVECTOR3 > _2flat_brushColors;

	// colors for rendering
	INTVECTOR4 _pcp_basic_layer_color;
	INTVECTOR4 _pcp_1flat_layer_color;
	INTVECTOR4 _pcp_2flat_layer_color;
	INTVECTOR4 _sc_basic_layer_color;
	INTVECTOR4 _sc_1flat_layer_color;
	std::vector<INTVECTOR4> _pcp_1flat_colormap;
	std::vector<INTVECTOR4> _pcp_2flat_colormap;

	std::vector<INTVECTOR4> _density_map_colormap;

	UINT64     _min_rec_num_thres_for_sample;
	int        _firstAxis_subspace;
	int        _pcp_selectedLayer;
	int        _pcp_compareMode; // 0 - normal mode (show current idx pt config); 1 - overlay diff mode (red- only 
	int        _sc_pen_width;
	int        _pcp_pen_width;
	float      _pcp_pen_alpha;
	float      _sc_pen_alpha;
	double     _pcp_render_scale; // the parameter controlling accumulation to color mapping for pcp rendering 
	int        _num_sample_per_update;
	double     _splom_render_scale;
	double     _pFlatIdPt_thres; // Threshold for the p-flat indexed point strength
	int        _trendData_stride;
	FLOATVECTOR3 _pcp_layer_weights; // weights for blending layers of the PCP
	double     _max_log_val; // the maximum possible log value (For now, use log(num_records))
	bool       _use_repeat_pcp; // Are we using repeating pcp (by default, no)

	// Volume rendering parameters
	float _volRenSampleRate; // sampling rate for volume rendering
	FLOATVECTOR3 _voxelRescaleFactors; // rescale factors of the volume
	VOL_RENDER_SHADING_TECH _vrShadeTech; // shading technique
	bool _volRenUseExtOptimization; // is use extinction optimization
	// Attribute names
	std::vector<std::string> _attrib_names; // names of all attributes
	// 
	std::vector<FLOATVECTOR3> _brushColormap; // a predefined colormap for brushes

	bool       _antialiasLines;  // draw anti-aliased lines???
	bool       _superSampleLines; // draw multiple samples for a point
	int        _superSampleNum; // sample number for supersampling

	bool       _processFullSpace; // Are we using the full space to find NN & eigenvectors ???
	int        _uncertain_samples_per_distr; // number of samplings to do to get different 
	// data information
	DAT_VOL_INFO  _datVolInfo; // the variable that records the dat vector file
	bool          _isVecData; // if it's a vector data
	bool          _isTrajData; // if it's a trajectory data
	// Tranformation related
	XFORM_OPTION  _xformMethod; // The transformation method for indexed points and/or polylines
	bool          _isFullPolyline; // Draw full line or line segment? Line segment by default.
	bool          _isScaleXform; // use scaling function in transformation? !!! DO NOT turn off other than to show the unscaled effects.
	NN_QUERY_METHOD _nnQueryMethod; // the method for neighborhood query

	// PCP render modes
	PC_RENDER_MODE _pcp_RenderMode;
	INTVECTOR2     _pcp_densPlotBufSize;
	int            _pcp_num_samples_skip_DS_LINE_MODE; // how many samples we skip when drawing the simplified pcp
	int            _pcp_num_samples_skip_LINE_IDXPT_MODE; // how many samples to skip when drawing line+idx pt mode?
	double         _pcp_downsample_max_log_val;
	//// colormap variables
	//std::vector<MyAdvColor>     _cmColormap; // the current colormap as an global parameter
	//std::vector<MyAdvColorNode> _cmCtrlPoints; // control points for the colormap

	//std::vector<QGradientStops>  _cm_pcp_tfs; // base layer colormaps 
	int                         _cmColormapLen;
	COLOR_SCHEME _visColorScheme;
	int          _num_segCurve_halfdist;
	// hdr variable
	bool          _isHDR;
	double        _boundScale; // the scaling factor at bounds
	// region of interest buffer for density map
	//std::vector<cv::Mat>       _densMapROILocs;
	bool          _isDrawMaxLocs;
	// density map params
	float         _densMapGamma;	   // gamma for density map
	float         _densMapCornerThres; // threshold for corner detection in density map
	//std::vector<QImage> _densMapRenderMask; // the render mask of density map

	INTVECTOR2    _pcp_densPlotAccumBufSize;
	INTVECTOR2    _pcp_densPlotDispBufSize;

	bool          _useMonteCarlo;
	// Downsampling configuration for Monte carlo sampling 
	INTVECTOR3    _downSampleDim;
	FLOATVECTOR3  _downSamplePortion;
	// Downsampling configuation for sampling on the grid
	float         _gridSampleRate;

	// array of outlier indices
	std::vector<UINT64> _outlierIds;

	bool          _runTiming;
	// normalization when loading data
	bool          _isNormalizeData;

	// clustering the raw data?
	bool          _doClustering;
	int           _numClusters;
	// label data
	std::vector<int> _labelData;
	// the directory of files
	std::string      _fileDir;
	// file name if only one file is loaded
	std::string      _fileName;

	// continuous indexed points?
	bool             _isIdxPtDensityMap;
	bool             _is2flatIdxPtDensityMap;
	std::string      _density1flatImgsPath; // path to the density images of 1-flat indexed points
	std::string      _density2flatImgsPath; // path to the density images of 2-flat indexed points
	bool             _isCompContIdxPts;
	FLOATVECTOR3     _spatialNeighborSize;
	int              _spatialNeighborNum;

	// compare with value-domain indexed points?
	bool             _isCompareValDomain;

	// Data in global domain
	std::vector<std::vector<float>> _pointData;
	std::vector<std::vector<float>> _ensembleData;

	// Line plot/ SFC linked views
	std::vector<FLOATVECTOR2>   _brushedSFCranges; // rectangular range selected in the SFC
	std::vector<FLOATVECTOR3>   _brushCols; // colors
			// for brushes linking the data and spatial domains
	int                             _volBlockSize; // the size of a basic volume block for analysis
	int                             _currOctreeLevel; // current octree level
	// ensemble data
	int                              _numRuns; // number of runs in the folder
	//std::vector<Vertex>             _vertices;// Particle data
	std::vector<VolumeData*>        _ensembleVols; // converted from the ensemble point data
		// for brushes
	std::vector<VolumeData*>        _brushVols; // one volume for a level
	// For the indexed points volumes
	std::vector<VolumeData*>        _idxPtVols; // one volume for an attribute, each records the x,y coordinates of the idx pt in the image plane
	// Set 2 of indexed points volumes
	std::vector<VolumeData*>        _idxPtVolsS2; // one volume for an attribute, each records the x,y coordinates of the idx pt in the image plane
	// 2-flats indexed points volumes
	std::vector<VolumeData*>        _idxPtVols2Flats; // one volume for an attribute, each records the x,y coordinates of the idx pt in the image plane

	
// data volume (multi-attributes)
	std::vector<VolumeData*>        _dataVols; // one volume for an attribute
	// Colors for indexed points differences
	FLOATVECTOR3                    _colorVdIdxPtOnly; // value domain computed indexed points only
	FLOATVECTOR3                    _colorBothIdxPts;  // overlapping indexed points for value and data domain computations

	// Query controls
	bool                            _isQueryPC;
	bool                            _isQuerySPLOM;
	// coordinate drawing control
	bool                            _drawCoordHor;

	// Eigen vector volumes
	std::vector<VolumeData*>        _majEigVols;
	std::vector<VolumeData*>        _secEigVols;

	// Explicit record of the number of dimensions of the original volume data
	int                             _Ndim;
};

#endif