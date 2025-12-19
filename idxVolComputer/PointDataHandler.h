#pragma once
#include "global.h"
#include <QImage>
#include "XPCPSample.h"
#include "opencv.hpp"
// Basically, it does what VolumeSamplers do without getting involved with volumes!
enum ROI_METHOD
{
	ROI_CORNER = 0,
	ROI_LOCAL_PERCENTILE_CORNER,
	ROI_LOCAL_MAX_CORNER
};
class PlotBuilder;
class PointDataHandler
{
public:
	PointDataHandler();
	PointDataHandler(UINT64 numSamples, int dataDim, int plotDispBufW, int plotDispBufH, const std::vector<XPCPSample>& pointData);
	virtual ~PointDataHandler();

	virtual void resetPlot();
	virtual bool updatePlot();
	virtual bool updatePlotFull();
	void updatePlotLayer(int layer); // Update plot on the selected layer
	void finishPlot();
	virtual bool updatePCPlot() { return true;  }
	virtual bool updateSCnPCPlots_from_streamline() { return true; }
	// scatterplot 3D
	bool updateSC3DforSubspace(int subspaceId); // Return if the update is successful
	const std::vector<FLOATVECTOR3>& getScatterplot3DSamples() const;
	const std::vector<FLOATVECTOR3>& getScatterplot3DSampleColors() const;

	QImage&  dispBuffer(DISPBUF_MODE dispMode);  // retrieve color buffer
	cv::Mat& accumBuffer(ACCUMBUF_MODE accumBufMode, int bufId = 0); // retrieve the 64-bit accumulation buffer
	void     writeHdrAccumBufFile(ACCUMBUF_MODE bufId);
	bool     savePCP1LayerToHdrFile(PC_RENDER_MODE pcMode, const QString& fileName);
	bool     savePCP2LayerToHdrFile(PC_RENDER_MODE pcMode, const QString& fileName);
	bool     colorizeLayerToHdrFile(ACCUMBUF_MODE accumBufMode);
	void     colormapAccumBuf(cv::Mat& colorBuf, const cv::Mat& accumBuf, int encodeChannel); // colorize 2-channel cv::Mat to 3-channel cv::Mat based on the "encodeChannel" content
	void     colormapDensityAccumBuf(cv::Mat& colorBuf, const cv::Mat& densityBuf, int subdimId); // colorize the density buf based on the density!

	void     clearAccumBuffer(int layer, int plotMode = 0);    // Clear accumulation buffer (0-, 1-, 2-flats)
	QImage&  createDispImageFromAccumBuffer(int plotMode = 0); // plotMode == 0: pcpPlot, plotMode == 1: splomPlot
	QImage&  createFocusContextDispImageFromAccumBuffer(int firstFocusAxis, int focusDim, int plotMode = 0); 
	QImage&  createFocusContextPCPDispImageRedraw(int firstFocusAxis, int focusDim); // Redraw polylines/or transformed polylines for PCP with given sub dimensions
	QImage&  createDispImageFrom1FlatAccumBuffer(int plotMode = 0);
	QImage&  createDispImageFrom2FlatAccumBuffer();
	//-------------------------------------------------
	// Nguyen-Rosen paper: needs 2 buffers for 1-flats
	QImage&  createDispImagePosCorr1FlatAccumBuffer();
	QImage&  createDispImageNegCorr1FlatAccumBuffer();
	//-------------------------------------------------
	QImage&  createDispImagePCPLineDensityBuilder(); // returns the _densPCP_compDispBuf
	QImage&  createFocusContextDispImagePCPLineDensityBuilder(int firstFocusAxis, int focusDim); // returns the focus/context version of _densPCP_compDispBuf
	void     redrawOneDispImagePCPLineDensityBuilder(int plotId); // redraw only one density plot
	void     setPCPLineDensityBuilderTF(int tfId); // set tf 

	void resetSampleCnt() { _sampleCnt = 0; resetPlot(); }
	void drawQueryResults();
	// Conversions between coordinates
	void Norm2InselbergCoords(float nx, float ny, float& x, float& y);
	void InselbergCoords2InselNorm(float x, float y, float& nx, float& ny);

	void DispBuffer2InselbergCoords(INTVECTOR2 scrPt, FLOATVECTOR2& pcpPt);
	void InselbergCoords2DispBuffer(FLOATVECTOR2 pcpPt, INTVECTOR2& scrPt);
	// PC line/curve 
	void drawPCline(const std::vector<float>& s, int subspaceId, int subspaceDim, std::vector<FLOATVECTOR2>& pcline, const XFORM_OPTION& xformMode);
	// Find regions of interest the density plot!
	void findDensityPlotROI();
	cv::Mat findDensityPlotROI(int dimId);
	void updateDensityMapRenderMasks();
protected:
	bool updatePlotFromLayerData(); 
	bool updatePlotFullFromLayerData();

protected:
	UINT64							_numTotalSamples;
	UINT64							_sampleCnt;
	int                             _valueDim; // dimension of the value in the data
	// Plot builders
	PlotBuilder*                   _pcpPlotBuilder;	  // the pcp plot
	PlotBuilder*                   _scPlotBuilder;	  // the 2D scatterplot
	PlotBuilder*                   _sc3DPlotBuilder;  // the 3D scatterplot
	PlotBuilder*                   _splomPlotBuilder; // the scatterplot matrix plot 
	// Data of XPCP samples
	std::vector<XPCPSample>        _xpcpData;
};

