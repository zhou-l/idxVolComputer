#include "PointDataHandler.h"
#include "PCPbuilder.h"
#include "PCPLineDensityBuilder.h"
#include "SplomPlotBuilder.h"
#include "ScatterPlot3DBuilder.h"
#include <QPainter>
#include  "PCimageProcessor.h"
#include <QMessageBox>
using namespace std;
#define IDX_PT_DEBUG 1
PointDataHandler::PointDataHandler():
_pcpPlotBuilder(NULL),
_scPlotBuilder(NULL),
_splomPlotBuilder(NULL),
_sc3DPlotBuilder(NULL),
_sampleCnt(0),
_valueDim(0)
{}

PointDataHandler::PointDataHandler(UINT64 numSamples, int dataDim, int plotDispBufW, int plotDispBufH, const vector<XPCPSample>& pointData)
{
	_xpcpData = pointData;
	_numTotalSamples = numSamples;
	_valueDim = dataDim;
	_sampleCnt = 0;

	////////////////////////////////////////////////////////////
	// PCP plot 
	switch (g_params.Pcp_RenderMode())
	{
	case PR_LINE_IDX_PT:
		_pcpPlotBuilder = new PCPbuilder(plotDispBufW, plotDispBufH);
		// Shall we flip Y????
		dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->flipY(true);
		break;
	case PR_DS_LINE_DENSITY:
		_pcpPlotBuilder = new PCPLineDensityBuilder(plotDispBufW, plotDispBufH);
		// Shall we flip Y????
		dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->flipY(true);
		cout << "PCPLineDensityBuilder created!" << endl;
		break;
	default:
		_pcpPlotBuilder = new PCPbuilder(plotDispBufW, plotDispBufH);
		// Shall we flip Y????
		dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->flipY(true);
		break;
	}
	/////////////////////////////////////

	// SC plot
	_scPlotBuilder = NULL;

	// Splom plot
	_splomPlotBuilder = new SplomPlotBuilder(g_splomDispBufWidth, g_splomDispBufHeight);

	// 3D scatterplot. For now, we don't use it.
	_sc3DPlotBuilder = NULL;
	//if (_valueDim == 3)
	int dataRange = 128;
	_sc3DPlotBuilder = new ScatterPlot3DBuilder(dataRange, dataRange, dataRange, dataRange, dataRange, dataRange);
}


PointDataHandler::~PointDataHandler()
{
	SAFE_DELETE(_pcpPlotBuilder);
	SAFE_DELETE(_scPlotBuilder);
	SAFE_DELETE(_sc3DPlotBuilder);
	SAFE_DELETE(_splomPlotBuilder);
}

void PointDataHandler::resetPlot()
{
	_sampleCnt = 0;
	if (g_params.Pcp_RenderMode() == PR_LINE_IDX_PT)
	{
		_pcpPlotBuilder->clearAccumBuffer(0);
		_pcpPlotBuilder->clearAccumBuffer(1);
		_pcpPlotBuilder->clearAccumBuffer(2);
	}
	else
	{
		_pcpPlotBuilder->clearAccumBuffer(0);
		dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->clearAllDensityPlotAccumBuffers();
	}

	_splomPlotBuilder->clearAccumBuffer(0);
	_splomPlotBuilder->clearAccumBuffer(1);
}

bool PointDataHandler::updatePlot()
{
	if (!g_sampled_raw.empty())
		return updatePlotFromLayerData();
	else
		return false;
}

bool PointDataHandler::updatePlotFull()
{
	if (!g_sampled_raw.empty())
		return updatePlotFullFromLayerData();
	else
		return false;
}


// The PC layer drawing function
void PointDataHandler::updatePlotLayer(int layer)
{
	// For density map mode.
	if (g_params.IsIdxPtDensityMap())
	{
		dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->drawDensity_pFlat_layer(layer);
		return;
	}
	//// Always draw all layers now!
	//if (!g_sampled_raw.empty())
	//{
	//	if (g_params.FirstAxis_subspace() == -1)
	//	{
	//		cout << "redraw background with full opacity!" << endl;
	//		for (size_t i = 0; i < g_sampled_raw.size(); i++)
	//			dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->drawPCPlayerSample(g_sampled_raw[i], 0);
	//	}
	//	else
	//	{
	//		int subspaceNum = g_params.Pcp_selectedLayer() == 1 ? 2 : 3;
	//		for (size_t i = 0; i < g_sampled_raw.size(); i++)
	//		{
	//			dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->drawPCP_basic_layersample(g_sampled_raw[i], g_params.FirstAxis_subspace(), subspaceNum);
	//		}
	//	}
	//	
	//}

	///////////////
	// For the regular mode with point-based indexed points
	//if (layer == 1 || layer == 3)
	{
		cout << "size of 1flat array " << g_1flat_list.size() << endl;
		if (!g_1flat_list.empty())
		{
#ifndef IDX_SCRPT_DEBUG
			for (size_t i = 0; i < g_1flat_list.size(); i++)
			{
					dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->drawPCP_pFlat_layersample(g_1flat_list[i], 1);
			}
#endif

#ifdef IDX_SCRPT_DEBUG
			ofstream ofScrPtDebug("ImgPlaneIdxPtDebug_compute.txt");
			ofScrPtDebug << "normalized Pt; screen Pt" << endl;

			for (size_t i = 0; i < g_1flat_list.size(); i++)
			{
				FLOATVECTOR2 normPt;
				INTVECTOR2   scrnPt;

				dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->drawPCP_pFlat_layersample_debug(g_1flat_list[i], 1, normPt, scrnPt);
				ofScrPtDebug << normPt << "; " << scrnPt << endl;
			}
			ofScrPtDebug.close();
#endif
	/*		return false;*/
		}


		cout << "size of 1flat array S2 " << g_1flat_listS2.size() << endl;
		if (g_params.Pcp_compare_mode() > 0) 
		{
			cout << "Draw 1flat S2 array with diff compare mode." << endl;
			if (!g_1flat_listS2.empty())
			{
				for (size_t i = 0; i < g_1flat_listS2.size(); i++)
				{
					dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->drawPCP_pFlat_layersample(g_1flat_listS2[i], 1, true, true);
				}
				/*		return false;*/
			}

		}
	}


	//if (layer == 2 || layer == 3) // Third pass
	{
		if (!g_2flat_list.empty())
		{
			for (size_t i = 0; i < g_2flat_list.size(); i++)
			{
				//for (int spaceId = 0; spaceId < numTwoFlatPerSample; spaceId++)
					dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->drawPCP_pFlat_layersample(g_2flat_list[i], 2);
			}
		}
	}
}

bool PointDataHandler::updateSC3DforSubspace(int subspaceId)
{
	if (_sc3DPlotBuilder == NULL)
		return false;
	int ssid = (subspaceId == 0) ? 0 : subspaceId - 1;
	if (g_sampled_raw.size() == 0)
		return false;
	else
	{
		if (g_sampled_raw[0].size() - 3 < ssid)
			return false;
	}
	dynamic_cast<ScatterPlot3DBuilder*>(_sc3DPlotBuilder)->clearSC3DBuffers();
	for (size_t cnt = 0; cnt < g_sampled_raw.size(); cnt++)
	{
		dynamic_cast<ScatterPlot3DBuilder*>(_sc3DPlotBuilder)->processSubspaceSample(g_sampled_raw[cnt], ssid);
	}
	dynamic_cast<ScatterPlot3DBuilder*>(_sc3DPlotBuilder)->drawQueryResults();
	return true;
}


const vector<FLOATVECTOR3>& PointDataHandler::getScatterplot3DSamples() const
{
	return dynamic_cast<ScatterPlot3DBuilder*>(_sc3DPlotBuilder)->getSamples();
}

const vector<FLOATVECTOR3>& PointDataHandler::getScatterplot3DSampleColors() const
{
	return dynamic_cast<ScatterPlot3DBuilder*>(_sc3DPlotBuilder)->getSampleColors();
}

// 
bool PointDataHandler::updatePlotFromLayerData()// update plot from already computed layer data. NOTE: this is where the first pass of plot drawing done!!!
{
	if (_numTotalSamples != g_sampled_raw.size())
		_numTotalSamples = g_sampled_raw.size();
	// Just do one random sample
	int pass = 1;
	if (g_params.Pcp_RenderMode() != PR_DS_LINE_DENSITY)
	{
		if (!g_1flat_list.empty())
			pass++;
		if (!g_2flat_list.empty())
			pass++;
	}
	
	if (_sampleCnt >= pass * _numTotalSamples)
	{
		if (g_params.Pcp_RenderMode() == PR_DS_LINE_DENSITY){
			
			// write out density buffer
			cv::Mat& buf = dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->densityPlotAccumBuffer(0);
			//PCimageProcessor::smoothAccumBuf(buf);
			//cv::Mat maxPosBuf = PCimageProcessor::getMaximalLocs(buf);
			////cv::Mat gradBuf = PCimageProcessor::localShapeFilter(buf);
			//g_params.DensMapROILocs(maxPosBuf);
		}
		return true; // we are done here!
	}
	int startCnt = _sampleCnt;
	const int num_sample_per_update = g_num_sample_per_update;
	int numOneFlatPerSample = _valueDim - 1;
	int numTwoFlatPerSample = _valueDim - 2;

	// For PCP downsampled data lines + full density mode: only one pass!
	if (g_params.Pcp_RenderMode() == PR_DS_LINE_DENSITY)
	{
		cout << "Draw With DS_LINE_DENSITY" << endl;
		if (_sampleCnt < _numTotalSamples) // First pass
		{
			if (!g_sampled_raw.empty()) // If we have _xpcpData from someone, only draw from the _xpcpData
			{
				for (; _sampleCnt < MIN(_numTotalSamples, startCnt + num_sample_per_update); _sampleCnt++)
				{
					dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->processSample(g_sampled_raw[_sampleCnt]);

					bool drawSample = false;
					if (_sampleCnt % g_pcp_num_samples_skip_DS_LINE_MODE == 0)
						drawSample = true; // draw every g_pcp_num_samples_skip_DS_LINE_MODE samples
					if (std::find(g_params.OutlierIds().begin(), g_params.OutlierIds().end(), _sampleCnt) != g_params.OutlierIds().end())
						drawSample = true; // outliers have to be drawn!

					if (drawSample) // Also draw data entry from the downsampled version of data
						dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->processOverviewSample(g_sampled_raw[_sampleCnt]);
					// Update other views
					if (_splomPlotBuilder)
						dynamic_cast<SplomPlotBuilder*>(_splomPlotBuilder)->processSample(g_sampled_raw[_sampleCnt]);
					if (_sc3DPlotBuilder) // For 3D scatterplot, plot only if we have three attributes, and select a subspace
					{
						if (g_attrib_names.size() >= 3)
							dynamic_cast<ScatterPlot3DBuilder*>(_sc3DPlotBuilder)->processSample(g_sampled_raw[_sampleCnt]);
					}
				}
				return false;
			}
		}
	}
	else
	{
		//-------------------------------------------------------------------------------------
		// pcp + p-flat indexed points mode:
		if (_sampleCnt < _numTotalSamples) // First pass
		{
			if (!g_sampled_raw.empty()) // If we have _xpcpData from someone, only draw from the _xpcpData
			{
				for (; _sampleCnt < MIN(_numTotalSamples, startCnt + num_sample_per_update); _sampleCnt++)
				{
					if (_sampleCnt % g_pcp_num_samples_skip_LINE_IDXPT_MODE == 0) // Also draw data entry from the downsampled version of data
					{
						if(!g_params.doClustering())
							dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->drawPCPlayerSample(g_sampled_raw[_sampleCnt], 0);
						else
						{
							dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->drawPCPlayerSample(g_sampled_raw[_sampleCnt], g_params.clusterLabelData()[_sampleCnt], 0);
							
						}
					}

					//dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->drawPCPlayerSample(g_sampled_raw[_sampleCnt], 0);
					// Update other views
					if (_splomPlotBuilder)
						dynamic_cast<SplomPlotBuilder*>(_splomPlotBuilder)->processSample(g_sampled_raw[_sampleCnt]);
					if (_sc3DPlotBuilder) // For 3D scatterplot, plot only if we have three attributes, and select a subspace
					{
						if (g_attrib_names.size() >= 3)
							dynamic_cast<ScatterPlot3DBuilder*>(_sc3DPlotBuilder)->processSample(g_sampled_raw[_sampleCnt]);
					}
				}
				return false;
			}
		}

		if (_sampleCnt >= _numTotalSamples && _sampleCnt < _numTotalSamples * 2) // Second pass
		{
			cout << "size of 1flat array " << g_1flat_list.size() << endl;
			UINT64 splomSampleCnt = _sampleCnt;
			if (!g_1flat_list.empty())
			{

				for (; _sampleCnt < MIN(2 * _numTotalSamples, startCnt + 2 * num_sample_per_update); _sampleCnt++)
				{
					for (int spaceId = 0; spaceId < numOneFlatPerSample; spaceId++){
						dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->drawPCP_pFlat_layersample(g_1flat_list[UINT64(_sampleCnt - _numTotalSamples) * UINT64(numOneFlatPerSample) + UINT64(spaceId)], 1);
					}
				}

			}

			cout << "size of 1flat array S2" << g_1flat_listS2.size() << endl;
			if (g_params.Pcp_compare_mode() > 0) 
			{
				cout << "draw 1flat array S2 with diff compare mode." << endl;
				if (!g_1flat_listS2.empty())
				{
					UINT64 scnt = splomSampleCnt;
					for (; scnt < MIN(2 * _numTotalSamples, startCnt + 2 * num_sample_per_update); scnt++)
					{
						for (int spaceId = 0; spaceId < numOneFlatPerSample; spaceId++) {
							dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->drawPCP_pFlat_layersample(g_1flat_listS2[UINT64(scnt - _numTotalSamples) * UINT64(numOneFlatPerSample) + UINT64(spaceId)], 1, true, true);
						}
					}
				}
			}
		
			// Paint on splom using major eigen vector
			if (!g_majEigData.empty())
			{
				for (; splomSampleCnt < MIN(2 * _numTotalSamples, startCnt + 2 * num_sample_per_update); splomSampleCnt += g_trendData_stride) // Draw g_trendData_stride samples
				{
					UINT64 ii = splomSampleCnt - _numTotalSamples;
					dynamic_cast<SplomPlotBuilder*>(_splomPlotBuilder)->draw1DTrend(g_sampled_raw[ii], g_majEigData[ii]);
				}
			}
			return false;
		}

		if (_sampleCnt >= 2 * _numTotalSamples) // Third pass
		{
			if (!g_2flat_list.empty())
			{
				UINT64 twoFlatCnt = UINT64(_sampleCnt - 2 * _numTotalSamples);
				UINT64 drawCntThisUpdate = 0;
				UINT64 totalTwoFlats = UINT64(g_2flat_list.size()) / UINT64(numTwoFlatPerSample);
				for (; twoFlatCnt < totalTwoFlats && drawCntThisUpdate < num_sample_per_update; twoFlatCnt++)
				{
					for (int spaceId = 0; spaceId < numTwoFlatPerSample; spaceId++)
						dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->drawPCP_pFlat_layersample(g_2flat_list[twoFlatCnt * UINT64(numTwoFlatPerSample) + UINT64(spaceId)], 2);
					drawCntThisUpdate++;
				}
				_sampleCnt += twoFlatCnt;
				if (twoFlatCnt >= g_2flat_list.size())
				{
					_sampleCnt = pass * _numTotalSamples;
					return true; // The number of 2-flats went wrong, but we are done here anyways! 
				}
				return false;
			}
		}
	}
	
	return false;
}

bool PointDataHandler::updatePlotFullFromLayerData()
{
	if (_numTotalSamples != g_sampled_raw.size())
		_numTotalSamples = g_sampled_raw.size();
	// Just do one random sample
	int pass = 1;
	if (g_params.Pcp_RenderMode() != PR_DS_LINE_DENSITY)
	{
		if (!g_1flat_list.empty())
			pass++;
		if (!g_2flat_list.empty())
			pass++;
	}

	if (_sampleCnt >= pass * _numTotalSamples)
	{
		if (g_params.Pcp_RenderMode() == PR_DS_LINE_DENSITY) {

			// write out density buffer
			cv::Mat& buf = dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->densityPlotAccumBuffer(0);
			//PCimageProcessor::smoothAccumBuf(buf);
			//cv::Mat maxPosBuf = PCimageProcessor::getMaximalLocs(buf);
			////cv::Mat gradBuf = PCimageProcessor::localShapeFilter(buf);
			//g_params.DensMapROILocs(maxPosBuf);
		}
		return true; // we are done here!
	}
	int startCnt = _sampleCnt;
	const int num_sample_per_update = g_num_sample_per_update;
	int numOneFlatPerSample = _valueDim - 1;
	int numTwoFlatPerSample = _valueDim - 2;

	// For PCP downsampled data lines + full density mode: only one pass!
	if (g_params.Pcp_RenderMode() == PR_DS_LINE_DENSITY)
	{
		cout << "Draw With DS_LINE_DENSITY" << endl;
		if (_sampleCnt < _numTotalSamples) // First pass
		{
			if (!g_sampled_raw.empty()) // If we have _xpcpData from someone, only draw from the _xpcpData
			{
				for (size_t i = 0; i < _numTotalSamples; i++)
				{
					dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->processSample(g_sampled_raw[i]);

					bool drawSample = false;
					if (i % g_pcp_num_samples_skip_DS_LINE_MODE == 0)
						drawSample = true; // draw every g_pcp_num_samples_skip_DS_LINE_MODE samples
					if (std::find(g_params.OutlierIds().begin(), g_params.OutlierIds().end(), i) != g_params.OutlierIds().end())
						drawSample = true; // outliers have to be drawn!

					if (drawSample) // Also draw data entry from the downsampled version of data
						dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->processOverviewSample(g_sampled_raw[i]);
					// Update other views
					if (_splomPlotBuilder)
						dynamic_cast<SplomPlotBuilder*>(_splomPlotBuilder)->processSample(g_sampled_raw[i]);
					if (_sc3DPlotBuilder) // For 3D scatterplot, plot only if we have three attributes, and select a subspace
					{
						if (g_attrib_names.size() >= 3)
							dynamic_cast<ScatterPlot3DBuilder*>(_sc3DPlotBuilder)->processSample(g_sampled_raw[i]);
					}
				}
				_sampleCnt = _numTotalSamples;
				return false;
			}
		}
	}
	else
	{
		//-------------------------------------------------------------------------------------
		// pcp + p-flat indexed points mode:
		if (_sampleCnt < _numTotalSamples) // First pass
		{
			if (!g_sampled_raw.empty()) // If we have _xpcpData from someone, only draw from the _xpcpData
			{

				for (; _sampleCnt < _numTotalSamples; _sampleCnt++)
				{
					if (_sampleCnt % g_pcp_num_samples_skip_LINE_IDXPT_MODE == 0) // Also draw data entry from the downsampled version of data
					{
						if (!g_params.doClustering())
							dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->drawPCPlayerSample(g_sampled_raw[_sampleCnt], 0);
						else
						{
							dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->drawPCPlayerSample(g_sampled_raw[_sampleCnt], g_params.clusterLabelData()[_sampleCnt], 0);
						}
					}

					//dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->drawPCPlayerSample(g_sampled_raw[_sampleCnt], 0);
					// Update other views
					if (_splomPlotBuilder)
						dynamic_cast<SplomPlotBuilder*>(_splomPlotBuilder)->processSample(g_sampled_raw[_sampleCnt]);
					if (_sc3DPlotBuilder) // For 3D scatterplot, plot only if we have three attributes, and select a subspace
					{
						if (g_attrib_names.size() >= 3)
							dynamic_cast<ScatterPlot3DBuilder*>(_sc3DPlotBuilder)->processSample(g_sampled_raw[_sampleCnt]);
					}
				}
				return false;
			}
		}

		if (_sampleCnt >= _numTotalSamples && _sampleCnt < _numTotalSamples * 2) // Second pass
		{
			cout << "size of 1flat array " << g_1flat_list.size() << endl;
			UINT64 splomSampleCnt = _sampleCnt;
			if (!g_1flat_list.empty())
			{
#ifdef IDX_SCRPT_DEBUG
				ofstream ofScrPtDebug("ImgPlaneIdxPtDebug_load.txt");
				ofScrPtDebug << "normalized Pt; screen Pt" << endl;
#endif
				for (; _sampleCnt < 2 * _numTotalSamples; _sampleCnt++)
				{
					for (int spaceId = 0; spaceId < numOneFlatPerSample; spaceId++) {
						FLOATVECTOR2 normPt;
						INTVECTOR2   scrnPt;
						dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->drawPCP_pFlat_layersample_debug(g_1flat_list[UINT64(_sampleCnt - _numTotalSamples) * UINT64(numOneFlatPerSample) + UINT64(spaceId)], 1, normPt, scrnPt);
#ifdef IDX_SCRPT_DEBUG
						ofScrPtDebug << normPt << "; " << scrnPt << endl;
#endif
					}
				}
#ifdef IDX_SCRPT_DEBUG
				ofScrPtDebug.close();
#endif
			}

			cout << "size of 1flat array S2" << g_1flat_listS2.size() << endl;
			if (g_params.Pcp_compare_mode() > 0)
			{
				cout << "draw 1flat array S2 with diff compare mode." << endl;
				if (!g_1flat_listS2.empty())
				{
					UINT64 scnt = splomSampleCnt;
					for (; scnt < 2 * _numTotalSamples; scnt++)
					{
						for (int spaceId = 0; spaceId < numOneFlatPerSample; spaceId++) {
							dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->drawPCP_pFlat_layersample(g_1flat_listS2[UINT64(scnt - _numTotalSamples) * UINT64(numOneFlatPerSample) + UINT64(spaceId)], 1, true, true);
						}
					}
				}
			}

			// Paint on splom using major eigen vector
			if (!g_majEigData.empty())
			{
				for (; splomSampleCnt < 2 * _numTotalSamples; splomSampleCnt += g_trendData_stride) // Draw g_trendData_stride samples
				{
					UINT64 ii = splomSampleCnt - _numTotalSamples;
					dynamic_cast<SplomPlotBuilder*>(_splomPlotBuilder)->draw1DTrend(g_sampled_raw[ii], g_majEigData[ii]);
				}
			}
			return false;
		}

		if (_sampleCnt >= 2 * _numTotalSamples) // Third pass
		{
			if (!g_2flat_list.empty())
			{
				UINT64 twoFlatCnt = UINT64(_sampleCnt - 2 * _numTotalSamples);
				UINT64 drawCntThisUpdate = 0;
				UINT64 totalTwoFlats = UINT64(g_2flat_list.size()) / UINT64(numTwoFlatPerSample);
				for (; twoFlatCnt < totalTwoFlats; twoFlatCnt++)
				{
					for (int spaceId = 0; spaceId < numTwoFlatPerSample; spaceId++)
						dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->drawPCP_pFlat_layersample(g_2flat_list[twoFlatCnt * UINT64(numTwoFlatPerSample) + UINT64(spaceId)], 2);
					drawCntThisUpdate++;
				}
				_sampleCnt += twoFlatCnt;
				if (twoFlatCnt >= g_2flat_list.size())
				{
					_sampleCnt = pass * _numTotalSamples;
					return true; // The number of 2-flats went wrong, but we are done here anyways! 
				}
				return false;
			}
		}
	}

	return false;
}

void PointDataHandler::drawPCline(const vector<float>& s, int subspaceId, int subspaceDim, vector<FLOATVECTOR2>& pcline, const XFORM_OPTION& xformMode)
{
	bool fullLine = true;
	dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->pcLine(s, subspaceId, subspaceDim, xformMode, fullLine, pcline);
}

void PointDataHandler::drawQueryResults()
{
	cout << "pointDataHandler: Draw query result" << endl;
	// Ask plot builders to draw query results
	if (_pcpPlotBuilder)
		_pcpPlotBuilder->drawQueryResults();
	if (_splomPlotBuilder)
		_splomPlotBuilder->drawQueryResults();
	if (_sc3DPlotBuilder && _valueDim >= 3){
		_sc3DPlotBuilder->drawQueryResults();
		
	}
}


QImage& PointDataHandler::dispBuffer(DISPBUF_MODE dispMode)
{
	switch (dispMode)
	{
	case DISPBUF_PCP_PLOT:
		return _pcpPlotBuilder->dispBuffer();
	case DISPBUF_PCP_1FLAT:
	case DISPBUF_PCP_NEG_1FLAT:
		return dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->oneFlatBuffer();
	case DISPBUF_PCP_2FLAT:
		return dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->twoFlatBuffer();
	case DISPBUF_SCPLOT:
		return _scPlotBuilder->dispBuffer();
	case DISPBUF_SPLOMPLOT:
		return _splomPlotBuilder->dispBuffer();
	case DISPBUF_SPLOM_1FLAT:
		return dynamic_cast<SplomPlotBuilder*>(_splomPlotBuilder)->oneDTrendBuffer();
	case DISPBUF_PCP_PLOT_SEL:
		return _pcpPlotBuilder->selBuffer();
	case DISPBUF_SCPLOT_SEL:
		return _scPlotBuilder->selBuffer();
	case DISPBUF_SPLOMPLOT_SEL:
		return _splomPlotBuilder->selBuffer();
	case DISPBUF_PCP_DENSITY:
		return dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->densityPlotBuffer();
	case DISPBUF_PCP_POS_1FLAT:
		return dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->focusPCPbuffer();
	default:
		return _pcpPlotBuilder->dispBuffer();
	}
}

cv::Mat& PointDataHandler::accumBuffer(ACCUMBUF_MODE accumBufMode, int bufId)
{
	switch (accumBufMode)
	{
	case ACCUMBUF_PCP_BASE:
		return dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->accumBuffer();
		break;
	case ACCUMBUF_PCP_1LAYER:
		if (g_params.Pcp_RenderMode() == PR_LINE_IDX_PT)
			return dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->accumBuffer1Flat();
		else
			return dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->densityPlotAccumBuffer(bufId);
		break;
	case ACCUMBUF_PCP_2LAYER:
		if (g_params.Pcp_RenderMode() == PR_LINE_IDX_PT)
			return dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->accumBuffer2Flat();
		else
			return dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->densityPlotAccumBuffer(bufId);
		break;
	case ACCUMBUF_SPLOM_BASE:
		return dynamic_cast<SplomPlotBuilder*>(_splomPlotBuilder)->accumBuffer();
		break;
	case ACCUMBUF_SPLOM_1LAYER:
		return dynamic_cast<SplomPlotBuilder*>(_splomPlotBuilder)->accumBuffer1DTrend();
		break;

	default:
		return _pcpPlotBuilder->accumBuffer();
	}
	return _pcpPlotBuilder->accumBuffer();
}

void PointDataHandler::writeHdrAccumBufFile(ACCUMBUF_MODE bufId)
{
	QString fileName;
	bool isOK = false;
	switch(bufId)
	{
	case ACCUMBUF_PCP_BASE:
		fileName = QString("accum_buf_pcp_base.hdr");
		isOK = PlotBuilder::saveAccumBufferToHdrFile(_pcpPlotBuilder->accumBuffer(), fileName);
		break;
	case ACCUMBUF_PCP_1LAYER:
		fileName = QString("accum_buf_pcp_layer1.hdr");
		isOK = savePCP1LayerToHdrFile(g_params.Pcp_RenderMode(), fileName);
		break;
	case ACCUMBUF_PCP_2LAYER:
		fileName = QString("accum_buf_pcp_layer2.hdr");
		isOK = savePCP2LayerToHdrFile(g_params.Pcp_RenderMode(), fileName);
		break;
	case ACCUMBUF_SPLOM_BASE:
		fileName = QString("accum_buf_splom_base.hdr");
		isOK = PlotBuilder::saveAccumBufferToHdrFile(dynamic_cast<SplomPlotBuilder*>(_splomPlotBuilder)->accumBuffer(), fileName);
		break;
	case ACCUMBUF_SPLOM_1LAYER:
		fileName = QString("accum_buf_splom_layer1.hdr");
		isOK = PlotBuilder::saveAccumBufferToHdrFile(dynamic_cast<SplomPlotBuilder*>(_splomPlotBuilder)->accumBuffer1DTrend(), fileName);
		break;
	}
}

bool PointDataHandler::savePCP1LayerToHdrFile(PC_RENDER_MODE pcmode, const QString& fileName)
{
	bool isOK = false;
	switch (pcmode)
	{
	case PR_LINE_IDX_PT:
		isOK = PlotBuilder::saveAccumBufferToHdrFile(dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->accumBuffer1Flat(), fileName);
		break;
	case PR_DS_LINE_DENSITY:
		// TODO: enable buffer id selection!
		// Currently, always use the first buffer!
		isOK = PlotBuilder::saveAccumBufferToHdrFile(dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->densityPlotAccumBuffer(0), fileName);
		break;
	default:
		isOK = PlotBuilder::saveAccumBufferToHdrFile(dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->accumBuffer1Flat(), fileName);
		break;
	}
	return isOK;
}

bool PointDataHandler::savePCP2LayerToHdrFile(PC_RENDER_MODE pcmode, const QString& fileName)
{
	bool isOK = false;
	switch (pcmode)
	{
	case PR_LINE_IDX_PT:
		isOK = PlotBuilder::saveAccumBufferToHdrFile(dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->accumBuffer2Flat(), fileName);
		break;
	case PR_DS_LINE_DENSITY:
		// TODO: enable buffer id selection!
		// Currently, always use the first buffer!
		isOK = PlotBuilder::saveAccumBufferToHdrFile(dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->densityPlotAccumBuffer(0), fileName);
		break;
	default:
		isOK = PlotBuilder::saveAccumBufferToHdrFile(dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->accumBuffer2Flat(), fileName);
		break;
	}
	return isOK;
}

bool PointDataHandler::colorizeLayerToHdrFile(ACCUMBUF_MODE bufMode)
{
	cv::Mat colorBuf;
	QString fileName;
	bool isOK = false;
	if (g_params.Pcp_RenderMode() == PR_LINE_IDX_PT)
	{
		cv::Mat& accumBuf = accumBuffer(bufMode);
		colorBuf = cv::Mat::zeros(accumBuf.rows, accumBuf.cols, CV_32FC3);
		// colormapping accum buffer
		colormapAccumBuf(colorBuf, accumBuf, 1);
	
		switch (bufMode)
		{
		case ACCUMBUF_PCP_1LAYER:
			fileName = QString("color_buf_pcp_layer1.hdr");
			break;
		case ACCUMBUF_PCP_2LAYER:
			fileName = QString("color_buf_pcp_layer2.hdr");
			break;
		case ACCUMBUF_SPLOM_BASE:
			fileName = QString("color_buf_splom_base.hdr");
			break;
		case ACCUMBUF_SPLOM_1LAYER:
			fileName = QString("color_buf_splom_layer1.hdr");
			break;
		}
		isOK = PlotBuilder::saveAccumBufferToHdrFile(colorBuf, fileName);
	}
	else if (g_params.Pcp_RenderMode() == PR_DS_LINE_DENSITY)
	{
		if (bufMode != ACCUMBUF_PCP_1LAYER)
		{
			cout << "Something went wrong, only the density layer can be rendered!" << endl;
			return isOK;
		}
		int num_bufs = dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->numDensityPlotAccumBuf();
		// the buffer can only be the density plot!
		for (int i = 0; i < num_bufs; i++)
		{
			cv::Mat& accumBuf = accumBuffer(bufMode, i);
			colorBuf = cv::Mat::zeros(accumBuf.rows, accumBuf.cols, CV_32FC3);
			// colormapping accum buffer
			colormapDensityAccumBuf(colorBuf, accumBuf, i);

			switch (bufMode)
			{
			case ACCUMBUF_PCP_1LAYER:
				fileName = QString("color_buf_pcp_density%1.hdr").arg(i);
				break;
			}
			isOK = PlotBuilder::saveAccumBufferToHdrFile(colorBuf, fileName);
			if (!isOK)
				return isOK;
		}
	}
	else
	{
		cout << "Render mode not supported yet!" << endl;
		return isOK;
	}
	return isOK;
}

void PointDataHandler::colormapAccumBuf(cv::Mat& colorBuf, const cv::Mat& accumBuf, int encodeChannel)
{
	float uc2fl = 1.0f / 255.0f;
	if (accumBuf.channels() == 2 && encodeChannel < accumBuf.channels())
	{
		float miny = numeric_limits<float>::max();
		float maxy = -numeric_limits<float>::max();
		cout << "Colormap the the buffer to HDR file." << endl;
		for (int r = 0; r < accumBuf.rows; r++)
		{
			for (int c = 0; c < accumBuf.cols; c++)
			{
				double y = accumBuf.at<cv::Vec2d>(r, c)[0];
				maxy = MAX(y, maxy);
				miny = MIN(y, miny);

				int sid = int(accumBuf.at<cv::Vec2d>(r, c)[encodeChannel]);
				FLOATVECTOR3 rgbCol;
				INTVECTOR4 iRGBCol = g_pcp_1flat_colormap[sid];
				rgbCol = FLOATVECTOR3(float(iRGBCol[0]) * uc2fl, float(iRGBCol[1])* uc2fl, float(iRGBCol[2]) * uc2fl);
				rgbCol = rgbCol * y; //scale by the value!
				colorBuf.at<cv::Vec3f>(r,c) = cv::Vec3f(rgbCol.x, rgbCol.y, rgbCol.z);
			}
		}
		cout << "Y min = " << miny << ", max = " << maxy << endl;
	}
}

void PointDataHandler::colormapDensityAccumBuf(cv::Mat& colorBuf, const cv::Mat& densityBuf, int subdimId)
{
	float uc2fl = 1.0f / 255.0f;
	{
		float miny = numeric_limits<float>::max();
		float maxy = -numeric_limits<float>::max();
		cout << "Colormap the the buffer to HDR file." << endl;
		for (int r = 0; r < densityBuf.rows; r++)
		{
			for (int c = 0; c < densityBuf.cols; c++)
			{
				double y = densityBuf.at<double>(r, c);
				maxy = MAX(y, maxy);
				miny = MIN(y, miny);

				FLOATVECTOR3 rgbCol;
				INTVECTOR4 iRGBCol = g_pcp_1flat_colormap[subdimId];
				rgbCol = FLOATVECTOR3(float(iRGBCol[0]) * uc2fl, float(iRGBCol[1])* uc2fl, float(iRGBCol[2]) * uc2fl);
				rgbCol = rgbCol * y; //scale by the value!
				colorBuf.at<cv::Vec3f>(r, c) = cv::Vec3f(rgbCol.x, rgbCol.y, rgbCol.z);
			}
		}
		cout << "Y min = " << miny << ", max = " << maxy << endl;
	}
}

void PointDataHandler::clearAccumBuffer(int layer, int plotMode)
{
	if (plotMode == 0)
	{
		dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->clearAccumBuffer(layer);
	}
	else
		dynamic_cast<SplomPlotBuilder*>(_splomPlotBuilder)->clearAccumBuffer(layer);
}

QImage& PointDataHandler::createDispImagePCPLineDensityBuilder()
{
	dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->renderAllDispBuffers();
	//cout << "Create disp Image with PCPLineDense" << endl;
	return dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->compositeMainDispBuffer();
}

void PointDataHandler::redrawOneDispImagePCPLineDensityBuilder(int plotId)
{
	//----------------------------------------------
	// Densitymap Rendering Methods:
	// global tf
	if (!g_params.IsDrawMaxLocs())
		dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->renderToDispBuffer(plotId);
	// try local tf
	//dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->renderToDispBuffer_localTF(plotId);
	// try local tf with maximal locations
	else
		dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->renderToDispBuffer_max_filt_localTF(plotId);
	//-----------------------------------------------------
	dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->compositeMainDispBuffer();
}

void PointDataHandler::setPCPLineDensityBuilderTF(int tfId)
{
	dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->setColorTF(g_cm_pcp_tfs[tfId], tfId);
	dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->setAlphaTF(g_alpha_pcp_tfs[tfId], tfId);
	//----------------------------------------------
	// Densitymap Rendering Methods:
	// global tf
	if (!g_params.IsDrawMaxLocs())
		dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->renderToDispBuffer(tfId); 
	// try local tf
	//dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->renderToDispBuffer_localTF(tfId);
	// try local tf with maximal locations
	else
		dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->renderToDispBuffer_max_filt_localTF(tfId);
	//-----------------------------------------------------
	dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->compositeMainDispBuffer();
}

QImage& PointDataHandler::createFocusContextDispImagePCPLineDensityBuilder(int firstFocusAxis, int focusDim)
{
	// DO NOT REDRAW display buffers???
	return dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->compositeMainDispBuffer(true, firstFocusAxis, focusDim - 1);
}

QImage& PointDataHandler::createDispImageFromAccumBuffer(int plotMode)
{
	// The drawing function for pcp as of June 1. 2016
	// Take log of the accumulation buffer and convert to the display buffer!
	cv::Mat& accumBuf = (plotMode == 0)? _pcpPlotBuilder->accumBuffer() : _splomPlotBuilder->accumBuffer();
	QImage&  dispBuf = (plotMode == 0)? _pcpPlotBuilder->dispBuffer() : _splomPlotBuilder->dispBuffer();

	// Clean dispBuf
	QPainter painter;
	painter.begin(&dispBuf);
	painter.fillRect(dispBuf.rect(), QColor(0, 0, 0, 0));
	painter.end();
	//////////////////

	//accumBuf = accumBuf + 1;
	double oneOverMaxLog = (g_params.Pcp_RenderMode() == PR_DS_LINE_DENSITY) ? 1.0 / g_pcp_downsample_max_log_val : 1.0 / g_params.Max_log_val();
	if (plotMode == 0 && g_params.doClustering() && g_params.numClusters() > 1)
	{
			vector<cv::Mat>& clusterBuffers = dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->clusterIdBuffer();
			bool chooselowClusterId = true;
			for (int r = 0; r < accumBuf.rows; r++)
			{
				for (int c = 0; c < accumBuf.cols; c++)
				{
					double val = log(accumBuf.at<double>(r, c) + 1);

					//double val = (plotMode == 0) ? log(accumBuf.at<double>(r, c) + 1) : log(accumBuf.at<double>(r, c) + 1) / log(2.0);
					double normVal = val * oneOverMaxLog;
					//double mappedVal = (plotMode == 1)? 1.0 - pow(normVal, 1.0 / (1.0 + g_splom_render_scale)) : 1.0 - normVal;
					double Y = pow(normVal, 1.0 / (1.0 + g_splom_render_scale));
					if (g_visColorScheme == CS_LIGHT)
						Y = 1.0 - Y;
					Y = MAX(MIN(1.0, Y), 0.0);


					// check cluster color
					int selClusterId = -1;
					for (int clusterId = 0; clusterId < g_params.numClusters(); clusterId++)
					{

						if (clusterBuffers[clusterId].at<unsigned char>(r, c) > 0)
						{
							selClusterId = clusterId;

							if (chooselowClusterId)
								break;
						}
					}

					if (selClusterId >= 0)
					{
						QColor rgbCol = QColor(g_clusterColors[selClusterId].x, g_clusterColors[selClusterId].y, g_clusterColors[selClusterId].z);
						QColor hsvCol = rgbCol.toHsv();
						int hsvVal = MIN(255, int(Y * 255.0f));

						hsvCol.setHsv(hsvCol.hsvHue(), hsvCol.hsvSaturation(), hsvVal);
						rgbCol = hsvCol.toRgb();
						dispBuf.setPixel(c, r, qRgb(rgbCol.red(), rgbCol.green(), rgbCol.blue()));
					}
					else
						dispBuf.setPixel(c, r, qRgb(int(255.0 * Y), int(255.0 * Y), int(255.0 * Y)));
				}
			}
		return dispBuf;
	}


	for (int r = 0; r < accumBuf.rows; r++)
	{
		for (int c = 0; c < accumBuf.cols; c++)
		{
			double val = log(accumBuf.at<double>(r, c) + 1);

			//double val = (plotMode == 0) ? log(accumBuf.at<double>(r, c) + 1) : log(accumBuf.at<double>(r, c) + 1) / log(2.0);
			double normVal = val * oneOverMaxLog;
			//double mappedVal = (plotMode == 1)? 1.0 - pow(normVal, 1.0 / (1.0 + g_splom_render_scale)) : 1.0 - normVal;
			double mappedVal = pow(normVal, 1.0 / (1.0 + g_splom_render_scale));
			if (g_visColorScheme == CS_LIGHT)
				mappedVal = 1.0 - mappedVal;
			mappedVal = MAX(MIN(1.0, mappedVal), 0.0);

			dispBuf.setPixel(c, r, qRgb(int(255.0 * mappedVal), int(255.0 * mappedVal), int(255.0 * mappedVal)));
		}
	}

	// Draw coordinates over!
	//dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->drawCoordinates();
	return dispBuf;
}

QImage& PointDataHandler::createFocusContextDispImageFromAccumBuffer(int firstFocusAxis, int focusDim, int plotMode)
{

	// The drawing function for pcp as of June 1. 2016
	// Take log of the accumulation buffer and convert to the display buffer!
	cv::Mat& accumBuf = (plotMode == 0) ? _pcpPlotBuilder->accumBuffer() : _splomPlotBuilder->accumBuffer();
	QImage&  dispBuf = (plotMode == 0) ? _pcpPlotBuilder->dispBuffer() : _splomPlotBuilder->dispBuffer();
	//accumBuf = accumBuf + 1;
	double oneOverMaxLog = (g_params.Pcp_RenderMode() == PR_DS_LINE_DENSITY) ? 1.0 / g_pcp_downsample_max_log_val : 1.0 / g_params.Max_log_val();
	for (int r = 0; r < accumBuf.rows; r++)
	{
		for (int c = 0; c < accumBuf.cols; c++)
		{
			double val = log(accumBuf.at<double>(r, c) + 1);
			FLOATVECTOR2 pcpPt;
			DispBuffer2InselbergCoords(INTVECTOR2(c, r), pcpPt);
			double normVal = val * oneOverMaxLog;
			double mappedVal = pow(normVal, 1.0 / (1.0 + g_splom_render_scale));
			if (g_visColorScheme == CS_LIGHT)
				mappedVal = 1.0 - mappedVal;
			int alpha = 255;
			if (pcpPt.x < float(firstFocusAxis) || pcpPt.x > float(firstFocusAxis + focusDim)) // use very low alpha color for context axes
			{
				//alpha = 50; // the default configuration
				alpha = 0; // only for showing ONE subdimension!
			}
			if (g_visColorScheme == CS_LIGHT)
				dispBuf.setPixel(c, r, qRgba(int(255.0 * mappedVal), int(255.0 * mappedVal), int(255.0 * mappedVal), alpha));
			else
			{
				int pixelVal = int(double(alpha) * mappedVal);
				dispBuf.setPixel(c, r, qRgba(pixelVal, pixelVal, pixelVal, 255));
			}
		}
	}
	return dispBuf;
}

QImage& PointDataHandler::createFocusContextPCPDispImageRedraw(int firstFocusAxis, int focusDim)
{
	//cv::Mat& accumBuf = _pcpPlotBuilder->accumBuffer();
	QImage&  dispBuf = _pcpPlotBuilder->dispBuffer();
/*
	for (int r = 0; r < accumBuf.rows; r++)
	{
		for (int c = 0; c < accumBuf.cols; c++)
		{
			double val = log(accumBuf.at<double>(r, c) + 1);
			FLOATVECTOR2 pcpPt;
			DispBuffer2InselbergCoords(INTVECTOR2(c, r), pcpPt);
			double normVal = val / g_params.Max_log_val();
			double mappedVal = 1.0 - pow(normVal, 1.0 / (1.0 + g_splom_render_scale));
			int alpha = 50;
			dispBuf.setPixel(c, r, qRgba(int(255.0 * mappedVal), int(255.0 * mappedVal), int(255.0 * mappedVal), alpha));
		}
	}
	*/// Clear focus buffer
	dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->clearAccumBuffer(-1);
	// Redraw ALL samples in subspace
	for (vector<vector<float>>::iterator IT = g_sampled_raw.begin(); IT != g_sampled_raw.end(); ++IT)
		dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->processSampleInFocus(*IT, firstFocusAxis, focusDim);

	cv::Mat& focusAccumBuf = dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->focusPCPAccumBuffer();
	QImage& focusDispBuf = dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->focusPCPbuffer();

	for (int r = 0; r < focusAccumBuf.rows; r++)
	{
		for (int c = 0; c < focusAccumBuf.cols; c++)
		{
			double val = log(focusAccumBuf.at<double>(r, c) + 1);
			FLOATVECTOR2 pcpPt;
			DispBuffer2InselbergCoords(INTVECTOR2(c, r), pcpPt);
			double normVal = val / g_params.Max_log_val();
			double mappedVal = 1.0 - pow(normVal, 1.0 / (1.0 + g_splom_render_scale));
			int alpha = 255;
			if (mappedVal == 0.0)
				alpha = 0;
			
			focusDispBuf.setPixel(c, r, qRgba(int(255.0 * mappedVal), int(255.0 * mappedVal), int(255.0 * mappedVal), alpha));
		}
	}
	// Blend the two display buffers
	// Blend Over dispBuf
	QPixmap tmp(dispBuf.width(), dispBuf.height());
	tmp.fill(Qt::white);
	QPainter tmpPainter(&tmp);

	tmpPainter.setCompositionMode(QPainter::CompositionMode_SourceOver);
	tmpPainter.setRenderHint(QPainter::Antialiasing);

	tmpPainter.setOpacity(1.0);
	tmpPainter.drawImage(0, 0, dispBuf);
	tmpPainter.setOpacity(0.9);
	tmpPainter.drawImage(0, 0, focusDispBuf);
	tmpPainter.end();
	dispBuf = tmp.toImage().copy();
	return dispBuf;
}

// NOTE:
// The actual display buffer drawing function where color coding is performed!!!
QImage& PointDataHandler::createDispImageFrom1FlatAccumBuffer(int plotMode)
{
	// The drawing function for pcp as of June 1. 2016
	// Take log of the accumulation buffer and convert to the display buffer!
	cv::Mat& accumBuf = (plotMode == 0) ? dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->accumBuffer1Flat() :
		dynamic_cast<SplomPlotBuilder*>(_splomPlotBuilder)->accumBuffer1DTrend();
	QImage&  dispBuf = (plotMode == 0) ? dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->oneFlatBuffer() :
		dynamic_cast<SplomPlotBuilder*>(_splomPlotBuilder)->oneDTrendBuffer();

	QColor rgbCol = Qt::darkGreen; // QColor(g_sc_1flat_layer_color.x, g_sc_1flat_layer_color.y, g_sc_1flat_layer_color.z);
	QColor hsvCol = rgbCol.toHsv();
	//QPainter painter(&dispBuf);

	if (accumBuf.channels() == 1)
	{
		cv::Mat logAccBuf(accumBuf);
		logAccBuf += 1;
		cv::log(accumBuf, logAccBuf);
		for (int r = 0; r < accumBuf.rows; r++)
		{
			for (int c = 0; c < accumBuf.cols; c++)
			{
				double val = logAccBuf.at<double>(r, c);//log(accumBuf.at<double>(r, c) + 1);
				double normVal = val / g_params.Max_log_val();
				float mappedVal = float(pow(normVal, 1.0 / (1.0 + g_splom_render_scale)));

				int hsvVal = (val == 0.0) ? 0 : MIN(255, int(mappedVal * 255.0f));

				hsvCol.setHsv(hsvCol.hsvHue(), hsvCol.hsvSaturation(), hsvVal);
				int alpha = (val == 0.0) ? 0 :  MIN(255, int(255.0f * val));
				hsvCol.setAlpha(alpha);
				rgbCol = hsvCol.toRgb();
				rgbCol.setAlpha(alpha);
				dispBuf.setPixel(c, r, hsvCol.toRgb().rgba());
			}
		}
	}
	else // We have 2 or more channels, do colormapping based on sid
	{
		FLOATVECTOR3 overlapCol = g_params.colorBothIdxPts();
		FLOATVECTOR3 vdOnlyCol = g_params.colorVdIdxPts();
		for (int r = 0; r < accumBuf.rows; r++)
		{
			for (int c = 0; c < accumBuf.cols; c++)
			{
				double val = log(accumBuf.at<cv::Vec2d>(r, c)[0] + 1) ;

				int sid = int(accumBuf.at<cv::Vec2d>(r, c)[1]);
				
				if (sid == 0 && g_params.Pcp_compare_mode() >= 1)
					continue; // skip if the id value is 0 
				int id = sid;
				//float mappedVal = val;
				double normVal = /*2.0 * */val / g_params.Max_log_val();
				double mappedVal = pow(normVal, 1.0 / (1.0 + g_pcp_render_scale)); // Gamma correction
				float factor = 0.5f;
				// NOTE: use ONLY the difference/overlapping colors
				float falpha = 1.0f;
			
				if (sid >= 0 && sid < BUF_OVERLAP_ID) // has only the original sid for the indexed points
				{
					id = sid - 1;
					if (id < 0)
						continue;
					if (g_params.Pcp_compare_mode() == 0)
						rgbCol = QColor(g_pcp_1flat_colormap[id].x, g_pcp_1flat_colormap[id].y, g_pcp_1flat_colormap[id].z);
				}				
				else if (sid >= BUF_OVERLAP_ID) // has overlapping indexed points
				{ 
					id = sid - BUF_OVERLAP_ID;
					if (g_params.Pcp_compare_mode() == 1) // diff compare mode
					{
						// blend with blue
						rgbCol = QColor(g_pcp_1flat_colormap[id].x * factor * (1.0f - falpha) + falpha * overlapCol.x,
							g_pcp_1flat_colormap[id].y * factor * (1.0f - falpha) + falpha * overlapCol.y,
							g_pcp_1flat_colormap[id].z * factor * (1.0f - falpha) + falpha * overlapCol.z);
					}
					else if (g_params.Pcp_compare_mode() == 2) // Set 2 mode
					{
						rgbCol = QColor(g_pcp_1flat_colormap[id].x, g_pcp_1flat_colormap[id].y, g_pcp_1flat_colormap[id].z);
					}
					else
					{
						rgbCol = Qt::red;
						//rgbCol = QColor(g_pcp_1flat_colormap[id].x, g_pcp_1flat_colormap[id].y, g_pcp_1flat_colormap[id].z); // Set 1 mode
					}
						
				}else
				{
					id = -(sid + 1);
					if (g_params.Pcp_compare_mode() == 1) // diff compare mode
					{
						// Set 2 only, blend with red
						rgbCol = QColor(g_pcp_1flat_colormap[id].x * factor * (1.0f - falpha) + falpha * vdOnlyCol.x,
							g_pcp_1flat_colormap[id].y * factor * (1.0f - falpha) + falpha * vdOnlyCol.y,
							g_pcp_1flat_colormap[id].z * factor * (1.0f - falpha) + falpha * vdOnlyCol.z);
					}
					else if (g_params.Pcp_compare_mode() == 2)// Set 2 mode
					{
						rgbCol = QColor(g_pcp_1flat_colormap[id].x, g_pcp_1flat_colormap[id].y, g_pcp_1flat_colormap[id].z);
					}
					else
						break;
				}				
				hsvCol = rgbCol.toHsv();
				int h, s, v;
				hsvCol.getHsv(&h, &s, &v);
				int alpha = 0;
				if (val != 0.0 && (g_params.FirstAxis_subspace() == -1 || g_params.FirstAxis_subspace() == id))
				{
					v = MAX(0, MIN(255, int(2.0 * mappedVal * (v + 50))));
					alpha = MAX(0, MIN(255, int(255.0 * mappedVal/*2.0 * val / g_params.Max_log_val()*/)));
				}
		
				hsvCol.setHsv(h, s, v, alpha);
				rgbCol = hsvCol.toRgb();
				dispBuf.setPixel(c, r, rgbCol.rgba());
			}
		}
	}
	//painter.end();
	return dispBuf;
}

QImage& PointDataHandler::createDispImagePosCorr1FlatAccumBuffer()
{
	// Abuse the focus pcp buffer for positive correlated samples!!!
	cv::Mat& accumBuf = dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->focusPCPAccumBuffer();
	QImage&  dispBuf = dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->focusPCPbuffer();

	QColor rgbCol(255, 127, 0);
	QColor hsvCol = rgbCol.toHsv();
	//QPainter painter(&dispBuf);
	double logmax = g_params.Max_log_val() - log(2); // as we have only half the values
	if (accumBuf.channels() == 1)
	{
		cv::Mat logAccBuf(accumBuf);
		logAccBuf += 1;
		cv::log(accumBuf, logAccBuf);
		for (int r = 0; r < accumBuf.rows; r++)
		{
			for (int c = 0; c < accumBuf.cols; c++)
			{
				double val = logAccBuf.at<double>(r, c);//log(accumBuf.at<double>(r, c) + 1);
				double normVal = val / logmax;
				float mappedVal = float(pow(normVal, 1.0 / (1.0 + g_splom_render_scale)));

				int hsvVal = (val == 0.0) ? 0 : MIN(255, int(mappedVal * 255.0f));

				hsvCol.setHsv(hsvCol.hsvHue(), hsvCol.hsvSaturation(), hsvVal);
				int alpha = (val == 0.0) ? 0 : MIN(255, int(255.0f * val));
				hsvCol.setAlpha(alpha);
				rgbCol = hsvCol.toRgb();
				rgbCol.setAlpha(alpha);
				dispBuf.setPixel(c, r, hsvCol.toRgb().rgba());
			}
		}
	}
	else // We have 2 or more channels, do colormapping based on sid
	{
		for (int r = 0; r < accumBuf.rows; r++)
		{
			for (int c = 0; c < accumBuf.cols; c++)
			{
				double val = log(accumBuf.at<cv::Vec2d>(r, c)[0] + 1);

				int sid = int(accumBuf.at<cv::Vec2d>(r, c)[1]);
				//rgbCol = QColor(g_pcp_1flat_colormap[sid].x, g_pcp_1flat_colormap[sid].y, g_pcp_1flat_colormap[sid].z);
				//float mappedVal = val;
				double normVal = /*2.0 * */val / logmax;
				double mappedVal = pow(normVal, 1.0 / (1.0 + g_pcp_render_scale)); // Gamma correction
				hsvCol = rgbCol.toHsv();
				int h, s, v;
				hsvCol.getHsv(&h, &s, &v);
				int alpha = 0;
				if (val != 0.0 && (g_params.FirstAxis_subspace() == -1 || g_params.FirstAxis_subspace() == sid))
				{
					v = MAX(0, MIN(255, int(2.0 * mappedVal * (v + 50))));
					alpha = MAX(0, MIN(255, int(255.0 * mappedVal/*2.0 * val / g_params.Max_log_val()*/)));
				}

				hsvCol.setHsv(h, s, v, alpha);
				rgbCol = hsvCol.toRgb();
				dispBuf.setPixel(c, r, rgbCol.rgba());
			}
		}
	}	
	dispBuf.save("PosCorrBuf.png", "PNG");
	return dispBuf;
}

QImage& PointDataHandler::createDispImageNegCorr1FlatAccumBuffer()
{
	// Abuse the focus pcp buffer for positive correlated samples!!!
	cv::Mat& accumBuf = dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->accumBuffer1Flat();
	QImage&  dispBuf = dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->oneFlatBuffer();

	QColor rgbCol(0,100,255);
	//QColor rgbCol(55, 184, 126);
	//QColor rgbCol = QColor(g_pcp_1flat_colormap[0].x, g_pcp_1flat_colormap[0].y, g_pcp_1flat_colormap[0].z);
	QColor hsvCol = rgbCol.toHsv();
	//QPainter painter(&dispBuf);
	double logmax = g_params.Max_log_val() - log(2); // as we have only have the values

	if (accumBuf.channels() == 1)
	{
		cv::Mat logAccBuf(accumBuf);
		logAccBuf += 1;
		cv::log(accumBuf, logAccBuf);
		for (int r = 0; r < accumBuf.rows; r++)
		{
			for (int c = 0; c < accumBuf.cols; c++)
			{
				double val = logAccBuf.at<double>(r, c);//log(accumBuf.at<double>(r, c) + 1);
				double normVal = val / g_params.Max_log_val();
				float mappedVal = float(pow(normVal, 1.0 / (1.0 + g_splom_render_scale)));

				int hsvVal = (val == 0.0) ? 0 : MIN(255, int(mappedVal * 255.0f));

				hsvCol.setHsv(hsvCol.hsvHue(), hsvCol.hsvSaturation(), hsvVal);
				int alpha = (val == 0.0) ? 0 : MIN(255, int(255.0f * val));
				hsvCol.setAlpha(alpha);
				rgbCol = hsvCol.toRgb();
				rgbCol.setAlpha(alpha);
				dispBuf.setPixel(c, r, hsvCol.toRgb().rgba());
			}
		}
		cout << "Here, path 1" << endl;
	}
	else // We have 2 or more channels, do colormapping based on sid
	{
		for (int r = 0; r < accumBuf.rows; r++)
		{
			for (int c = 0; c < accumBuf.cols; c++)
			{
				double val = log(accumBuf.at<cv::Vec2d>(r, c)[0] + 1);

				int sid = int(accumBuf.at<cv::Vec2d>(r, c)[1]);
				//rgbCol = QColor(g_pcp_1flat_colormap[sid].x, g_pcp_1flat_colormap[sid].y, g_pcp_1flat_colormap[sid].z);
				//float mappedVal = val;
				double normVal = /*2.0 * */val / logmax;
				double mappedVal = pow(normVal, 1.0 / (1.0 + g_pcp_render_scale)); // Gamma correction
				hsvCol = rgbCol.toHsv();
				int h, s, v;
				hsvCol.getHsv(&h, &s, &v);
				int alpha = 0;
				if (val != 0.0 && (g_params.FirstAxis_subspace() == -1 || g_params.FirstAxis_subspace() == sid))
				{
					v = MAX(0, MIN(255, int(2.0 * mappedVal * (v + 50))));
					alpha = MAX(0, MIN(255, int(255.0 * mappedVal/*2.0 * val / g_params.Max_log_val()*/)));
				}

				hsvCol.setHsv(h, s, v, alpha);
				rgbCol = hsvCol.toRgb();
				dispBuf.setPixel(c, r, rgbCol.rgba());
			}
		}
		cout << "Here, path 2" << endl;
	}
	dispBuf.save("NegCorrBuf.png", "PNG");
	return dispBuf;
}

QImage& PointDataHandler::createDispImageFrom2FlatAccumBuffer()
{
	// The drawing function for pcp as of June 1. 2016
	// Take log of the accumulation buffer and convert to the display buffer!
	cv::Mat& accumBuf = dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->accumBuffer2Flat();
	QImage&  dispBuf = dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->twoFlatBuffer();
	QColor rgbCol = QColor(g_pcp_1flat_layer_color.x, g_pcp_1flat_layer_color.y, g_pcp_1flat_layer_color.z);
	QColor hsvCol = rgbCol.toHsv();

	if (accumBuf.channels() == 1)
	{
		for (int r = 0; r < accumBuf.rows; r++)
		{
			for (int c = 0; c < accumBuf.cols; c++)
			{
				double val = log(accumBuf.at<double>(r, c) + 1);

				float mappedVal = val /*/ g_pcp_render_scale*/;

				int hsvVal = (val == 0.0) ? 0 : MIN(255, int(mappedVal * 255.0f));

				hsvCol.setHsv(hsvCol.hsvHue(), hsvCol.hsvSaturation(), hsvVal);
				int alpha = (val == 0.0) ? 0 : 255;// MIN(255, int(255.0f * val));
				hsvCol.setAlpha(alpha);
				rgbCol = hsvCol.toRgb();
				dispBuf.setPixel(c, r, hsvCol.toRgb().rgba());
			}
		}
	}
	else // We have 2 or more channels, do colormapping based on sid
	{
		for (int r = 0; r < accumBuf.rows; r++)
		{
			for (int c = 0; c < accumBuf.cols; c++)
			{
				double val = log(accumBuf.at<cv::Vec2d>(r, c)[0] + 1);
				int sid = int(accumBuf.at<cv::Vec2d>(r, c)[1]);
				double normVal =/*2.0 **/ val / g_params.Max_log_val();
				double mappedVal = pow(normVal, 1.0 / (1.0 + g_pcp_render_scale)); // Gamma correction
				// Now: draw point using 1-flat colormap
				rgbCol = QColor(g_pcp_1flat_colormap[sid].x, g_pcp_1flat_colormap[sid].y, g_pcp_1flat_colormap[sid].z);
				hsvCol = rgbCol.toHsv();
				int h, s, v;
				hsvCol.getHsv(&h, &s, &v);
				int alpha = 0;
				if (val != 0.0 && (g_params.FirstAxis_subspace() == -1 || g_params.FirstAxis_subspace() == sid))
				{
					v = MAX(0, MIN(255, int(2.0 * mappedVal * (v + 50))));
					alpha = MAX(0, MIN(255, int(255.0 * mappedVal)));
				}
				//rgbCol.setAlpha(alpha);
				hsvCol.setHsv(h, s, v, alpha);
				//hsvCol.setHsv(h, MIN(255, s + 50), v, alpha);
				rgbCol = hsvCol.toRgb();
				dispBuf.setPixel(c, r, rgbCol.rgba());
			}
		}
	}
	return dispBuf;
}

void PointDataHandler::Norm2InselbergCoords(float x, float y, float& nx, float& ny)
{
	dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->Norm2InselbergCoords(nx, x);
	y = ny;
}

void PointDataHandler::InselbergCoords2InselNorm(float nx, float ny, float& x, float& y)
{
	dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->InselbergCoords2InselNorm(x, nx);
}

void PointDataHandler::DispBuffer2InselbergCoords(INTVECTOR2 scrPt, FLOATVECTOR2& pcpPt)
{
	dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->DispBuffer2InselbergCoords(scrPt, pcpPt);
}

void PointDataHandler::InselbergCoords2DispBuffer(FLOATVECTOR2 pcpPt, INTVECTOR2& scrPt)
{
	dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->InselbergCoords2DispBuffer(pcpPt, scrPt);
}

cv::Mat PointDataHandler::findDensityPlotROI(int dimId)
{
	cv::Mat& buf = dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->densityPlotAccumBuffer(dimId);
	cv::Mat& mask = dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->densityPlotTFMask(dimId);
	// Try different ROI finding techniques
	ROI_METHOD roi_tech =  ROI_LOCAL_PERCENTILE_CORNER;
	//ROI_METHOD roi_tech = ROI_CORNER;
	cv::Mat maxLoc;
	switch (roi_tech)
	{
	case ROI_CORNER:
		maxLoc = PCimageProcessor::findCorners(buf, mask);
		break;
	case ROI_LOCAL_PERCENTILE_CORNER:
		maxLoc = PCimageProcessor::findCornersLocalPercentile(buf, mask);
		break;
	case ROI_LOCAL_MAX_CORNER:
		maxLoc = PCimageProcessor::getMaximalLocs(buf, mask);
		break;
	default:
		maxLoc = PCimageProcessor::findCornersLocalPercentile(buf, mask);
		break;
	}
	//cv::Mat maxLoc = PCimageProcessor::findCornersLocalPercentile(buf, mask);
	// find local maximum and then evaluate eigen values
	return maxLoc;
}

void PointDataHandler::findDensityPlotROI()
{
	QMessageBox pleaseWait;
	pleaseWait.setText(QString("Please wait...computing ROI"));
	pleaseWait.show();
	vector<cv::Mat> roi = g_params.DensMapROILocs();
	if (g_params.FirstAxis_subspace() != -1)
	{
		pleaseWait.setText(QString("Please wait...computing ROI for subdimension %1").arg(g_params.FirstAxis_subspace()));
		
		roi[g_params.FirstAxis_subspace()] = findDensityPlotROI(g_params.FirstAxis_subspace());

	}
	else
	{
		for (int i = 0; i < int(roi.size()); i++)
		{
			pleaseWait.setText(QString("Please wait...computing ROI for subdimension %1").arg(i));
			roi[i] = findDensityPlotROI(i);

		}
	}
	g_params.DensMapROILocs(roi);
	pleaseWait.close();
}

void PointDataHandler::updateDensityMapRenderMasks()
{
	vector<QImage> renderMask = g_params.DensMapRenderMask();
	// draw to rendermask
	cv::Mat cornerMask;
	float radius = 15.0f;
	vector<cv::Mat> roi = g_params.DensMapROILocs();
	if (g_params.FirstAxis_subspace() == -1)
	{
		// change all bufs!
	
		for (size_t i = 0; i < roi.size(); i++)
		{
			if (g_params.DensMapCornerThres() == 0){ // set all pass
				renderMask[i].fill(Qt::white);
				continue;
			}

			renderMask[i].fill(0);
			QPainter painter(&renderMask[i]); 
			painter.setPen(Qt::NoPen);
			painter.setCompositionMode(QPainter::CompositionMode_SourceOver);

			roi[i].convertTo(cornerMask, CV_32FC1);
			cv::resize(cornerMask, cornerMask, cv::Size(g_params.Pcp_densPlotDispBufSize().x, g_params.Pcp_densPlotDispBufSize().y));
			for (int r = 0; r < cornerMask.rows; r++)
			{
				for (int c = 0; c < cornerMask.cols; c++)
				{
					if (cornerMask.at<float>(r, c) >= 255.0f * g_params.DensMapCornerThres())
					{
						QRadialGradient rg(QPointF(c, r), radius);
						rg.setColorAt(0, QColor(255, 255, 255, 255));
						rg.setColorAt(0.65, QColor(0, 0, 0, 0));
						painter.setBrush(QBrush(rg));
						painter.drawEllipse(QPoint(c, r), int(radius), int(radius));
					}
				}
			}
			painter.end();
		
		}
	}
	else
	{
		// change the selected buff
		int id = g_params.FirstAxis_subspace();
		if (g_params.DensMapCornerThres() == 0){ // set all pass
			renderMask[id].fill(Qt::white);
		}
		else
		{
			renderMask[id].fill(0);
			cv::resize(roi[id], cornerMask, cv::Size(g_params.Pcp_densPlotDispBufSize().x, g_params.Pcp_densPlotDispBufSize().y));
			QPainter painter(&renderMask[id]);
			painter.setPen(Qt::NoPen);
			painter.setCompositionMode(QPainter::CompositionMode_SourceOver);

			for (int r = 0; r < cornerMask.rows; r++)
			{
				for (int c = 0; c < cornerMask.cols; c++)
				{
					if (cornerMask.at<float>(r, c) >= 255.0f * g_params.DensMapCornerThres())
					{
						QRadialGradient rg(QPointF(c, r), radius);
						rg.setColorAt(0, QColor(255, 255, 255, 255));
						rg.setColorAt(0.65, QColor(0, 0, 0, 0));
						painter.setBrush(QBrush(rg));
						painter.drawEllipse(QPoint(c, r), int(radius), int(radius));
					}
				}
			}
			painter.end();
			// write out painter image
			renderMask[id].save("cornerThresMap.png", "PNG");
		}
	}


	g_params.DensMapRenderMask(renderMask);
}