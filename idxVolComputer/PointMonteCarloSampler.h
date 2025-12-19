#pragma once
//#include "DirectionalMonteCarloSampler.h"
#include "sampling.h"
#include <QImage>
#include "XPCPSample.h"
#include "global.h"

class PointMonteCarloSampler : public DirectionalMonteCarloSampler
{
public:
	PointMonteCarloSampler(const std::vector<VolumeData*>& volData, UINT64VECTOR3 volDim, std::string& strOutputFile);
	virtual ~PointMonteCarloSampler();

	bool init(float portion_voxels, int neighborSize);
	bool init(INTVECTOR3 sample_dim, int neighborSize);
	void sample(float portion_voxels, int neighborSize);
	virtual bool updatePlot(); 
	void finishPlot(bool buildAllStreamlines); // Process after drawing the plot
	
	virtual bool updatePCPlot();
	virtual bool updateSCnPCPlots_from_streamline();

	// Sample voxel with interpolation [NOTE: without normalization!!!]
	static std::vector<float> sampleVoxelNormCoords(FLOATVECTOR3 c_nl, const std::vector<VolumeData*>& volList);
	// With normalization
	std::vector<float> sampleVoxelWithInterpNormCoords(FLOATVECTOR3 c_nl); // use normalized coordinates
	// Without interpolation
	static std::vector<float> readVoxelNormalCoords(FLOATVECTOR3 c_nl, const std::vector<VolumeData*>& volList);

	void testPCP2Flat(); // Test the construction of 2-flat in the PCP

	//-------------------------------------------------------
	// Switch to the value domain!!!
	// Build the value domain by Monte-Carlo sampling the spatial domain!
	virtual std::vector<std::vector<float>> buildValueDomain(float portion_voxels, std::vector<std::vector<float>>& sampledData);
	virtual std::vector<std::vector<float>> sampleSpatialDomain(float portion_voxels, std::vector<std::vector<float>>& sampledData);
	// build the value domain by down-sampling the spatial domain (skipping certain voxels per axis)! Will cause artifact intentionally!
	std::vector<std::vector<float>> buildValueDomain_regDownsample_volDimSamplePos(INTVECTOR3 sampleDim, std::vector<std::vector<float>>& sampledData);

	std::vector<std::vector<float>> buildValueDomain_regDownsample_volDimSamplePos(FLOATVECTOR3 samplePortion, std::vector<std::vector<float>>& sampledData);
	// Vector-based volume data 
	std::vector<std::vector<float>> buildVecValueDomain_regDownsample_volDimSamplePos(FLOATVECTOR3 sampleDim, std::vector<std::vector<float>>& sampledData);

	vector<vector<float>> PointMonteCarloSampler::buildVecValueDomain_regDownsample_volDimSamplePos(INTVECTOR3 sampleDim, vector<vector<float>>& sampledData);
	void setXpcpData(std::vector<XPCPSample>& xpcpData)
	{
		_xpcpData = xpcpData;
	}
	//-------------------------------------------------------
private:
	void sampleCubeNeighborhood(UINT64VECTOR3 c, int neighborSize); // Sample discrete location with neighborhood
	void sampleCubeNeighborhoodNormCoords(FLOATVECTOR3 c_nl, int neighborSize1D); // Sample continuous normalized location with neighborhood
	// ----------------------------
	// Functions for update plot
	bool updatePlotFromLayerData();			// update plot from already computed layer data
	bool updatePlotFromVolumeSampling();    // update plot from directly sampling the volume
	//------------------------------
	// Sampling without neighborhood
	std::vector<float> sampleVoxelWithInterp(FLOATVECTOR3 p);// override functions

	///////////////////////////////////////////////////////////////
	// Process multivariate values within the neighborhood
	// TODO: We should really put the processing codes in some statistical analyzer class 
	void calcCovMatrix(const std::vector<float>& centralVal, const std::vector<std::vector<float>>& valList); // Calculate a covariance matrix and plot it at centralVal
	bool calcPlane(const std::vector<float>& v1, const std::vector<float>& v2, const std::vector<float>& p,  PlaneCoeffs& Pl); // Calculate the plane coefficient from vectors v1, and v2. Return if the vectors are valid?
	void calcLinearRegScatterplot(FLOATVECTOR2 centralVal);
	//////////////////////////////////////////////////////////////

	void accumSampleDensity(const FLOATVECTOR2& val);
	void processSampleScatterplot(const std::vector<float>& vals);
private:
	float  _sample_portion;
	int    _neighborSize; 
};