#pragma once
#include "VolumeData.h"
#include "PointDataHandler.h"

class MTRand;
class Film;
class PlotBuilder;
struct DirectionalSample
{
	UINT64VECTOR3 _pos; // sample location
	FLOATVECTOR3  _dir; // sample direction as a vector
	float         _len; // length of the sample direction

	friend std::ostream& operator<<(std::ostream& os, const DirectionalSample& obj);
	friend std::istream& operator>>(std::istream& is, DirectionalSample& obj);
};

class DirectionalMonteCarloSampler : public PointDataHandler
{
public:
	DirectionalMonteCarloSampler(const std::vector<VolumeData*>& volData, UINT64VECTOR3 volDim, const std::string& strOutputFile);
	~DirectionalMonteCarloSampler(void);

	void sample(float portion_voxels, int lines_per_voxel);
	void rayTracingSampling(int num_dir, int sample_img_width, int sample_img_height);

	//-------------------------------------------------------
	// Switch to the value domain!!!
	// Build the value domain by sampling the spatial domain!
	// Also returns sample locations
	virtual std::vector<std::vector<float>> buildValueDomain(float portion_voxels, std::vector<std::vector<float>>& sampledData) = 0;

	// Triliniear interpolation of volume's value at any given location [NOTE: No normalization!]
	static std::vector<float> trilinear(FLOATVECTOR3 p, const std::vector<VolumeData*>& volList);

protected:
	// NOTE: the reason for having the sample generation function separately is that the samples may be reused afterward, and we don't have to generate them again.

	void generatesamples(std::vector<DirectionalSample>& samples, float portion_voxels, int lines_per_voxel);
	void sampleVolume(const std::vector<DirectionalSample>& samples);
	void sampleLine(UINT64VECTOR3 start, FLOATVECTOR3 dir, float len);
	void sampleLineWithInterp(FLOATVECTOR3 start, FLOATVECTOR3 dir, float len); // sample line with interplolation
	void sampleVoxel(UINT64VECTOR3 p);
	void sampleVoxelWithInterp(FLOATVECTOR3 p); // with interpolation
	bool processCheckpoint(); // Stop here and do some check-ups

	// Directional sampling a unit sphere
	void generateSphereDirectionalSamples(std::vector<FLOATVECTOR3>& samples, UINT64 num_samples);
	void raySampleVolume(FLOATVECTOR3 eyePos); // Shoot parallel rays through the volume with direction viewDir
	//void sampleFullshpereDirection(UINT64 num_samples);
	FLOATVECTOR3 uniformSampleHemisphere(float u1, float u2);
	FLOATVECTOR3 uniformSampleFullsphere(float u1, float u2);

	float uniformHemispherePdf();
	float uniformFullspherePdf();

	// draw lines from one voxel's data value to another in the histogram 
	void lineToSampleVoxel(UINT64VECTOR3 p); // without interpolation
	void lineToSampleVoxelWithInterp(FLOATVECTOR3 p); // with interpolation
	// A 2D lineto function, essentially the same as sampleLine (in 3D).
	void lineTo2D(UINT64VECTOR2 start, UINT64VECTOR2 end);

	void incrScatterplotElemCnt(UINT64VECTOR2 p);
	void incrScatterplotElemCnt(FLOATVECTOR2 fp);
	// A skeleton function where the plot builder transform the sample into plot
	void processSample(const std::vector<float>& vals);
	void processSampleLineTo(const std::vector<float>& vals);
	void showPlot();
	void clearWorkingPlotBuf();
	
	void writeResult();
	// linear inteprolation
	void lerp(float x, float a);
	void bilinear();
	std::vector<float> trilinear(FLOATVECTOR3 p); // compute trilinear interpolation, with full range 3-vec voxel p; Assume each voxel has M-value, and output an array of float for the M-value vector


protected:
	std::vector<VolumeData*>       _volData;
	UINT64VECTOR3                  _volDim;
	double                         _volValNormalizer; // denom for normalizing the volume
	
	std::vector<DirectionalSample> _samples;
	UINT64*                        _scatterplot;

	UINT64VECTOR2                  _scatterplotDim;
	std::string                    _strOutputFile;

	UINT64VECTOR2                  _lastPosScatterplot; // records the last voxel's location in the scatterplot
	FLOATVECTOR2                   _fLastPosScatterplot; // normalized float version of the location
	std::vector<float>             _lastSampleVal; // last sample value (multivariate)
	MTRand*                        _rng;

	int                            _rayCheckpointSteps; // How many steps on the ray should we stop and check, e.g., determine if current sample has abnormal distribution?
	float                          _sampleStepSize; // the step size for ray sampling, in the unit of voxel (i.e., 0.5 is half voxel)
	FLOATMATRIX4                   _worldToTex; // world to tex transformation
	// Raytracing stuff
	Film*                          _film;
	bool                           _isRedraw; // indicates whether current drawing is redraw? If it's redraw, don't do statistics!!
};

