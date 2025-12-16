#pragma once
#include "XPCPSample.h"
#include "UncertainXPCPSample.h"
#include "Eigen/Core"
#include "global.h"
#include "MyStatistics.h"

// Handles forward and backward mapping of the high-dimensional data (infovis)
// Compute the 1-flat and 2-flat in the data domain (No more spatial domain!)
class ANNkd_tree;
class PCPInselberg;
class VolumeData;

enum PR_METHOD // Raw data processing method
{
	PR_SUBSPACE = 0,  // build from taking sub 3D spaces
	PR_FROM_LOWER_DIM // build from 1-flat
}; 

class XPCPdataAnalyzer
{
public:
	XPCPdataAnalyzer();
	~XPCPdataAnalyzer();
	
	////////////////////////////////////////////////
	// Static helper functions
	// Transformation of 1-flat indexed points!
	static void xform(double xin, double yin, int subspaceid, double& xout, double& yout); // the real coordinates input
	static void xform_idxpt_complete_line_desc(double c1, double c2, double c3, double& xout, double& yout); // tranformation using the homogeneous line coordinates with the origin at 0.5 (For mirror symmetry)
	static void xform_idxpt_complete_line_desc_noScaling(double c1, double c2, double c3, double& xout, double& yout); // tranformation using the homogeneous line coordinates without scaling with the origin at 0.5
	static void xform_idxpt_complete_line_desc_origin0(double c1, double c2, double c3, double& u, double& v); // tranformation using the homogeneous line coordinates with the origin at 0 (For indexed points)
	static void xform_idxpt_complete_line_desc_origin0_noScaling(double c1, double c2, double c3, double& xout, double& yout); // tranformation using the homogeneous line coordinates without scaling with the origin at 0 
	static void xform_nguyen_rosen_tvcg(double u, double v, double x0, double y0, double& q, double& r, bool& rotated);
	static void lineCoord2pointCoord2D(double c1, double c2, double c3, double& x, double& y); // convert line coordinate (3D) to 2D point coordinates
																							   // cluster HD data
	static void cluster(const std::vector<std::vector<float>>& rawData, std::vector<int>& labelData, int numK = 3);
	// Interpolation method(s) for eigenvectors
	// The interpolation scheme: 
	//Wang, Y., Jiang, F. & Liu, Y. Spectrum-sine interpolation framework for DTI processing. Med Biol Eng Comput 60, 279¨C295 (2022). https://doi.org/10.1007/s11517-021-02471-2

	static void eigenVecInterp1D(FLOATVECTOR3& vt, const FLOATVECTOR3& v1, const FLOATVECTOR3& v2, float t);
	static void eigenVecInterp1D(Eigen::VectorXf& vt, const Eigen::VectorXf& v1, const Eigen::VectorXf& v2, float t);
	static void eigenVecInterpND(FLOATVECTOR3& vt, const std::vector<FLOATVECTOR3>& v, const std::vector<float>& t);
	static void eigenVecInterpND(Eigen::VectorXf& vt, const std::vector<Eigen::VectorXf>& v, const std::vector<float>& t);
	// Interpolation of eigenvalues
	static void eigenValInterpND(float& lamt, const std::vector<float>& lambda, const std::vector<float>& t);

	// Matrix-based eigenvec,val interpolation
	// Mt, Lambt are outputs
	// M1, M2 are 3*3 matrices of eigenvectors (Three 3D vectors)
	// Lamb1, Lamb2 are 3D vectors (Three eigenvalues)
	static void matEigInterp1D(std::vector<Eigen::VectorXf>& Mt, Eigen::VectorXf& Lambt, const Eigen::MatrixXf& M1, const Eigen::MatrixXf& M2, const Eigen::VectorXf& Lamb1, const Eigen::VectorXf& Lamb2, float t);

	static void matEigInterpND(std::vector<Eigen::VectorXf>& Mt, Eigen::VectorXf& Lambt, const std::vector<Eigen::MatrixXf>& Ms, const std::vector<Eigen::VectorXf>& Lambs, const std::vector<float>& t);
	/////////////////////////////////////////////
	 


	// Raw Data Preprocessing
	void buildKDTree(const std::vector<std::vector<float>>& rawData, ANNkd_tree** pKdTree);
	////////////////////////////////////////////////////////
	// Process raw data in full dimension
	void processRawData(const std::vector<std::vector<float>>& rawData, int num_nearest_neighbors, PR_METHOD method = PR_SUBSPACE);

	// Process raw data in vector form: no need to compute local fittings anymore!
	void processVecRawData(const std::vector<std::vector<float>>& rawData, const std::vector<std::vector<float>>& samplePos);

	// Process raw data in trajectory form: no need to compute local fittings: use tangent for the 1-flat, and connect all 1-indexed points!
	void processTrajRawData(const std::vector<std::vector<float>>& rawData);

	// Process uncertain point data (distribution-based raw data)
	void processUncertainRawData(const std::vector<GaussianXd>& distrPointData, int num_nearest_neighbors);


	// Process raw data in the spatial domain with the neighborhood for continuous indexed points
	bool analyzeSpatialRawDataContIdxPts(const std::vector<VolumeData*>& volData, const std::vector<std::vector<float>>& samplePos, UINT64VECTOR3 idxVolDim, int numSamplesInNeighborhood = 50, FLOATVECTOR3 neighborR=FLOATVECTOR3(5.0f,5.0f,5.0f), bool writeOutFile = false); 
	// Volume loading framework for all types of volumes+

	bool processDataVol_load(const std::vector<VolumeData*>& volList, const std::vector<VolumeData*>& neighborInfoVolList, const std::vector<std::vector<float>>& samplePos);

	// Load precomptued indexed points data from volumes
	bool processDataVol_loadIdxVol(const std::vector<VolumeData*>& volList, const std::vector<VolumeData*>& idxVolData, const std::vector<VolumeData*>& idx2FlatsVolData, const std::vector<std::vector<float>>& samplePos);
	// Load precomptued eigenvector data from volumes
	bool processDataVol_loadEigVecVol(const std::vector<VolumeData*>& volList, const std::vector<VolumeData*>& majEigVols, const std::vector<VolumeData*>& secEigVols, const std::vector<std::vector<float>>& samplePos);
	////////////////////////////////////////////////////////
	// Process raw data in 2D or 3D subspaces
	void buildKDTreeInSubSpaces(const std::vector<std::vector<float>>& rawData, int subspaceId, int subspaceDim, ANNkd_tree** pKdTree);
	void processRawDataInSubSpaces(const std::vector<std::vector<float>>& rawData, int num_nearest_neighbors);
	// Accessor
	std::vector<XPCPSample>& xPCPData() { return g_xpcpData; }

	// make XPCPdata from the raw data. Compute 1-flat and 2-flat in the data domain.
	void mapToRawData(); // map brushed region back to the raw data

	std::vector<UINT64> detectOutliers(const std::vector<std::vector<float>>& rawData);

	// Eigenvectors interpolation test
	void testEigenVectorsInterp();


private:
	//*********************************************************
	// functions for computing p-flats

	// functions computing p-flats from eigenvectors
	void calc_p_flats_subSpaceMethod(const Eigen::VectorXf& mu, const Eigen::VectorXf& majEig, const Eigen::VectorXf& secEig, XPCPSample& xpcp_tuple); // construct 1-flat 2-flat from lower dimension
	// Use interactions of lower dimensional flats to compute p-flats
	void calc_p_flats_fromLowerDimMethod(const Eigen::VectorXf& mu, const Eigen::VectorXf& majEig, const Eigen::VectorXf& secEig, XPCPSample& xpcp_tuple); // construct 1-flat 2-flat using 3D subspaces.
	// Calculate 1-flats and 2-flats from 0-flats: p0: central point, p1: vPtMajEig, p2: vPtSecEig
	void calc_p_flats_fromLowFlats(const Eigen::VectorXf& mu, const vector<Point>& pcp_p0_no_repeat, const vector<Point>& pcp_p1_no_repeat, const vector<Point>& pcp_p2_no_repeat, XPCPSample& xpcp_tuple);

	// function computing p-flats from vector data: use sample location (pos) and vector value (vecVal). 
	void calc_p_flats_fromVecData(const vector<float>& vecVal, const vector<float>& pos, XPCPSample& xpcp_tuple);
	// compute p-flats using trajectory data: given p1 (current point), p2 (next point in time) 
	void calc_p_flats_fromTrajData(const vector<float>& p1, const vector<float>& p2, XPCPSample& xpcp_tuple);
	//**********************************************************

	void calc1Flat(const Eigen::VectorXf& mu, const Eigen::VectorXf& majEig, int ssId, FLOATVECTOR2& pFlat);
	void calc2Flat(const Eigen::VectorXf& mu, const Eigen::VectorXf& majEig, const Eigen::VectorXf& secEig, int ssId, FLOATVECTOR2& pFlat);

	void buildKDTree_for_Selection();
	// reduce duplicate data for 1-flats and 2-flats for better performance
	void removeDuplicateRecords(const std::vector<std::vector<float>>& rawData);
	// compute the strength of the normal vector
	float computeNormalVecStrength(const Eigen::Vector3f& vT1, const Eigen::Vector3f& vT2, const Eigen::Vector3f& vNorm, int method);

	// transform 1-flat indexed points to hyperbolic space
	void transformOneFlatIdxPts(std::vector<XPCPSample>& xpcpData);

	// Normalization of p-flats in subspaces
	void nomralizePFlatsInSubspaces(std::vector<XPCPSample>& xpcpData, int normMethod = 0); 
	// Compute the pdf for p-flats with a p = 1 or 2 based on the "mode" variable.
	// The inout variable "pdf" has "numBin" items. 
	// In "pdf", each FLOATVECTOR3 item records <x: min_val in bin, y: max_val in bin, z: sample count in bin>
	void computeSubspacePDF(std::vector<FLOATVECTOR3>& pdf, int subspaceId, int mode, int numBin, const std::vector<XPCPSample>& xpcpData); 
	// Compute CDF from a PDF
	void computeCDF(std::vector<FLOATVECTOR3>& cdf, const std::vector<FLOATVECTOR3>& pdf); 
	// Set the percentage of a sample in the whole population for each subspace
	void findSetSubspaceSamplePercentage(std::vector<XPCPSample>& xpcpData, int subspaceId, int mode, std::vector<FLOATVECTOR3>& cdf);



	//************************************************************
	// Functions for uncertain point data
	void uncertain_getPointCenters(const std::vector<GaussianXd>& distrData, std::vector<std::vector<float>>& pointCenterData);

	// Sample the neighborhood: for all the "num_nn" neighboring distributions, draw "num_samples_per_distr" samples each from a distribution
	void uncertain_sampleNeighborhood(const std::vector<GaussianXd>& distrData, int* nnIdx, int num_nn, std::vector<Eigen::MatrixXd>& pointNeighbors, int num_samples_per_distr);

	// draw "num_samples" samples from a single distribution
	void uncertain_drawSamplesPointDistrib(const GaussianXd& pointDistr, int num_samples, Eigen::MatrixXd& samples);

	// get all neighbors from the sampling with "sample_id" 
	void uncertain_get_neighbors_one_sample(const std::vector<Eigen::MatrixXd>& sampled_neighbors, int sample_id, std::vector<Eigen::VectorXd>& neighbors_one_sample);

	// get neighbors from all samples we get from the Monte-Carlo process
	void uncertain_get_neighbors_from_all_samples(const std::vector<Eigen::MatrixXd>& all_samples, std::vector<Eigen::VectorXd>& neighbors);

	// uncertain version of computing p-flats
	// Feb 24: need to think about the output data structure.
	void uncertain_calc_p_flats_subSpaceMethod(const Eigen::VectorXd& mu, const Eigen::VectorXd& majEig, const Eigen::VectorXd& secEig, UncertainXPCPSample& xpcp_tuple, int sampleId);

	//************************************************************

	static bool loadScaleFactors(const std::string& fileName, std::vector<double>& scaleLocs, std::vector<double>& scaleFactors);

private:
	int                       _dimRawData; // dimension of raw data
	int                       _dimXPCPdata; // dimension fo the XPCP data
	ANNkd_tree*               _annKdTree_rawData;
	PCPInselberg*             _pcp;
	// Record Min/Max for p-flat strength for each subspace
	std::vector<FLOATVECTOR2> _1flatsMinMaxPerSubspace;
	std::vector<FLOATVECTOR2> _2flatsMinMaxPerSubspace;
	UINT64                    _numSamples; // number of samples in the raw data
};

